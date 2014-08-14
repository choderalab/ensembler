def get_modeller_version():
    '''Hacky attempt to get Modeller version by regex searching the installation directory.
    '''
    import modeller
    import re
    regex = re.compile('/modeller-[0-9.]{2,6}/')
    match = re.search(regex, modeller.__file__)
    if match != None:
        version = match.group()[10:-1]
        return version
    else:
        return ''

def build_models(process_only_these_targets=None, process_only_these_templates=None, verbose=False):
    r'''Uses the build_model method to build homology models for a given set of
    targets and templates.

    MPI-enabled.
    '''
    import sys
    import os
    import Bio
    import Bio.SeqIO
    import mpi4py.MPI
    comm = mpi4py.MPI.COMM_WORLD 
    rank = comm.rank
    size = comm.size

    targets_dir = os.path.abspath('targets')
    templates_dir = os.path.abspath('templates')
    models_dir = os.path.abspath('models')
    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')

    targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )
    templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )
    ntemplates = len(templates)

    for target in targets:
        if process_only_these_targets and target.id not in process_only_these_targets: continue

        models_target_dir = os.path.join(models_dir, target.id)
        if rank == 0:
            print "========================================================================="
            print "Working on target '%s'" % (target.id)
            print "========================================================================="
            if not os.path.exists(models_target_dir):
                os.mkdir(models_target_dir)

        for template_index in range(rank, ntemplates, size):
            template = templates[template_index]
            if process_only_these_templates and template.id not in process_only_these_templates: continue
            build_model(target, template, verbose)

        comm.Barrier()

        # ========
        # Metadata
        # ========

        if rank == 0:
            import yaml
            import msmseeder
            import subprocess
            datestamp = msmseeder.core.get_utcnow_formatted()
            nsuccessful_models = subprocess.check_output(['find', models_target_dir, '-name', 'model.pdb']).count('\n')

            with open('meta.yaml') as meta_file:
                metadata = yaml.load(meta_file)
            with open('targets/meta.yaml') as meta_file:
                metadata.update(yaml.load(meta_file))
            with open('templates/meta.yaml') as meta_file:
                metadata.update(yaml.load(meta_file))

            metadata['build_models'] = {
                'target_id': target.id,
                'datestamp': datestamp,
                'nsuccessful_models': nsuccessful_models,
                'python_version': sys.version.split('|')[0].strip(),
                'python_full_version': sys.version,
                'msmseeder_version': msmseeder.__version__,
                'msmseeder_commit': msmseeder.core.get_src_git_commit_hash(),
                'modeller_version': get_modeller_version(),
                'biopython_version': Bio.__version__
            }

            metadata = msmseeder.core.ProjectMetadata(metadata)
            meta_filepath = os.path.join(models_target_dir, 'meta.yaml')
            metadata.write(meta_filepath)

    if rank == 0:

        print 'Done.'

def build_model(target, template, verbose=False):
    r'''Uses Modeller to build a homology model for a given target and
    template.

    Will not run Modeller if the output files already exist.

    Parameters
    ----------
    target : BioPython SeqRecord
    template : BioPython SeqRecord
        Must be a corresponding .pdb template file with the same ID in the
        templates/structures directory.
    '''
    # align target and template
    import os, tempfile, shutil, gzip
    import Bio.pairwise2
    import Bio.SubsMat.MatrixInfo

    templates_dir = os.path.abspath('templates')
    models_dir = os.path.abspath('models')
    models_target_dir = os.path.join(models_dir, target.id)
    model_dir = os.path.join(models_target_dir, template.id)
    aln_filename = os.path.join(model_dir, 'alignment.pir')
    seqid_filename = os.path.join(model_dir, 'sequence-identity.txt')
    model_pdbfilename = os.path.join(model_dir, 'model.pdb')
    restraint_filename_gz = os.path.join(model_dir, 'restraints.rsr.gz')
    current_dir = os.getcwd() 

    # Skip model-building if files already exist.
    files_to_check = [model_dir, model_pdbfilename, seqid_filename, aln_filename, restraint_filename_gz]
    files_are_present = [os.path.exists(filename) for filename in files_to_check]
    if all(files_are_present):
        if verbose: print "Output files already exist for target '%s' // template '%s'; files were not overwritten." % (target.id, template.id)
        return

    print "-------------------------------------------------------------------------"
    print "Modelling '%s' => '%s'" % (target.id, template.id)
    print "-------------------------------------------------------------------------"

    # Conduct alignment
    matrix = Bio.SubsMat.MatrixInfo.gonnet
    gap_open = -10
    gap_extend = -0.5
    aln = Bio.pairwise2.align.globalds(str(target.seq), str(template.seq), matrix, gap_open, gap_extend)

    # Create temp dir for modelling, and chdir
    temp_dir = tempfile.mkdtemp()
    try:
        os.chdir(temp_dir)

        # Write Modeller-format PIR alignment file
        tmp_aln_filename = 'aligned.pir'
        contents = "Target-template alignment by clustal omega\n"
        contents += ">P1;%s\n" % target.id
        contents += "sequence:%s:FIRST:@:LAST :@:::-1.00:-1.00\n" % target.id
        contents += aln[0][0] + '*\n'
        contents += ">P1;%s\n" % template.id
        contents += "structureX:%s:FIRST:@:LAST : :undefined:undefined:-1.00:-1.00\n" % template.id
        contents += aln[0][1] + '*\n'
        outfile = open('aligned.pir', 'w')
        outfile.write(contents)
        outfile.close()

        # Run Modeller
        import modeller
        import modeller.automodel
        modeller.log.none()
        env = modeller.environ()
        env.io.atom_files_directory = [os.path.join(templates_dir, 'structures')]

        a = modeller.automodel.allhmodel(env,
                                         # file with template codes and target sequence
                                         alnfile  = tmp_aln_filename,
                                         # PDB codes of the template
                                         knowns   = template.id,
                                         # code of the target
                                         sequence = target.id)
        a.make()                            # do homology modeling

        tmp_model_pdbfilename = a.outputs[0]['name']
        target_model = modeller.model(env, file=tmp_model_pdbfilename)

        # Create directory to place final models in
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

        target_model.write(file=model_pdbfilename)

        # Write sequence identity.
        with open(seqid_filename, 'w') as seqid_file:
            seqid_file.write('%.1f\n' % target_model.seq_id)

        # Copy alignment            
        shutil.move(tmp_aln_filename, aln_filename)

        # Copy restraints.
        with open('%s.rsr' % target.id, 'r') as rsrfile:
            with gzip.open(restraint_filename_gz, 'wb') as rsrgzfile:
                rsrgzfile.write(rsrfile.read())

        # Clean up temp dir
        os.chdir(current_dir)

        if os.path.getsize(model_pdbfilename) < 1:
            raise Exception, 'Output PDB file is empty. Could be a filesystem error.'

        text  = "---------------------------------------------------------------------------------\n"
        text += 'Successfully modeled target %s on template %s.\n' % (target.id, template.id)
        text += "Sequence identity was %.1f%%.\n" % (target_model.seq_id)
        return text

    finally:
        shutil.rmtree(temp_dir)

def sort_by_sequence_identity(process_only_these_targets=None, verbose=False):
    '''Compile sorted list of templates by sequence identity.
    Runs serially.
    '''
    import os
    import numpy
    import Bio.SeqIO
    import mpi4py.MPI
    comm = mpi4py.MPI.COMM_WORLD 
    rank = comm.rank

    if rank == 0:
        targets_dir = os.path.abspath("targets")
        templates_dir = os.path.abspath("templates")
        models_dir = os.path.abspath("models")

        targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
        targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )
        templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')
        templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )

        # ========
        # Compile sorted list by sequence identity
        # ========

        for target in targets:
            
            # Process only specified targets if directed.
            if process_only_these_targets and (target.id not in process_only_these_targets): continue

            models_target_dir = os.path.join(models_dir, target.id)
            if not os.path.exists(models_target_dir): continue

            print "-------------------------------------------------------------------------"
            print "Compiling template sequence identities for target %s" % (target.id)
            print "-------------------------------------------------------------------------"

            # ========
            # Build a list of valid models
            # ========

            if verbose: print "Building list of valid models..."
            valid_templates = list()
            for template in templates:
                model_filename = os.path.join(models_target_dir, template.id, 'model.pdb')
                if os.path.exists(model_filename):
                    valid_templates.append(template)

            nvalid = len(valid_templates)
            if verbose: print "%d valid models found" % nvalid

            # ========
            # Sort by sequence identity
            # ========

            if verbose: print "Sorting models in order of decreasing sequence identity..."
            seqids = numpy.zeros([nvalid], numpy.float32)
            for (template_index, template) in enumerate(valid_templates):
                model_seqid_filename = os.path.join(models_target_dir, template.id, 'sequence-identity.txt')
                with open(model_seqid_filename, 'r') as model_seqid_file:
                    firstline = model_seqid_file.readline().strip()
                seqid = float(firstline)
                seqids[template_index] = seqid
            sorted_seqids = numpy.argsort(-seqids)

            # ========
            # Write templates sorted by sequence identity
            # ========

            seq_ofilename = os.path.join(models_target_dir, 'sequence-identities.txt')
            with open(seq_ofilename, 'w') as seq_ofile:
                for index in sorted_seqids:
                    template = valid_templates[index]
                    identity = seqids[index]
                    seq_ofile.write('%-40s %6.1f\n' % (template.id, identity))

            # ========
            # Metadata
            # ========

            import sys
            import yaml
            import msmseeder
            datestamp = msmseeder.core.get_utcnow_formatted()

            meta_filepath = os.path.join(models_target_dir, 'meta.yaml')
            with open(meta_filepath) as meta_file:
                metadata = yaml.load(meta_file)

            metadata['sort_by_sequence_identity'] = {
                'target_id': target.id,
                'datestamp': datestamp,
                'python_version': sys.version.split('|')[0].strip(),
                'python_full_version': sys.version,
                'msmseeder_version': msmseeder.__version__,
                'msmseeder_commit': msmseeder.core.get_src_git_commit_hash(),
                'biopython_version': Bio.__version__
            }

            metadata = msmseeder.core.ProjectMetadata(metadata)
            meta_filepath = os.path.join(models_target_dir, 'meta.yaml')
            metadata.write(meta_filepath)

    comm.Barrier()
    if rank == 0:
        print 'Done.'

def cluster_models(process_only_these_targets=None, verbose=False):
    '''Cluster models based on RMSD, and filter out non-unique models as
    determined by a given cutoff.

    Runs serially.
    '''
    import os, glob
    import Bio.SeqIO
    import mdtraj
    import mpi4py.MPI
    comm = mpi4py.MPI.COMM_WORLD 
    rank = comm.rank

    if rank == 0:
        targets_dir = os.path.abspath("targets")
        templates_dir = os.path.abspath("templates")
        models_dir = os.path.abspath("models")

        targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
        targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )
        templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')
        templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )

        cutoff = 0.06 # Cutoff for RMSD clustering (nm)

        for target in targets:
            if process_only_these_targets and (target.id not in process_only_these_targets): continue

            models_target_dir = os.path.join(models_dir, target.id)
            if not os.path.exists(models_target_dir): continue

            # =============================
            # Construct a mdtraj trajectory containing all models
            # =============================

            print 'Building a list of valid models...'

            model_pdbfilenames = []
            valid_templateIDs = []
            for t, template in enumerate(templates):
                model_dir = os.path.join(models_target_dir, template.id)
                model_pdbfilename = os.path.join(model_dir, 'model.pdb')
                if not os.path.exists(model_pdbfilename):
                    continue
                model_pdbfilenames.append(model_pdbfilename)
                valid_templateIDs.append(template.id)

            print 'Constructing a trajectory containing all valid models...'

            traj = mdtraj.load(model_pdbfilenames)

            # =============================
            # Clustering
            # =============================

            print 'Conducting RMSD-based clustering...'

            # Remove any existing unique_by_clustering files
            for f in glob.glob( models_target_dir+'/*_PK_*/unique_by_clustering' ):
                os.unlink(f)

            # Each template will be added to the list uniques if it is further than
            # 0.2 Angstroms (RMSD) from the nearest template.
            uniques=[]
            min_rmsd = []
            for (t, templateID) in enumerate(valid_templateIDs):
                model_dir = os.path.join(models_target_dir, templateID)

                # Add the first template to the list of uniques
                if t==0:
                    uniques.append(templateID)
                    with open( os.path.join(model_dir, 'unique_by_clustering'), 'w') as unique_file: pass
                    continue

                # Cluster using CA atoms
                CAatoms = [a.index for a in traj.topology.atoms if a.name == 'CA']
                rmsds = mdtraj.rmsd(traj[0:t], traj[t], atom_indices=CAatoms, parallel=False)
                min_rmsd.append( min(rmsds) )

                if min_rmsd[-1] < cutoff:
                    continue
                else:
                    uniques.append( templateID )
                    # Create a blank file to say this template was found to be unique
                    # by clustering
                    with open( os.path.join(model_dir, 'unique_by_clustering'), 'w') as unique_file: pass

            with open( os.path.join(models_target_dir, 'unique-models.txt'), 'w') as uniques_file:
                for u in uniques:
                    uniques_file.write(u+'\n')
                print '%d unique models (from original set of %d) using cutoff of %.3f nm' % (len(uniques), len(valid_templateIDs), cutoff)

            # ========
            # Metadata
            # ========

            import sys
            import yaml
            import msmseeder
            import mdtraj.version
            datestamp = msmseeder.core.get_utcnow_formatted()

            meta_filepath = os.path.join(models_target_dir, 'meta.yaml')
            with open(meta_filepath) as meta_file:
                metadata = yaml.load(meta_file)

            metadata['cluster_models'] = {
                'target_id': target.id,
                'datestamp': datestamp,
                'python_version': sys.version.split('|')[0].strip(),
                'python_full_version': sys.version,
                'msmseeder_version': msmseeder.__version__,
                'msmseeder_commit': msmseeder.core.get_src_git_commit_hash(),
                'biopython_version': Bio.__version__,
                'mdtraj_version': mdtraj.version.short_version,
                'mdtraj_commit': mdtraj.version.git_revision
            }

            metadata = msmseeder.core.ProjectMetadata(metadata)
            meta_filepath = os.path.join(models_target_dir, 'meta.yaml')
            metadata.write(meta_filepath)

    comm.Barrier()
    if rank == 0:
        print 'Done.'

