def build_models(process_only_these_targets=None, process_only_these_templates=None):
    import os
    import Bio.SeqIO
    import mpi4py
    from mpi4py import MPI
    comm = MPI.COMM_WORLD 
    rank = comm.rank
    size = comm.size

    targets_dir = os.path.abspath('targets')
    templates_dir = os.path.abspath('templates')
    models_dir = os.path.abspath('models')
    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')

    targets = Bio.SeqIO.parse(targets_fasta_filename, 'fasta')
    templates = list( Bio.SeqIO.parse(templates_fasta_filename, 'fasta') )
    ntemplates = len(templates)

    for target in targets:
        if process_only_these_targets and target.id not in process_only_these_targets: continue

        if rank == 0:
            print "========================================================================="
            print "Modelling '%s'" % (target.id)
            print "========================================================================="

        models_target_dir = os.path.join(models_dir, target.id)
        if not os.path.exists(models_target_dir):
            os.mkdir(models_target_dir)

        for template_index in range(rank, ntemplates, size):
            template = templates[template_index]
            if process_only_these_templates and template.id not in process_only_these_templates: continue
            build_model(target, template)

def build_model(target, template):
    # align target and template
    import os, tempfile, shutil, gzip
    import Bio.pairwise2
    import Bio.SubsMat.MatrixInfo

    print "-------------------------------------------------------------------------"
    print "Modelling '%s' => '%s'" % (target.id, template.id)
    print "-------------------------------------------------------------------------"

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
        text  = "---------------------------------------------------------------------------------\n"
        text += "Output files already exist for target %s, template %s; files were not overwritten.\n" % (target.id, template.id)
        print text
        return

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

