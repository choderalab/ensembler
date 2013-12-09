# Take model refined in OpenMM, make Gromacs topology, and run a brief energy minimzation and MD equilibration
#
# Daniel L. Parton <partond@mskcc.org> - 13 Mar 2013
#
# PREREQUISITES
#
# TODO Change to amber99sb-ildn-star?
#
# =============================
# IMPORTS
# =============================

import sys,re,os,collections
from subprocess import call
import choderalab as clab

# =============================
# PARAMETERS
# =============================

#process_only_these_targets = ['EGFR_HUMAN_P00533_PK0']
#process_only_these_templates = ['EGFR_HUMAN_P00533_PK0_1M14_A'] # XXX tmp

process_only_these_targets = ['BRAF_HUMAN_P15056_PK0', 'MET_HUMAN_P08581_PK0', 'ERBB2_HUMAN_P04626_PK0', 'AKT1_HUMAN_P31749_PK0', 'JAK2_HUMAN_O60674_PK0', 'TITIN_HUMAN_Q8WZ42_PK0', 'ALK_HUMAN_Q9UM73_PK0', 'FGFR2_HUMAN_P21802_PK0', 'RET_HUMAN_P07949_PK0', 'CDK2_HUMAN_P24941_PK0', 'KIT_HUMAN_P10721_PK0', 'MK14_HUMAN_Q16539_PK0']
process_only_these_templates = ['BRAF_HUMAN_P15056_PK0_4EHG_A', 'MET_HUMAN_P08581_PK0_3I5N_A', 'ERBB2_HUMAN_P04626_PK0_3PP0_A', 'AKT1_HUMAN_P31749_PK0_4EJN_A', 'JAK2_HUMAN_O60674_PK0_2B7A_A', 'TITIN_HUMAN_Q8WZ42_PK0_1TKI_A', 'ALK_HUMAN_Q9UM73_PK0_2YFX_A', 'FGFR2_HUMAN_P21802_PK0_2PZ5_A', 'RET_HUMAN_P07949_PK0_2IVT_A', 'CDK2_HUMAN_P24941_PK0_3EID_A', 'CSF1R_HUMAN_P07333_PK0_2I0Y_A', 'MK14_HUMAN_Q16539_PK0_3C5U_A']

ff = 'amber99sb-ildn'
water_model = 'tip3p'

if '-gmx_multiprocessing' in sys.argv:
    gmx_multiprocessing = True
else:
    gmx_multiprocessing = False

# =============================
# PATHS
# =============================

targets_directory = os.path.abspath("targets") # target sequences for modeling
templates_directory = os.path.abspath("templates") # template structures for use in modeling
models_directory = os.path.abspath("models")

# Gromacs bin
gmx_bin = '/usr/local/gromacs/4.5.5/bin'
pdb2gmx_path = os.path.join(gmx_bin, 'pdb2gmx')
editconf_path = os.path.join(gmx_bin, 'editconf')
make_ndx_path = os.path.join(gmx_bin, 'make_ndx')
grompp_path = os.path.join(gmx_bin, 'grompp')
mdrun_path = os.path.join(gmx_bin, 'mdrun') # NB: not using mpi version of mdrun

# =============================
# READ TEMPLATE AND TARGET INDICES
# =============================

targets_index_filename = os.path.join(targets_directory, 'targets.txt')
with open(targets_index_filename, 'r') as infile:
    targets = [ line.strip() for line in infile ]
print '%d target sequences' % len(targets)

templates_index_filename = os.path.join(templates_directory, 'templates.txt')
with open(templates_index_filename, 'r') as infile:
    templates = [ line.strip() for line in infile ]
print '%d template structures' % len(templates)

# =============================
# SET UP INPUT DATA (TPR) FOR EACH TARGET AND TEMPLATE
# =============================

for target in targets:

    # Process only specified targets if directed.
    if process_only_these_targets and (target not in process_only_these_targets): continue

    target_directory = os.path.join(models_directory, target)
    if not os.path.exists(target_directory): continue

    for template in templates:
        # Process only specified templates if directed.
        if process_only_these_templates and (template not in process_only_these_templates): continue

        model_directory = os.path.join(models_directory, target, template)
        if not os.path.exists(model_directory):
            continue
        # Create subdirectory to be used for gromacs simulations
        gmx_directory = os.path.join(model_directory, 'gromacs')
        if not os.path.exists(gmx_directory):
            os.makedirs(gmx_directory)

        stdout_path = os.path.join(gmx_directory, 'tmp-stdout')
        stderr_path = os.path.join(gmx_directory, 'tmp-stderr')
        stdin_path = os.path.join(gmx_directory, 'tmp-stdin')

        # Now check if OpenMM explicit refinement has completed correctly.
        # Note that failures are expected for some models, in which case ignore.
        omm_refined_path = os.path.join(model_directory, 'explicit-refined.pdb')
        if not os.path.exists(omm_refined_path):
            continue

        tpr_prod_path = os.path.join(gmx_directory, 'prod.tpr')
        if not os.path.exists(tpr_prod_path):
            # ============
            # Create topology file
            # (Need to do this for each model, since histidine titration state may differ)
            # ============
            print '------------------------------------------------'
            print 'Creating topology file for target %s...' % target

            # Various files which will be created and later deleted
            omm_protein_path = os.path.join(gmx_directory, 'tmp-omm-protein.pdb')
            omm_nonprotein_path = os.path.join(gmx_directory, 'tmp-omm-nonprotein.pdb')
            gro_nonprotein_path = os.path.join(gmx_directory, 'tmp-nonprotein.gro')
            gro_post_pdb2gmx_path = os.path.join(gmx_directory, 'tmp-post-pdb2gmx.gro')
            top_post_pdb2gmx_path = os.path.join(gmx_directory, 'tmp-post-pdb2gmx.top')
            posre_path = os.path.join(gmx_directory, 'tmp-posre.itp')
            gro_protein_path = os.path.join(gmx_directory, 'tmp-protein.gro')
            gro_system_path = os.path.join(gmx_directory, 'tmp-system.gro')
            top_path = os.path.join(gmx_directory, 'tmp-topol.top')
            index_path = os.path.join(gmx_directory, 'tmp-index.ndx')
            em_mdp_path = os.path.join(gmx_directory, 'tmp-em.mdp')
            tpr_em_path = os.path.join(gmx_directory, 'tmp-em.tpr')
            mdout_path = os.path.join(gmx_directory, 'tmp-mdout.mdp')
            trr_em_path = os.path.join(gmx_directory, 'tmp-em.trr')
            gro_em_path = os.path.join(gmx_directory, 'tmp-em.gro')
            energy_em_path = os.path.join(gmx_directory, 'tmp-em.edr')
            log_em_path = os.path.join(gmx_directory, 'tmp-em.log')
            eq_mdp_path = os.path.join(gmx_directory, 'tmp-eq.mdp')
            tpr_eq_path = os.path.join(gmx_directory, 'tmp-eq.tpr')
            mdout_path = os.path.join(gmx_directory, 'tmp-mdout.mdp')
            trr_eq_path = os.path.join(gmx_directory, 'tmp-eq.trr')
            gro_eq_path = os.path.join(gmx_directory, 'tmp-eq.gro')
            energy_eq_path = os.path.join(gmx_directory, 'eq.edr')
            log_eq_path = os.path.join(gmx_directory, 'tmp-eq.log')
            xtc_eq_path = os.path.join(gmx_directory, 'tmp-eq.xtc')
            snap_eq_path = os.path.join(gmx_directory, 'tmp-eq.cpt')
            prod_mdp_path = os.path.join(gmx_directory, 'tmp-prod.mdp')

            # Take the input file (from OpenMM) and split the protein and
            # non-protein parts into separate .pdb files.
            # For the protein part, translate the hydrogen-naming from OpenMM
            # style to Gromacs style.
            with open(omm_refined_path, 'r') as omm_refined:

                # First get histidine assignments, as the resnames must be changed to
                # work with pdb2gmx. This method can account for PDB wrapping.
                # This is a fast method for grepping HIS lines into a variable.
                his_lines = '\t'.join(line for line in omm_refined if ('HD1 HIS' in line or 'HE2 HIS' in line))
                prev_resid = ''
                h_res_iter = -1
                his_data = list()
                for h in his_lines.split('\t'):
                    atomname = h[11:16]
                    resid = h[22:27]
                    if resid != prev_resid:
                        h_res_iter+=1
                        his_data.append({'resid':resid,'atomnames':list()})
                    his_data[h_res_iter]['atomnames'].append(atomname)
                    prev_resid = h[22:27]

                his_assignments = list()
                for h in range(len(his_data)):
                    resid = his_data[h]['resid']
                    if len(his_data[h]['atomnames']) == 2:
                        resname = 'HIP'
                    elif '  HD1' in his_data[h]['atomnames']:
                        resname = 'HID'
                    elif '  HE2' in his_data[h]['atomnames']:
                        resname = 'HIE'
                    his_assignments.append( {'resid':resid , 'resname':resname})

                print 'Assigning histidines based on atoms present, as follows:'
                print his_assignments
                omm_refined.seek(0)

                prev_his_resid = 0
                h_res_iter = -1
                nonprotein_counts = collections.OrderedDict() # needs to be an ordered dict, so that the nonprotein molecules will be written out in the topol file in the correct order
                # This is the loop where we write out the protein and nonprotein parts
                with open(omm_protein_path, 'w') as omm_protein:
                    with open(omm_nonprotein_path, 'w') as omm_nonprotein:
                        for line in omm_refined.readlines():
                            try:
                                if line[0:5] == 'CRYST':
                                    omm_nonprotein.write(line)
                                elif line[0:4] == 'ATOM':
                                    if line[17:20] in clab.gmx.aa3codes:
                                        # For histidines, must rename to HID, HIE or HIP
                                        # This method can account for resid wrapping in the PDB file
                                        if line[17:20] == 'HIS':
                                            resid = line[22:27]
                                            if resid != prev_his_resid:
                                                h_res_iter += 1
                                            line = line.replace('HIS', his_assignments[h_res_iter]['resname'])
                                            prev_his_resid = resid
                                        # Get updated hydrogen names
                                        Hmap_line = clab.gmx.Hmap_pdb_omm2gmx(line)
                                        # Write to protein file
                                        omm_protein.write(Hmap_line)
                                        continue
                                    # Sort out ion naming
                                    elif line[12:20] == 'Na    NA':
                                        line = line.replace('Na    NA','NA    NA')
                                    elif line[12:20] == 'Cl    CL':
                                        line = line.replace('Cl    CL','CL    CL')
                                    # Take counts of nonprotein atoms in dict. resnames are the keys
                                    try:
                                        nonprotein_counts[line[17:20]] += 1
                                    except KeyError:
                                        nonprotein_counts[line[17:20]] = 1
                                    # Write to nonprotein file
                                    omm_nonprotein.write(line)
                            except IndexError:
                                continue

            print nonprotein_counts
            nwaters = nonprotein_counts['HOH'] / 3

            # Now run pdb2gmx on the protein .pdb file
            stdout_file = open(stdout_path, 'w')
            stderr_file = open(stderr_path, 'w')
            # pdb2gmx
            print 'Running pdb2gmx...'
            retcode = call([pdb2gmx_path, '-f', omm_protein_path, '-o', gro_post_pdb2gmx_path, '-p', top_post_pdb2gmx_path, '-i', posre_path, '-ff', ff, '-water', water_model], stdout=stdout_file, stderr=stderr_file)
            stdout_file.close()
            stderr_file.close()
            if retcode != 0:
                # Clean up any # files
                for f in os.listdir('.'):
                    if re.search('#tmp*',f):
                        os.remove(f)
                raise Exception('pdb2gmx failed while working on:\ntarget: %s\ntemplate: %s\n' % (target, template))

            # Read in the .gro file which was output by pdb2gmx, and convert hydrogen names back to OpenMM style
            with open(gro_post_pdb2gmx_path, 'r') as gro_post_pdb2gmx_file:
                with open(gro_protein_path, 'w') as gro_protein_file:
                    for line in gro_post_pdb2gmx_file.readlines():
                        new_line = re.sub('HI[DEP]','HIS', clab.gmx.Hmap_gro_gmx2omm(line))
                        gro_protein_file.write(new_line)

            # Read in the .top file which was output from pdb2gmx, and convert hydrogen names to OpenMM style
            with open(top_post_pdb2gmx_path, 'r') as ifile:
                lines = ifile.readlines()
            # Write out again with altered naming
            with open(top_path, 'w') as top_file:
                for line in lines:
                    # Add HOH params (Gromacs expects water to be named SOL)
                    if line == '[ moleculetype ]\n':
                        top_file.write(clab.gmx.top_hoh_text)
                    new_line = re.sub('HI[DEP]','HIS', clab.gmx.Hmap_top_gmx2omm(line))
                    top_file.write(new_line)
                for n in nonprotein_counts.keys():
                    nmols = nonprotein_counts[n]
                    if n == 'HOH':
                        nmols /= 3
                    top_file.write('%s    %d\n' % (n.strip(),nmols))

            # Convert the nonprotein.pdb file to .gro format, and then put together with the protein .gro
            with open(stdout_path, 'w') as stdout_file:
                with open(stderr_path, 'w') as stderr_file:
                    retcode = call([editconf_path, '-f', omm_nonprotein_path, '-o', gro_nonprotein_path], stdout=stdout_file, stderr=stderr_file)
            with open(gro_system_path, 'w') as gro_system_file:
                with open(gro_protein_path, 'r') as gro_protein_file:
                    gro_protein = gro_protein_file.readlines()
                    gro_system_file.write(gro_protein[0]) # Header
                    with open(gro_nonprotein_path, 'r') as gro_nonprotein_file:
                        gro_nonprotein = gro_nonprotein_file.readlines()
                        natoms = len(gro_protein)-3 + len(gro_nonprotein)-3
                        gro_system_file.write('%5d\n' % natoms) # No. atoms
                        for a in range(2,len(gro_protein)-1):
                            gro_system_file.write(gro_protein[a]) # Protein coords
                        for a in range(2,len(gro_nonprotein)-1):
                            gro_system_file.write(gro_nonprotein[a]) # Nonprotein coords
                        gro_system_file.write(gro_nonprotein[-1]) # Box

            # editconf to reset resid numbering
            with open(stdout_path, 'w') as stdout_file:
                with open(stderr_path, 'w') as stderr_file:
                    retcode = call([editconf_path, '-f', gro_system_path, '-o', gro_system_path, '-resnr', '1'], stdout=stdout_file, stderr=stderr_file)

            # ============
            # Set up and run energy minimization
            # ============
            # Make .ndx file
            with open(stdout_path, 'w') as stdout_file:
                with open(stderr_path, 'w') as stderr_file:
                    with open(stdin_path, 'w') as stdin_file:
                        stdin_file.write('q\n')
                    with open(stdin_path, 'r') as stdin_file:
                        # make_ndx
                        print 'Running make_ndx...'
                        retcode = call([make_ndx_path, '-f', gro_system_path, '-o', index_path], stdin=stdin_file, stdout=stdout_file, stderr=stderr_file)
                        if retcode != 0:
                            # Clean up any # files
                            for f in os.listdir(gmx_directory):
                                if re.search('#',f):
                                    os.remove(os.path.join(gmx_directory,f))
                            raise Exception('make_ndx failed while working on:\ntarget: %s\ntemplate: %s\n' % (target, template))

            # Make tmp em.mdp file
            with open(em_mdp_path, 'w') as em_mdp_file:
                em_mdp_file.write(clab.gmx.em_mdp_contents)

            # grompp em.tpr file
            if not os.path.exists(tpr_em_path):
                print 'grompp: preparing binary input file for energy minimzation...'
                with open(stdout_path, 'w') as stdout_file:
                    with open(stderr_path, 'w') as stderr_file:
                        retcode = call([grompp_path, '-f', em_mdp_path, '-c', gro_system_path, '-n', index_path, '-p', top_path, '-o', tpr_em_path, '-po', mdout_path], stdout=stdout_file, stderr=stderr_file)
                        if retcode != 0:
                            # Clean up any # files
                            for f in os.listdir(gmx_directory):
                                if re.search('#',f):
                                    os.remove(os.path.join(gmx_directory,f))
                            raise Exception('grompp failed while working on:\ntarget: %s\ntemplate: %s\n' % (target, template))

            # Run energy minimization
            if not os.path.exists(gro_em_path):
                print 'mdrun: Running energy minimzation...'
                with open(stdout_path, 'w') as stdout_file:
                    with open(stderr_path, 'w') as stderr_file:
                        gmx_call_args = [mdrun_path, '-s', tpr_em_path, '-v', '-o', trr_em_path, '-c', gro_em_path, '-e', energy_em_path, '-g', log_em_path]
                        if gmx_multiprocessing:
                            retcode = call(gmx_call_args, stdout=stdout_file, stderr=stderr_file)
                        else:
                            gmx_call_args += ['-nt', '1']
                            retcode = call(gmx_call_args, stdout=stdout_file, stderr=stderr_file)
                        if retcode != 0:
                            # Clean up any # files
                            for f in os.listdir(gmx_directory):
                                if re.search('#',f):
                                    os.remove(os.path.join(gmx_directory,f))
                            raise Exception('mdrun failed while working on:\ntarget: %s\ntemplate: %s\n' % (target, template))

            # ============
            # Set up and run MD equilibration
            # ============
            # Make tmp eq.mdp file
            with open(eq_mdp_path, 'w') as eq_mdp_file:
                eq_mdp_file.write(clab.gmx.eq_mdp_contents)

            # grompp eq.tpr file
            if not os.path.exists(tpr_eq_path):
                print 'grompp: preparing binary input file for MD equilibration...'
                with open(stdout_path, 'w') as stdout_file:
                    with open(stderr_path, 'w') as stderr_file:
                        retcode = call([grompp_path, '-f', eq_mdp_path, '-c', gro_em_path, '-n', index_path, '-p', top_path, '-o', tpr_eq_path, '-po', mdout_path], stdout=stdout_file, stderr=stderr_file)
                        if retcode != 0:
                            # Clean up any # files
                            for f in os.listdir(gmx_directory):
                                if re.search('#',f):
                                    os.remove(os.path.join(gmx_directory,f))
                            raise Exception('grompp failed while working on:\ntarget: %s\ntemplate: %s\n' % (target, template))

            # Run equilibration (1 core, 1 thread)
            if not os.path.exists(gro_eq_path):
                print 'mdrun: running equilibration...'
                with open(stdout_path, 'w') as stdout_file:
                    with open(stderr_path, 'w') as stderr_file:
                        gmx_call_args = [mdrun_path, '-s', tpr_eq_path, '-v', '-o', trr_eq_path, '-x', xtc_eq_path, '-c', gro_eq_path, '-e', energy_eq_path, '-g', log_eq_path, '-cpo', snap_eq_path]
                        if gmx_multiprocessing:
                            retcode = call(gmx_call_args, stdout=stdout_file, stderr=stderr_file)
                        else:
                            gmx_call_args += ['-nt', '1']
                            retcode = call(gmx_call_args, stdout=stdout_file, stderr=stderr_file)
                        if retcode != 0:
                            # Clean up any # files
                            for f in os.listdir(gmx_directory):
                                if re.search('#',f):
                                    os.remove(os.path.join(gmx_directory,f))
                            raise Exception('mdrun failed while working on:\ntarget: %s\ntemplate: %s\n' % (target, template))

            # ============
            # Produce the input file for a MD production run
            # ============
            # Make tmp prod.mdp file
            with open(prod_mdp_path, 'w') as prod_mdp_file:
                prod_mdp_file.write(clab.gmx.prod_mdp_contents)

            # grompp prod.tpr file
            if not os.path.exists(tpr_prod_path):
                print 'grompp: preparing binary input file for MD production run...'
                with open(stdout_path, 'w') as stdout_file:
                    with open(stderr_path, 'w') as stderr_file:
                        retcode = call([grompp_path, '-f', prod_mdp_path, '-c', gro_eq_path, '-n', index_path, '-p', top_path, '-o', tpr_prod_path, '-po', mdout_path], stdout=stdout_file, stderr=stderr_file)
                        if retcode != 0:
                            # Clean up any # files
                            for f in os.listdir(gmx_directory):
                                if re.search('#',f):
                                    os.remove(os.path.join(gmx_directory,f))
                            raise Exception('grompp failed while working on:\ntarget: %s\ntemplate: %s\n' % (target, template))

            # TODO Save ns/day from std_err file

            # Clean up tmp files and any remaining # files
            for f in os.listdir(gmx_directory):
                if re.search('tmp-',f):
                    os.remove(os.path.join(gmx_directory,f))
                elif re.search('^#',f):
                    os.remove(os.path.join(gmx_directory,f))
            print 'Done.'

