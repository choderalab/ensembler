import os
import yaml
import mdtraj
from ensembler.refinement import refine_implicit_md, refine_explicit_md
from simtk import unit
from ensembler.core import default_project_dirnames
from ensembler.tests.integrationtest_utils import integrationtest_context
from nose.plugins.attrib import attr


@attr('unit')
def test_refine_implicit_md_short():
    with integrationtest_context(set_up_project_stage='clustered'):
        #targetid = 'EGFR_HUMAN_D0'
        #templateid = 'KC1D_HUMAN_D0_4KB8_D'
        targetid = 'MYB_HUMAN_D0'
        templateid = 'MYB_MOUSE_D0_1GUU'
        refine_implicit_md(
            process_only_these_targets=[targetid],
            process_only_these_templates=[templateid],
            sim_length=2.0*unit.femtosecond,
            nsteps_per_iteration=1,
            minimization_steps=1,
            loglevel='debug'
        )
        implicit_metadata_filepath = os.path.join(
            default_project_dirnames.models, targetid, 'refine_implicit_md-meta0.yaml'
        )
        implicit_model_filepath = os.path.join(
            default_project_dirnames.models, targetid, templateid, 'implicit-refined.pdb.gz'
        )
        implicit_energies_filepath = os.path.join(
            default_project_dirnames.models, targetid, templateid, 'implicit-energies.txt'
        )
        implicit_log_filepath = os.path.join(
            default_project_dirnames.models, targetid, templateid, 'implicit-log.yaml'
        )

        assert all(map(
            os.path.exists,
            [implicit_model_filepath, implicit_energies_filepath, implicit_log_filepath]
        ))
        with open(implicit_log_filepath) as implicit_log_file:
            implicit_log = yaml.load(implicit_log_file)
        assert implicit_log.get('finished') is True
        assert implicit_log.get('successful') is True
        assert implicit_log.get('ph') == 8.0
        assert os.path.exists(implicit_metadata_filepath)
        with open(implicit_metadata_filepath) as implicit_metadata_file:
            implicit_metadata = yaml.load(implicit_metadata_file)
        #assert implicit_metadata.get('refine_implicit_md').get('custom_residue_variants') == {
        #    'EGFR_HUMAN_D0': {49: 'ASH'}
        assert implicit_metadata.get('refine_implicit_md').get('custom_residue_variants') == {
            'MYB_HUMAN_D0': {47: 'ASH'}
        }
        implicit_model_traj = mdtraj.load_pdb(implicit_model_filepath)
        resis = [resi for resi in implicit_model_traj.top.residues]
        #resi49 = resis[49]
        #resi49_atom_strings = [str(atom) for atom in resi49.atoms]
        #assert 'ASP50-HD2' in resi49_atom_strings
        resi47 = resis[47]
        resi47_atom_strings = [str(atom) for atom in resi47.atoms]
        assert 'ASP48-HD2' in resi47_atom_strings


@attr('slow')
def test_refine_explicit_md_short():
    with integrationtest_context(set_up_project_stage='solvated'):
        targetid = 'EGFR_HUMAN_D0'
        templateid = 'KC1D_HUMAN_D0_4KB8_D'
        refine_explicit_md(
            process_only_these_targets=[targetid],
            process_only_these_templates=[templateid],
            sim_length=2.0*unit.femtosecond,
            nsteps_per_iteration=1,
            verbose=True
        )
        explicit_metadata_filepath = os.path.join(
            default_project_dirnames.models, targetid, 'refine_explicit_md-meta0.yaml'
        )
        explicit_model_filepath = os.path.join(
            default_project_dirnames.models, targetid, templateid, 'explicit-refined.pdb.gz'
        )
        explicit_energies_filepath = os.path.join(
            default_project_dirnames.models, targetid, templateid, 'explicit-energies.txt'
        )
        explicit_log_filepath = os.path.join(
            default_project_dirnames.models, targetid, templateid, 'explicit-log.yaml'
        )

        assert all(map(
            os.path.exists,
            [explicit_model_filepath, explicit_energies_filepath, explicit_log_filepath]
        ))
        with open(explicit_log_filepath) as explicit_log_file:
            explicit_log = yaml.load(explicit_log_file)
        assert explicit_log.get('finished') is True
        assert explicit_log.get('successful') is True
        explicit_model_traj = mdtraj.load_pdb(explicit_model_filepath)
