import os
import shutil
import ensembler.cli_commands

helpstring = """\
Conduct a testrun of the entire ensembler pipeline in the current directory,
using the target protein EGFR_HUMAN, and one template structure from AURKB_HUMAN.
"""
docopt_helpstring = ''


def dispatch(args):

    ensembler.cli_commands.init.dispatch(
    args={
    '--project_dir': '.',
    }
    )

    ensembler.cli_commands.gather_targets.dispatch(
    args={
    '--gather_from': 'uniprot',
        '--query': 'mnemonic:EGFR_HUMAN',
        '--uniprot_domain_regex': '^Protein kinase',
    }
    )

    structure_paths = ','.join([
        os.path.join(ensembler.core.src_toplevel_dir, 'tests', 'integration_test_resources', 'structures', 'pdb'),
        os.path.join(ensembler.core.src_toplevel_dir, 'tests', 'integration_test_resources', 'structures', 'sifts'),
    ])

    ensembler.cli_commands.gather_templates.dispatch(
        args={
            '--gather_from': 'uniprot',
            '--query': 'mnemonic:AURKB_HUMAN',
            '--uniprot_domain_regex': '^Protein kinase',
            '--structure_paths': structure_paths,
            '--no-loopmodel': True,
            '--overwrite_structures': False,
        }
    )

    ensembler.cli_commands.align.dispatch(
        args={
            '--targets': False,
            '--templates': 'AURKB_HUMAN_D0_4AF3_A',
            '--verbose': False,
        }
    )

    ensembler.cli_commands.build_models.dispatch(
        args={
            '--targets': False,
            '--templates': 'AURKB_HUMAN_D0_4AF3_A',
            '--verbose': False,
        }
    )

    ensembler.cli_commands.cluster.dispatch(
        args={
            '--targets': False,
            '--verbose': False,
        }
    )

    ensembler.cli_commands.refine_implicit.dispatch(
        args={
            '--targetsfile': False,
            '--targets': False,
            '--templates': 'AURKB_HUMAN_D0_4AF3_A',
            '--verbose': False,
            '--openmm_platform': False,
            '--gpupn': False,
            '--simlength': '1.0',
        }
    )

    ensembler.cli_commands.solvate.dispatch(
        args={
            '--targetsfile': False,
            '--targets': False,
            '--templates': 'AURKB_HUMAN_D0_4AF3_A',
            '--verbose': False,
        }
    )

    ensembler.cli_commands.refine_explicit.dispatch(
        args={
            '--targetsfile': False,
            '--targets': False,
            '--templates': 'AURKB_HUMAN_D0_4AF3_A',
            '--verbose': False,
            '--openmm_platform': False,
            '--gpupn': False,
            '--simlength': '1.0',
        }
    )

    ensembler.cli_commands.package_models.dispatch(
        args={
            '--targetsfile': False,
            '--targets': False,
            '--templates': 'AURKB_HUMAN_D0_4AF3_A',
            '--nFAHclones': False,
            '--archiveFAHproject': False,
        }
    )