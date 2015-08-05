import os
from nose.plugins.attrib import attr
from ensembler.packaging import package_for_fah
from ensembler.core import default_project_dirnames
from ensembler.tests.integrationtest_utils import integrationtest_context


@attr('unit')
def test_package_for_fah():
    with integrationtest_context(set_up_project_stage='refined_explicit'):
        package_for_fah(
            process_only_these_targets=['EGFR_HUMAN_D0'],
            process_only_these_templates=[
                'KC1D_HUMAN_D0_4HNF_A',
                'KC1D_HUMAN_D0_4KB8_D'
            ]
        )
        packaged_project_base_path = os.path.join(
            default_project_dirnames.packaged_models,
            'fah-projects',
            'EGFR_HUMAN_D0'
        )
        assert os.path.exists(packaged_project_base_path)
        assert os.path.exists(os.path.join(
            packaged_project_base_path,
            'RUN0'
        ))
        assert os.path.exists(os.path.join(
            packaged_project_base_path,
            'RUN1'
        ))
        target_filenames = [
            'system.xml',
            'integrator.xml',
        ]

        for target_filename in target_filenames:
            assert os.path.exists(os.path.join(
                packaged_project_base_path,
                target_filename
            ))

        run_filenames = [
            'template.txt',
            'system.pdb',
            'protein.pdb',
            'sequence-identity.txt',
            'state0.xml',
        ]

        for run_id in range(2):
            for run_filename in run_filenames:
                assert os.path.exists(os.path.join(
                    packaged_project_base_path,
                    'RUN{}'.format(run_id),
                    run_filename
                ))
