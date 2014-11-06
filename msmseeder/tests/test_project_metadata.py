import os
import yaml
import msmseeder


def test_project_metadata_init():
    project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
    test_data = {'test_field': 'test_value'}
    project_metadata.add_data(test_data)
    assert project_metadata.data == {'init': test_data}


def test_project_metadata_init_write():
    with msmseeder.utils.enter_temp_dir():
        project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
        test_data = {'test_field': 'test_value'}
        project_metadata.add_data(test_data)
        project_metadata.write()
        assert os.path.exists('meta0.yaml')
        assert yaml.load(open('meta0.yaml')) == {
            'init': {
                'test_field': 'test_value',
                'iteration': 0
            }
        }


def test_project_metadata_init_write_twice():
    with msmseeder.utils.enter_temp_dir():
        project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
        test_data = {'test_field': 'test_value'}
        project_metadata.add_data(test_data)
        project_metadata.write()
        project_metadata.write()
        assert os.path.exists('meta0.yaml')
        assert yaml.load(open('meta0.yaml')) == {
            'init': {
                'test_field': 'test_value',
                'iteration': 0
            }
        }
        assert os.path.exists('meta1.yaml')
        assert yaml.load(open('meta1.yaml')) == {
            'init': {
                'test_field': 'test_value',
                'iteration': 1
            }
        }


def test_project_metadata_add_prev_metadata():
    with msmseeder.utils.enter_temp_dir():
        os.mkdir('targets')
        init_project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
        test_data = {'test_field': 'test_value'}
        init_project_metadata.add_data(test_data)
        init_project_metadata.write()
        gather_targets_project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_targets')
        gather_targets_project_metadata.add_data(test_data)
        assert gather_targets_project_metadata.data == {
            'init': {
                'test_field': 'test_value',
                'iteration': 0
            },
            'gather_targets': {
                'test_field': 'test_value',
                'iteration': 0
            }
        }


def test_project_metadata_gather_targets_write():
    with msmseeder.utils.enter_temp_dir():
        os.mkdir('targets')
        init_project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
        test_data = {'test_field': 'test_value'}
        init_project_metadata.add_data(test_data)
        init_project_metadata.write()
        gather_targets_project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_targets')
        gather_targets_project_metadata.add_data(test_data)
        gather_targets_project_metadata.write()
        assert yaml.load(open(os.path.join('targets', 'meta0.yaml'))) == {
            'init': {
                'test_field': 'test_value',
                'iteration': 0
            },
            'gather_targets': {
                'test_field': 'test_value',
                'iteration': 0
            }
        }


def test_project_metadata_gather_templates_write():
    with msmseeder.utils.enter_temp_dir():
        os.mkdir('targets')
        os.mkdir('templates')
        init_project_metadata = msmseeder.core.ProjectMetadata(project_stage='init')
        test_data = {'test_field': 'test_value'}
        init_project_metadata.add_data(test_data)
        init_project_metadata.write()
        gather_targets_project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_targets')
        gather_targets_project_metadata.add_data(test_data)
        gather_targets_project_metadata.write()
        gather_templates_project_metadata = msmseeder.core.ProjectMetadata(project_stage='gather_templates')
        gather_templates_project_metadata.add_data(test_data)
        gather_templates_project_metadata.write()
        assert yaml.load(open(os.path.join('templates', 'meta0.yaml'))) == {
            'init': {
                'test_field': 'test_value',
                'iteration': 0
            },
            'gather_targets': {
                'test_field': 'test_value',
                'iteration': 0
            },
            'gather_templates': {
                'test_field': 'test_value',
                'iteration': 0
            }
        }