import ensembler
from nose.plugins.attrib import attr


@attr('unit')
def test_project_metadata():
    project_metadata = ensembler.core.ProjectMetadata(project_stage='init')