import ensembler.tests
import functools
import nose
import os
from pkg_resources import resource_filename


def expected_failure(test):
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise nose.SkipTest
        else:
            raise AssertionError(
                'A failure was expected, but this test appeared to pass. You may want to remove the expected_failure decorator.'
            )
    return inner


def get_installed_resource_filename(relative_path):
    """Get the full path to one of the reference files shipped for testing.
    In the source distribution, these files are in ``ensembler/tests/resources`` and
    ``ensembler/tests/example_project``, but on installation, they're moved to somewhere
    in the user's python site-packages directory.
    This function uses the pkg_resources package to find the file within the installation directory
    structure.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the ``ensembler/tests`` folder).

    Examples
    --------
    get_installed_resource_filename('example_project/meta0.yaml')
    """

    fn = resource_filename(ensembler.tests.__name__, relative_path)

    if not os.path.exists(fn):
        raise ValueError(
            "Sorry! {0} does not exist."
            "If you just added it, you'll have to re-install".format(relative_path)
        )

    return fn