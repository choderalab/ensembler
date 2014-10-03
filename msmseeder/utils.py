import os
import functools
import nose


def notify_when_done(fn):
    def print_done(*args, **kwargs):
        fn(*args, **kwargs)
        print 'Done'
    return print_done


def create_dir(dirpath):
    """
    :param dirpath: str
    """
    try:
        os.makedirs(dirpath)
        print 'Created directory "%s"' % dirpath
    except OSError as e:
        if e.errno == 17:
            print 'Directory "%s" already exists - will not overwrite' % dirpath
        else:
            raise


def expected_failure(test):
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise nose.SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner