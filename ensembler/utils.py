import contextlib
import os
import logging
import functools
import shutil
import tempfile
from ensembler.core import mpistate

logger = logging.getLogger('info')


def nonefn():
    return None


def mpirank0only(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        if mpistate.rank == 0:
            fn(*args, **kwargs)
    return wrapper


def mpirank0only_and_end_with_barrier(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        if mpistate.rank == 0:
            fn(*args, **kwargs)
        mpistate.comm.Barrier()
    return wrapper


def notify_when_done(fn):
    @functools.wraps(fn)
    def print_done(*args, **kwargs):
        fn(*args, **kwargs)
        log_done()
    return print_done


@mpirank0only
def log_done():
    logger.info('Done.')


def create_dir(dirpath, quiet=True):
    """
    :param dirpath: str
    """
    try:
        os.makedirs(dirpath)
        if not quiet:
            logger.info('Created directory "%s"' % dirpath)
    except OSError as e:
        if e.errno == 17:
            logger.debug('Directory "%s" already exists - will not overwrite' % dirpath)
        else:
            raise


def file_exists_and_not_empty(filepath):
    if os.path.exists(filepath):
        if os.path.getsize(filepath) > 0:
            return True
    return False


def loglevel_setter(logger, loglevel):
    if loglevel is not None:
        loglevel_obj = getattr(logging, loglevel.upper())
        logger.setLevel(loglevel_obj)


@contextlib.contextmanager
def enter_temp_dir():
    """Create a temporary directory, enter, yield, exit, rmdir; used as context manager."""
    temp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(temp_dir)
    yield temp_dir
    os.chdir(cwd)
    shutil.rmtree(temp_dir)


def debug_method(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            fn(*args, **kwargs)
        except Exception as e:
            print e
            import traceback
            print traceback.format_exc()
            import ipdb; ipdb.set_trace()
    return wrapper


def set_arg_with_default(arg, default_arg):
    if arg is None:
        arg = default_arg
    return arg