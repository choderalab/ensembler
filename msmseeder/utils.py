import os
import logging
import functools


def notify_when_done(fn):
    @functools.wraps(fn)
    def print_done(*args, **kwargs):
        try:
            import mpi4py.MPI
            comm = mpi4py.MPI.COMM_WORLD
            rank = comm.rank
            if rank == 0:
                fn(*args, **kwargs)
                logger = logging.getLogger('info')
                logger.info('Done.')
        except ImportError:
            fn(*args, **kwargs)
            logger = logging.getLogger('info')
            logger.info('Done.')

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


def file_exists_and_not_empty(filepath):
    if os.path.exists(filepath):
        if os.path.getsize(filepath) > 0:
            return True
    return False


def loglevel_setter(logger, loglevel):
    if loglevel is not None:
        loglevel_obj = getattr(logging, loglevel.upper())
        logger.setLevel(loglevel_obj)