import os


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


def file_exists_and_not_empty(filepath):
    if os.path.exists(filepath):
        if os.path.getsize(filepath) > 0:
            return True
    return False