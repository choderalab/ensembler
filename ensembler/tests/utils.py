import functools
import nose


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