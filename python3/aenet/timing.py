import timeit

t0 = 0.0
dt = 0.0

def timing(func):
    def wrap(*args, **kwargs):
        t0 = timeit.default_timer()
        result = func(*args, **kwargs)
        dt = timeit.default_timer() - t0
        print("call to `{}' - elapsed time (s): {}".format(
            func.__name__, dt))
        return result
    return wrap

def timethis(this):
    t0 = timeit.default_timer()
    result = this
    dt = timeit.default_timer() - t0
    print("elapsed time (s): {}".format(dt))
    return result

class Tick(object):
    def __init__(self):
        self.tick = False
    def __call__(self):
        if self.tick:
            dt = timeit.default_timer() - self.t0
            print("elapsed time since tick (s): {}".format(dt))
            self.tick = False
        else:
            self.t0 = timeit.default_timer()
            self.tick = True
