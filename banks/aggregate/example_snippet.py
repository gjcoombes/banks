# This snippet shows the code I commonly use time functions
# Also see for timing and profiling in IPython

import time

def func_with_internal_timer(n):
    """The disadvantage here is that you need the timer code inside each func."""
    start = time.time()

    # Do the work of the function here e.g.
    time.sleep(n)

    elapsed = time.time() - start
    print("{} finished in {:2.3f} seconds".format("func_with_internal_timer", elapsed))

# Or we can make a decorator function which wraps any function you 'decorate'

def timer(func):
    """This is the decorator function """
    def inner(*args, **kwargs):
        start = time.time()
        out = func(*args, **kwargs)
        elapsed = time.time() - start
        print("{} finished in {:2.3f} seconds".format(func.__name__, elapsed))
        return out
    return inner

# Now we can annotate any with the decorator syntax: @decorator

@timer
def ordinary_func(n):
    time.sleep(n)

func_with_internal_timer(1)
ordinary_func(1)

### IPython has really useful magic functions for profiling like
### %time, %timeit, %prun, %lprun
### See more here http://pynash.org/2013/03/06/timing-and-profiling.html
