"""
Dear Gael,

I've been playing around with joblib lately. Amazing piece of software!

One of the most obvious applications, I thought, would be to
build very large data structures and have joblib managing the storage
of them for me so that I can pretend they are in memory and work
conveniently with them.

As I didn't see any direct support for this, I reasoned as follows.
The dump and load functions demand that I store one single
data structure in one file and keep track of the filename. That's
not suitable for large data as they have to be in memory and the
filename is beyond my interest.

I therefore though of making  the retrieve function below, with synonym
save, such that I for each piece of the data I can call save, joblib will
save the piece for me, and at any time I get a memory mapped (?) array
back and attach it in, e.g., a big dictionary that holds numerous such
views to arrays.

The amazing thing is that I can store the data in joblib and start
many interactive sessions or programs manipulating the data without
any need to recompute (i.e., restore pieces). All my big data are just
there!! Totally amazing :-)

Is this the right way to use joblib? Or is the same functionality available
in better ways?
"""

from joblib import Memory
memory = Memory(cachedir='tmp', verbose=0)

@memory.cache(ignore=['data'])
def retrieve(name, data=None):
    print 'joblib save of', name
    return data

save = retrieve

def create_big_data(n):
    """
    Create many pieces of a big data structure, one by one (here 2).
    Make joblib save the pieces.
    """
    import numpy
    big1 = numpy.linspace(0, 1, n+1)
    big2 = numpy.linspace(1, 2, n+1)

    save('a', big1)
    save('b', big2)

def retrieve_big_data():
    """
    Grab references to the pieces of a big data structure, all pieces
    stored by joblib, and collect them in some appropriate data
    structure, here a dict.
    """
    a = retrieve('a')
    b = retrieve('b')
    data = {'arrays': {'a': a, 'b': b}}
    return data

def compute(data):
    data['arrays']['b'][-10:] = 0.1
    return data

if __name__ == '__main__':
    create_big_data(1000)  # create pieces and store in joblib
    data = retrieve_big_data()
    data = compute(data)
    print data['arrays']['b'][-12:]

# Try to do different things with b and see if retrieve is ever
# called.
