from joblib import Memory
memory = Memory(cachedir='tmp', verbose=0)

@memory.cache(ignore=['data'])
def retrieve(name, data=None):
    print 'joblib save of', name
    return data

save = retrieve

def create_big_data(n):
    import numpy
    big1 = numpy.linspace(0, 1, n+1)
    big2 = numpy.linspace(1, 2, n+1)

    save('a', big1)
    save('b', big2)

    a = retrieve('a')
    b = retrieve('b')
    data = {'arrays': {'a': a, 'b': b}}
    return data

def compute(data):
    data['arrays']['b'][-10:] = 0.1
    return data

data = create_big_data(1000)  # create pieces and store in joblib
data = compute(data)
print data['arrays']['b'][-12:]


