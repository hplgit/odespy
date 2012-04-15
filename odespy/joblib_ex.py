from joblib import Memory
memory = Memory(cachedir='tmp', verbose=0)

@memory.cache(ignore=['data'])
def retrieve(name, data=None):
    print 'Inside retrieve'
    return data

save = retrieve

save('a', bigarray)
a = retrieve('a')

# Could make some empty dict with all keys and map nested keys to a name
# used for retrieve, can just be a tuple of keys

