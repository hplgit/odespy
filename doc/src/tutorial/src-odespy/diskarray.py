import numpy as np
class DiskArray:
    def __init__(self, size, in_memory=3):
        if isinstance(size, int):
            stored_size = tuple(size)
            self.idxmap = np.zeros(size, dtype=np.int)
        elif isinstance(size, (tuple,list)):
            stored_size = list(size)
            stored_size[0] = in_memory
            self.idxmap = np.zeros(size[0], dtype=np.int)
        self.array = np.zeros(stored_size)
        self.stored = self.array.size
        # shelve

    def _next(self, idx):
        idx = tuple(idx)  # can be a slice!
        idx0 = idx[0]
        if idx0 < self.stored-1:
            return idx0 + 1
        elif:



