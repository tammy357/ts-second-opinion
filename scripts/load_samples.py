"""Load GoFlow samples_all.pkl without torch: stub torch/PyG classes with numpy."""
import pickle, io, zipfile, sys
import numpy as np

DT = {'FloatStorage': np.float32, 'DoubleStorage': np.float64, 'LongStorage': np.int64,
      'IntStorage': np.int32, 'BoolStorage': np.bool_, 'ByteStorage': np.uint8,
      'HalfStorage': np.float16, 'ShortStorage': np.int16, 'CharStorage': np.int8,
      'BFloat16Storage': np.uint16}

class StorageType:
    def __init__(self, name): self.name = name

class InnerUnpickler(pickle.Unpickler):
    def __init__(self, f, zf, prefix):
        super().__init__(f); self.zf = zf; self.prefix = prefix
    def find_class(self, mod, name):
        if mod == 'torch' and name.endswith('Storage'): return StorageType(name)
        if mod == 'torch' and name in ('float32','float64','int64','bool'): return name
        if mod == 'collections' and name == 'OrderedDict':
            import collections; return collections.OrderedDict
        if mod == 'torch._utils' and name == '_rebuild_tensor_v2': return rebuild_tensor
        return super().find_class(mod, name)
    def persistent_load(self, pid):
        typename, stype, key, location, numel = pid
        raw = self.zf.read(f'{self.prefix}/data/{key}')
        name = stype.name if isinstance(stype, StorageType) else str(stype)
        return np.frombuffer(raw, dtype=DT.get(name, np.float32))

class LegacyUnpickler(pickle.Unpickler):
    def __init__(self, f):
        super().__init__(f); self.pids = {}
    def find_class(self, mod, name):
        if mod == 'torch' and name.endswith('Storage'): return StorageType(name)
        if mod.startswith('torch'): return lambda *a, **k: None
        return super().find_class(mod, name)
    def persistent_load(self, pid):
        # ('storage', storage_type, root_key, location, numel, view_metadata)
        self.pids[pid[2]] = pid[1].name if isinstance(pid[1], StorageType) else 'FloatStorage'
        return ('__storage__', pid[2])

def load_from_bytes(b):
    if b[:2] == b'PK':
        zf = zipfile.ZipFile(io.BytesIO(b))
        pk = [n for n in zf.namelist() if n.endswith('data.pkl')][0]
        prefix = pk.rsplit('/', 1)[0]
        return InnerUnpickler(io.BytesIO(zf.read(pk)), zf, prefix).load()
    f = io.BytesIO(b)
    pickle.load(f); pickle.load(f); pickle.load(f)  # magic, protocol, sys_info
    up = LegacyUnpickler(f); obj = up.load()
    keys = pickle.load(f)
    out = {}
    for k in keys:
        n = int(np.frombuffer(f.read(8), dtype=np.int64)[0])
        dt = np.dtype(DT.get(up.pids.get(k, 'FloatStorage'), np.float32))
        out[k] = np.frombuffer(f.read(n * dt.itemsize), dtype=dt)
    if isinstance(obj, tuple) and obj and obj[0] == '__storage__': return out[obj[1]]
    return obj

def rebuild_tensor(storage, offset, size, stride, *args):
    if len(size) == 0: return storage[offset]
    itemsize = storage.itemsize
    return np.lib.stride_tricks.as_strided(storage[offset:], shape=tuple(size),
                                           strides=tuple(s*itemsize for s in stride)).copy()

class Obj:
    def __init__(self, *a, **k): self.__dict__['_store'] = {}
    def __setstate__(self, state):
        self.__dict__['_state'] = state
    def __getattr__(self, k):
        st = self.__dict__.get('_state', {})
        for v in (st.values() if isinstance(st, dict) else []):
            if isinstance(v, Obj) and k in v.__dict__.get('_state', {}): return v._state[k]
            if isinstance(v, dict) and k in v: return v[k]
        raise AttributeError(k)
    def keys(self):
        st = self.__dict__.get('_state', {})
        out = []
        for v in st.values():
            if isinstance(v, Obj): out += list(v.__dict__.get('_state', {}).keys())
            elif isinstance(v, dict): out += list(v.keys())
        return out

class OuterUnpickler(pickle.Unpickler):
    def find_class(self, mod, name):
        if mod == 'torch.storage' and name == '_load_from_bytes': return load_from_bytes
        if mod == 'torch._utils' and name == '_rebuild_tensor_v2': return rebuild_tensor
        if mod.startswith('torch_geometric') or name in ('DataEdgeAttr','DataTensorAttr'): return Obj
        if mod.startswith('torch'): return lambda *a, **k: None
        return super().find_class(mod, name)

def load(path):
    with open(path, 'rb') as f: return OuterUnpickler(f).load()

if __name__ == '__main__':
    d = load(sys.argv[1])
    print(type(d), len(d) if hasattr(d, '__len__') else '')
    x = d[0] if isinstance(d, (list, tuple)) else d
    print(type(x))
    if isinstance(x, Obj):
        st = x._state
        print('state keys', list(st.keys()) if isinstance(st, dict) else type(st))
        for k in x.keys():
            v = getattr(x, k)
            print(' ', k, type(v).__name__, getattr(v, 'shape', ''), (v if not hasattr(v,'shape') else ''))
