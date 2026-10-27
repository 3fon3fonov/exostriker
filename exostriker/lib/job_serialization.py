"""Fast calculation snapshots, retaining dill support for custom session objects."""
import pickle
import io
import dill


class _SessionPickler(pickle.Pickler):
    def reducer_override(self, obj):
        if getattr(obj, '__module__', None) == '__main__':
            raise pickle.PicklingError('Interactive definitions require dill')
        return NotImplemented


def _pickle_dump(value, stream):
    _SessionPickler(stream, protocol=pickle.HIGHEST_PROTOCOL).dump(value)


def dumps(value):
    try:
        stream = io.BytesIO()
        _pickle_dump(value, stream)
        return stream.getvalue()
    except Exception:
        stream.close()
        return dill.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)


def dump(value, stream):
    start = stream.tell()
    try:
        _pickle_dump(value, stream)
    except Exception:
        stream.seek(start)
        stream.truncate()
        dill.dump(value, stream, protocol=pickle.HIGHEST_PROTOCOL)


load = dill.load
loads = dill.loads
