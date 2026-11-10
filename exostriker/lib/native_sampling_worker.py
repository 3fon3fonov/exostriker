"""Install immutable likelihood inputs once per sampling worker.

Task messages contain parameter vectors, rather than the complete data context.
"""
_CALLBACK = None


def initialize(callback, args):
    global _CALLBACK
    _CALLBACK = (callback, args)


def evaluate(parameters):
    if _CALLBACK is None:
        raise RuntimeError('Native sampling worker was not initialized.')
    callback, args = _CALLBACK
    return callback(parameters, *args)
