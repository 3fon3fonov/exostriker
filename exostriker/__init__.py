__version__ = "0.96.8"
__author__ = "Trifon Trifonov"

import importlib
from pathlib import Path
import sys

# Legacy numerical modules use absolute imports from the bundled lib directory.
# Make them available without importing Qt or constructing the GUI's global fit.
_lib_path = str(Path(__file__).resolve().parent / 'lib')
if _lib_path not in sys.path:
    sys.path.insert(0, _lib_path)


def __getattr__(name):
    if name in ('gui', 'lib'):
        module = importlib.import_module(f'{__name__}.{name}')
        globals()[name] = module
        return module
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
