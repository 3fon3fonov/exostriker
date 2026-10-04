"""Load the native extension and retain Windows compiler DLL directories."""
import importlib
import os
from pathlib import Path
import sys
import subprocess
import runpy
import threading
import types

_DLL_HANDLES = []
# Exo-Striker can import this package both as swift and exostriker.lib.swift.
# The same shared object must have the same lock under either import name.
_registry = sys.modules.setdefault('_exostriker_swift_locks', types.ModuleType('_exostriker_swift_locks'))
_locks = _registry.__dict__.setdefault('locks', {})
SESSION_LOCK = _locks.setdefault(str(Path(__file__).resolve().parent), threading.Lock())


_BUILD_LOCK = _locks.setdefault(str(Path(__file__).resolve().parent) + ':build', threading.Lock())
_MODULE_NAME = runpy.run_path(str(Path(__file__).with_name('_build_info.py')))['build_identity']()


def load_native(build_if_missing=True):
    if sys.platform == 'win32' and not _DLL_HANDLES:
        for directory in dict.fromkeys(os.environ.get('PATH', '').split(os.pathsep)):
            if directory and Path(directory).is_dir():
                try:
                    _DLL_HANDLES.append(os.add_dll_directory(directory))
                except OSError:
                    pass
    package = __package__
    name = '.' + _MODULE_NAME if package else _MODULE_NAME
    with _BUILD_LOCK:
        try:
            return importlib.import_module(name, package)
        except (ImportError, OSError):
            if not build_if_missing:
                raise
            builder = Path(__file__).resolve().parents[2] / 'compile_swift.py'
            print('Building SWIFT for the current Exo-Striker version and sources...')
            try:
                subprocess.run([sys.executable, str(builder)], check=True)
            except (OSError, subprocess.CalledProcessError) as exc:
                raise ImportError(f'Cannot build the current SWIFT library: {builder}') from exc
            importlib.invalidate_caches()
            try:
                return importlib.import_module(name, package)
            except (ImportError, OSError) as exc:
                raise ImportError('SWIFT build did not produce an importable current library') from exc
