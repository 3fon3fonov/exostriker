"""Build the in-memory SWIFT/F2PY extension from source.

Run: python exostriker/compile_swift.py [--debug] [--build-dir PATH]
Requires NumPy >= 1.26, Meson, Ninja, a C compiler and GNU Fortran >= 10.
"""
import argparse
import importlib.machinery
import os
from pathlib import Path
import re
import runpy
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
DEST = ROOT / 'lib' / 'swift'
CORE_DIRS = ('anal', 'coord', 'discard', 'io', 'mvs', 'orbel', 'obl',
             'util', 'helio', 'symba5', 'rmvs')
# Used only by tests/build_swift_reference.py, never by the library build.
DRIVERS = {
    'geninit': 'symba_f/geninit_j3_in_days.f',
    'symba': 'symba_f/swift_symba5_j.f',
    'mvs': 'mvs_f/swift_mvs_j.f',
    'mvs_gr': 'mvs_f/swift_mvs_j_GR.f',
    'follow_symba': 'symba_f/follow_symba2.f',
    'follow_mvs': 'mvs_f/follow2.f',
}


def fixed_statements(source):
    """Group fixed-form continuation lines, preserving statement labels."""
    statements = []
    for raw in source.expandtabs(8).splitlines():
        if not raw.strip() or raw[0] in 'cC*!#':
            continue
        line = raw.ljust(6)
        if line[5] not in (' ', '0'):
            if not statements:
                raise ValueError('Orphan Fortran continuation')
            statements[-1].append(raw)
        else:
            statements.append([raw])
    return statements


def statement_text(lines):
    return ' '.join(line[6:].strip().split('!')[0] for line in lines).strip().lower()


def expand_includes(source, parent):
    def replace(match):
        include = (parent / match.group(1)).resolve()
        return '\n' + expand_includes(include.read_text().expandtabs(8), include.parent)
    return re.sub(r"include\s+'([^']+)'", replace, source, flags=re.I)


def prepare_core(path, fc):
    source = path.read_text().expandtabs(8)
    if path.suffix == '.F':
        source = subprocess.check_output([fc, '-E', '-cpp', '-P', '-D_RECUR_SUB',
                                          str(path)], text=True)
    statements = fixed_statements(source)
    output = []
    for lines in statements:
        text = statement_text(lines)
        if text.startswith('include') and 'swift.inc' in text:
            if path.stem == 'drift_kepu_stumpff':
                output.append('      use, intrinsic :: ieee_arithmetic,\n'
                              '     &    only: ieee_is_finite\n')
            output.append('      use swift_memory_support, only: swift_status,\n'
                          '     &    swift_poll, swift_remove\n')
        if re.match(r'(open|read|close|stop)\b', text):
            raise ValueError(f'Unexpected process/file operation in {path}: {text}')
        if text.startswith('write(') or text.startswith('write ('):
            # Numerical diagnostics are returned as status, not console I/O.
            continue
        if text.startswith('call util_exit('):
            output.append('      swift_status=3\n      return\n')
            continue
        if path.stem == 'symba5_step_recur' and text.startswith('dtl ='):
            output.append('      call swift_poll()\n'
                          '      if(swift_status.ne.0) return\n'
                          '      if(ireci.gt.64) then\n'
                          '         swift_status=4\n         return\n      endif\n')
        if path.stem == 'drift_kepu_stumpff' and text.startswith('do while'):
            # Infinity never shrinks under repeated division: reject it before
            # the original range-reduction loop instead of hanging Python.
            output.append('      if(.not.ieee_is_finite(x)) then\n'
                          '         swift_status=4\n         return\n      endif\n')
        if path.stem == 'discard_mass_reorder5' and text.startswith('do i=ip,'):
            output.append('      call swift_remove(ip)\n')
        if path.stem == 'discard_massive5' and text.startswith('do i=1,nbod'):
            if '      isperih=0\n' not in output:
                output.append('      isperih=0\n')
        output.append('\n'.join(lines)+'\n')
        if re.search(r'\bcall\s+', text):
            output.append('      if(swift_status.ne.0) return\n')
    return expand_includes(''.join(output), path.parent)


def selected_sources():
    """Only the numerical dependency closure; no original I/O or exit code."""
    core = ROOT / 'source' / 'swift_j'
    candidates = {}
    for folder in CORE_DIRS:
        if folder == 'io':
            continue
        for path in (core / folder).rglob('*'):
            if path.suffix in ('.f', '.F') and not path.stem.endswith('_CPP') and path.stem != 'util_exit':
                candidates[path.stem.lower()] = path
    selected = {}
    pending = ['\n'.join(line.split('!')[0] for line in (DEST / 'swift.f95').read_text().splitlines())]
    while pending:
        text = pending.pop()
        tokens = set(re.findall(r'\b(\w+)\s*\(', text.lower()))
        for name in sorted(tokens & candidates.keys() - selected.keys()):
            path = candidates[name]
            selected[name] = path
            # Strip comments for dependency discovery. Recursive declaration
            # preprocessing isn't needed: the file is already selected by name.
            pending.append('\n'.join(line for line in path.read_text().splitlines()
                                     if line and line[0] not in 'cC*!#'))
    return [selected[name] for name in sorted(selected)]


def compile_python_extension(build_dir=None, debug=False):
    import numpy as np
    if np.lib.NumpyVersion(np.__version__) < '1.26.0':
        raise RuntimeError('The SWIFT Meson build requires NumPy >= 1.26')
    fc, ar = os.environ.get('FC', 'gfortran'), os.environ.get('AR', 'ar')
    for command in (fc, ar, 'meson', 'ninja'):
        if shutil.which(command) is None:
            raise RuntimeError(f'Required build tool not found: {command}')
    if build_dir is None:
        with tempfile.TemporaryDirectory(prefix='exostriker-swift-build-') as temp:
            return compile_python_extension(temp, debug)
    build = Path(build_dir).resolve()
    build.mkdir(parents=True, exist_ok=True)
    flags = ['-O2', '-std=legacy', '-ffixed-line-length-none', '-ffree-line-length-none',
             '-fallow-argument-mismatch']
    if sys.platform != 'win32':
        flags.append('-fPIC')
    if debug:
        flags += ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    support = build / 'swift_support.f90'
    shutil.copyfile(DEST / support.name, support)
    objects = []
    # Compile the support module first so all numerical files can USE it.
    support_obj = build / 'swift_support.o'
    subprocess.run([fc, *flags, '-c', support.name, '-o', support_obj.name], cwd=build, check=True)
    objects.append(str(support_obj))
    sources = selected_sources()
    (build / 'numerical_sources.txt').write_text('\n'.join(str(p.relative_to(ROOT)) for p in sources)+'\n')
    for path in sources:
        target = build / (path.stem + '.f')
        target.write_text(prepare_core(path, fc))
        obj = target.with_suffix('.o')
        subprocess.run([fc, *flags, '-c', target.name, '-o', obj.name], cwd=build, check=True)
        objects.append(str(obj))
    # swift.f95 contains the actual array-based integration entry points.
    target = build / 'swift.f95'
    target.write_text(expand_includes((DEST / 'swift.f95').read_text(), DEST))
    obj = build / 'swift.o'
    subprocess.run([fc, *flags, '-c', target.name, '-o', obj.name], cwd=build, check=True)
    objects.append(str(obj))
    archive = build / 'libswift_memory.a'
    archive.unlink(missing_ok=True)
    for offset in range(0, len(objects), 25):
        subprocess.run([ar, 'rcs', str(archive), *objects[offset:offset+25]], check=True)
    module_name = runpy.run_path(str(DEST / '_build_info.py'))['build_identity']()
    signature = (DEST / 'swift.pyf').read_text().replace('_swift_memory', module_name)
    (build / 'swift.pyf').write_text(signature)
    env = os.environ.copy()
    env['FC'] = fc
    env['FFLAGS'] = ' '.join(flags)
    env['LDFLAGS'] = env.get('LDFLAGS', '') + ' ' + shlex.quote('-L'+build.as_posix())
    if sys.platform == 'win32':
        env.setdefault('CC', 'gcc')
    subprocess.run([sys.executable, '-m', 'numpy.f2py', '-c', 'swift.pyf',
                    '-lswift_memory', '--backend', 'meson'], cwd=build, env=env, check=True)
    extensions = [p for p in build.iterdir() if p.name.startswith(module_name)
                  and any(p.name.endswith(s) for s in importlib.machinery.EXTENSION_SUFFIXES)]
    if len(extensions) != 1:
        raise RuntimeError(f'Expected one extension, found {extensions}')
    subprocess.run([sys.executable, '-c',
                    'import runpy,sys; runpy.run_path(sys.argv[1])["load_native"](build_if_missing=False)',
                    str(DEST / '_load.py')], cwd=build, env=env, check=True)
    target = DEST / extensions[0].name
    with tempfile.NamedTemporaryFile(prefix='._swift_memory-', dir=DEST, delete=False) as staged:
        staged_path = Path(staged.name)
    try:
        shutil.copy2(extensions[0], staged_path)
        os.replace(staged_path, target)
    finally:
        staged_path.unlink(missing_ok=True)
    print(f'Built in-memory SWIFT library: {target}')
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', help='Retain generated sources and objects here')
    parser.add_argument('--debug', action='store_true', help='Enable Fortran runtime checks')
    args = parser.parse_args()
    compile_python_extension(args.build_dir, args.debug)
