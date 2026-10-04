"""Native build identity, independent of Exo-Striker's GUI imports."""
import ast
import hashlib
from pathlib import Path
import re


def build_identity(root=None):
    root = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    tree = ast.parse((root / '__init__.py').read_text())
    version = next(ast.literal_eval(node.value) for node in tree.body
                   if isinstance(node, ast.Assign)
                   and any(isinstance(target, ast.Name) and target.id == '__version__'
                           for target in node.targets))
    sources = [root / 'compile_swift.py', root / 'lib/swift/_build_info.py']
    sources.extend(path for path in (root / 'lib/swift').iterdir()
                   if path.suffix in ('.f95', '.f90', '.pyf', '.inc'))
    sources.extend(path for path in (root / 'source/swift_j').rglob('*')
                   if path.is_file() and path.suffix in ('.f', '.F', '.inc'))
    digest = hashlib.sha256()
    for path in sorted(sources):
        digest.update(path.relative_to(root).as_posix().encode() + b'\0')
        digest.update(path.read_bytes() + b'\0')
    tag = re.sub(r'[^a-zA-Z0-9_]', '_', str(version))
    return '_swift_memory_v' + tag + '_' + digest.hexdigest()[:20]
