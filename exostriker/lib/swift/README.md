# In-memory SWIFT Python library

`RV_mod.functions.run_stability()` calls `swift.integrate()`, backed by the
versioned F2PY extension `_swift_memory_v<version>_<source-hash>`. SyMBA, MVS and MVS+GR now execute **inside the
calling Python process**. Initial conditions, integration state, orbital
histories, energy and angular momentum all stay in memory. Integration uses
no executables, Python workers, temporary directories or intermediate files.

The build now compiles **`swift.f95`**. Its initializer and Jacobi/MVS/GR
mathematical helpers come from the original prototype. The unfinished
file-based integration and follower routines have been replaced by a corrected
array interface. The accompanying original SWIFT numerical routines are linked
into the same Python extension; they are not separate installed programs.

## Build and use

Use the Python environment that runs Exo-Striker, from the repository root:

```sh
python -m pip install "numpy>=1.26" meson ninja
python exostriker/compile_swift.py
python -m unittest discover -s tests -p test_swift_library.py -v
```

After a pip upgrade or Git pull, restart Exo-Striker. At startup SWIFT checks
the Exo-Striker version and a fingerprint of its Fortran sources, includes,
interface and build script. A matching binary is reused; a missing or outdated
one is rebuilt with the running Python interpreter. Git source changes trigger
a rebuild even when the application version stays the same. Direct calls to
`integrate()` also build a missing current library. A failed build is reported;
the loader never falls back to an outdated binary. Automatic builds require
the compiler/tools listed above and write access to the installed package.

Binary names include the version, source fingerprint and Python ABI suffix.
Older versioned binaries and the previous `_swift_memory`, `_swift_native`
and `swift.cpython-*.so` binaries are ignored. The fingerprint is checked once
per process so integrations themselves keep their no-file-I/O behavior.

```python
from exostriker.lib.swift import integrate

result = integrate(
    1.0,  # stellar mass, solar masses
    [[0.001, 1.0, 0.05, 90.0, 10.0, 0.0, 20.0]],
    timemax=100.0, timestep=2.0, output_step=1.0, integrator="symba",
)
energy = result["energy"]
orbit = result["orbits"][0]
```

For GUI-free use from a checkout, put `exostriker/lib` on Python's import path
and `from swift import integrate`. The existing top-level `exostriker` package
imports GUI dependencies. A runnable example is `swift_test.py` in this folder.

## Data interface

Planet rows contain mass in solar masses, Jacobi semimajor axis in AU,
eccentricity, inclination, argument of periastron, node and mean anomaly.
Angles are degrees. Supply planets in increasing period order; `run_stability`
sorts enabled GUI slots and fills the existing `evol_*` fields in that order.
The interface supports 1–50 planets; the GUI still has its existing nine slots.

`timemax` and `output_step` are years; `timestep` is days. Fixed output steps
use `365.25 * timemax / fixed_steps`. The original SWIFT cadence is preserved:
output occurs after a completed integration step crossing a requested output
time, and the loop can overshoot `timemax` by one step. A single step produces
at most one output even if it crosses several requested output times.

Energy columns: time (days), total energy, Lx, Ly, Lz.
Orbit columns: time (days), body ID, a (AU), e, inclination, node, argument of
periastron, mean anomaly, pericentre, apocentre, planet/star mass ratio.
All arrays are double precision, including orbital elements that the standalone
programs used to round through float32 binary files and formatted text.
The adapter converts times to years. Initial and final recorded samples are
retained. Discarded histories end at their last valid sample. Merged planets
follow the surviving body, matching the original SyMBA follower convention.
The initial frame always contains the actual initial conditions: unlike the
legacy file follower, a merger tagged at time zero does not retroactively
remap that frame to a different planet.

`save_dir=PATH` is an optional **export of completed results**: `energy.out`,
`pl_N.out` and `result.npz`. The directory must not exist. This is the only
integration API option that writes files. There are no `param.in`, `pl.in`,
`bin.dat`, discard files or workers even when exporting.
In `run_stability`, `remove_stab_save_dir=False` selects this export under
`stability/stab_save_dir` relative to the caller, or an absolute supplied path.
The default performs no output file I/O and does not change the working directory.

## State, errors and timeouts

`swift_begin` validates and initializes a session from arrays. `swift_advance`
runs bounded groups of steps and returns batches of samples with validity masks
and explicit filled counts. Python concatenates these in memory. This avoids
allocating a guessed full-run output array up front, and handles early discards
and timeouts without exposing uninitialized output rows.

The legacy numerical routines have saved work arrays. A process-wide lock
serializes entire sessions, including when the package is imported under its
two supported names. Cached accelerations and body tracking are reset at each
new run. Concurrent callers are safe but **do not integrate simultaneously in
one process**. F2PY releases the GIL during native calls so other Python threads
can run. The private extension entry points must not be called outside this
wrapper while another session is active.

Reachable `util_exit` calls become status returns with propagation through the
calling routines. The linked numerical subset has no Fortran process-exit or
file-I/O code. Non-finite states are reported as errors; excessive encounter
recursion is stopped, and non-finite Stumpff arguments are rejected before
range reduction can loop indefinitely.

Timeouts are **cooperative**: Fortran checks its monotonic system clock every
64 outer steps and on each encounter-recursion entry. Python checks between
batches as well. A timeout raises `TimeoutError`; numerical problems raise
`SwiftError`. A run returns only completed histories, and a failed run does not
update the GUI fit object's evolution arrays. The next run reinitializes state.
Timeout includes waiting for the session lock, but excludes compilation and
optional result export. It is not an operating-system hard kill deadline.

## Build implementation and platforms

Requirements: Python 3.9+, NumPy 1.26+, Meson, Ninja, a C compiler, GNU Fortran
10+ and GNU-compatible `ar`. `FC`, `CC` and `AR` can select tool executables.

- Linux: install GCC/GFortran and Python development headers.
- macOS: install Xcode command-line tools and Homebrew GCC. Set `FC` to the
  versioned GFortran executable if needed. Match Python's arm64/x86_64 architecture.
- Windows: use native MinGW-w64 GCC/GFortran and `ar`, matching Python's
  architecture, on PATH. Use `CC=gcc` and `FC=gfortran`. Do not mix MSVC archives
  or a Cygwin compiler into this build. The loader registers PATH directories
  for dependent compiler DLLs; those runtime DLLs must also be present when
  launching Exo-Striker.

The build always compiles sources, never bundled precompiled archives. It
selects the numerical dependency closure, creates library-safe copies of the
fixed-form sources in the build directory, and compiles `swift_support.f90`
and `swift.f95`. The original numerical source files remain available for
independent standalone comparisons. The generated core copies remove console
writes, add error propagation/timeout checks, and hook body compaction into
in-memory tracking. Mathematical stepping equations are retained. The explicit
`swift.pyf` exposes only the two validated session entry points.

`--debug` enables runtime/bounds checks. `--build-dir PATH` retains generated
sources, `numerical_sources.txt`, objects and the archive for inspection.
Extension filenames use Python's actual platform suffix, and builds support
paths with spaces. Rebuild after changing Python, NumPy's ABI, architecture,
compiler or sources. This does not build distributable binary wheels or bundle
compiler runtime DLLs.

The GitHub Actions workflow builds and tests Linux, macOS and Windows on Python
3.11 and 3.13. macOS and Windows still need actual CI runs before they can be
claimed as verified. Local validation is on Linux x86_64, Python 3.12 and GNU
Fortran 13.3: all 17 tests pass with NumPy 1.26.4 and 2.5.3, with checked and
optimized builds exercised. The 100-year example also completes successfully.

## Verification

Tests cover the three integrators, conservation, CLI numerical comparisons,
multiple output batches, single planets, all planets discarded, mergers and
body remapping, mixed repeated/concurrent sessions, timeouts and recovery,
numerical-error recovery, exports, invalid inputs and the actual GUI adapter.
A dedicated audit blocks Python process and file-creation APIs during integration
and verifies that the working directory stays empty.

For an independent build of the original standalone sources:

```sh
python tests/build_swift_reference.py /tmp/swift-ref-build /tmp/swift-reference
SWIFT_REFERENCE_DIR=/tmp/swift-reference python -m unittest discover -s tests -p test_swift_library.py -v
```

The reference initializes the original GR driver's uninitialized counter.
Comparison tolerances account for its float32/text output precision; the new
arrays retain full double precision. These checks validate integration and
selected physical scenarios, not every long chaotic encounter or merger chain.
The original MVS behavior is retained: its discard routine applies to test
particles, and this interface initializes no test particles. Arbitrary-orbit
runs now use the same library API. MVS_GR stops when a planet becomes
non-elliptic because its orbit-averaged GR correction requires a bound ellipse.
This returns the completed output samples with `completed=False` and an
explanatory `status`, rather than raising a generic numerical error. These
histories end at the last output sample before the stop, not the exact event
time. Successful runs return `completed=True` and `status='finished'`.
