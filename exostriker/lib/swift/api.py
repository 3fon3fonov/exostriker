"""In-process, in-memory SWIFT integration backed by swift.f95/F2PY.

A lock serializes access to legacy Fortran SAVE arrays. The native calls release
Python's GIL; timeout checks return through Fortran, never terminating Python.
"""
import math
from pathlib import Path
import threading
import time

import numpy as np

if __package__:
    from ._load import load_native, SESSION_LOCK
else:  # Direct execution of the small swift_test.py example.
    from _load import load_native, SESSION_LOCK

BUILD = Path(__file__).resolve().parents[2] / "compile_swift.py"
_LOCK = SESSION_LOCK


class SwiftError(RuntimeError):
    """A SWIFT build, numerical or interface error."""


def is_available():
    """Check the array API without launching a worker or executable."""
    try:
        native = load_native(build_if_missing=False)
        return hasattr(native, 'swift_begin') and hasattr(native, 'swift_advance')
    except (ImportError, OSError):
        return False


def _check_status(status):
    if status == 2:
        raise TimeoutError('SWIFT exceeded its time limit')
    messages = {1: 'Invalid native SWIFT arguments or session state',
                3: 'SWIFT could not complete a Kepler drift or encounter step',
                4: 'SWIFT encountered a singular/non-finite state or excessive encounter recursion',
                5: 'MVS_GR stopped because a planet became non-elliptic'}
    if status:
        raise SwiftError(messages.get(status, f'SWIFT error {status}'))


def integrate(stellar_mass, planets, *, timemax=3000.0, timestep=10.0,
              timeout_sec=1000.0, integrator="symba", output_step=1.0,
              use_fixed_steps=False, fixed_steps=1000, gr_step=1,
              save_dir=None):
    """Return energy and orbital histories as double-precision NumPy arrays.

    Planet rows: mass (solar masses), Jacobi a (AU), e, inclination, argument
    of periastron, node, mean anomaly (degrees), sorted by increasing period.
    timemax/output_step are years; timestep is days. There is no file I/O or
    subprocess during integration. save_dir optionally exports completed arrays
    to a new directory. Timeout includes waiting for the native-session lock.

    Energy columns: time (days), total energy, Lx, Ly, Lz. Orbit columns:
    time (days), id, a, e, i, node, omega, M, pericentre, apocentre, mass ratio.
    Discarded histories end at the last valid sample; merged bodies follow
    the survivor, as in the original SyMBA follower.
    """
    planets = np.asarray(planets, dtype=float)
    if planets.ndim != 2 or planets.shape[1] != 7 or not 1 <= len(planets) <= 50:
        raise ValueError("planets must have shape (1..50, 7)")
    if not np.all(np.isfinite(planets)):
        raise ValueError("Planet parameters must be finite")
    if np.any(planets[:, :2] <= 0) or np.any((planets[:, 2] < 0) | (planets[:, 2] >= 1)):
        raise ValueError("Require positive planet masses and a, and 0 <= e < 1")
    if integrator not in ("symba", "mvs", "mvs_gr"):
        raise ValueError("integrator must be symba, mvs or mvs_gr")
    values = dict(stellar_mass=stellar_mass, timemax=timemax, timestep=timestep,
                  timeout_sec=timeout_sec)
    for name, value in values.items():
        if not math.isfinite(float(value)) or float(value) <= 0:
            raise ValueError(f"{name} must be finite and positive")
    duration = float(timemax) * 365.25
    if not math.isfinite(duration):
        raise ValueError("timemax is too large")
    if use_fixed_steps:
        if not math.isfinite(float(fixed_steps)) or int(fixed_steps) != fixed_steps or fixed_steps <= 0:
            raise ValueError("fixed_steps must be a positive integer")
        interval = duration / int(fixed_steps)
    else:
        interval = float(output_step) * 365.25
    if not math.isfinite(interval) or interval <= 0:
        raise ValueError("output_step must be finite and positive")
    if integrator == "mvs_gr" and (not math.isfinite(float(gr_step)) or int(gr_step) != gr_step
                                    or not 1 <= gr_step <= 2147483647):
        raise ValueError("gr_step must be an integer between 1 and 2147483647")
    folder = Path(save_dir).resolve() if save_dir is not None else None
    if folder is not None and folder.exists():
        raise FileExistsError(folder)
    try:
        native = load_native()
    except (ImportError, OSError) as exc:
        raise SwiftError(f"Build the in-memory SWIFT library with: python {BUILD}") from exc
    deadline = time.monotonic() + float(timeout_sec)
    # Copies prevent another Python thread from modifying caller-owned inputs
    # while the GIL is released in Fortran.
    parameters = np.array(planets, dtype=np.float64, order='F', copy=True)
    if not _LOCK.acquire(timeout=min(float(timeout_sec), threading.TIMEOUT_MAX)):
        raise TimeoutError('SWIFT timed out waiting for another integration')
    try:
        remaining = deadline-time.monotonic()
        if remaining <= 0:
            raise TimeoutError('SWIFT exceeded its time limit')
        mode = {'symba': 1, 'mvs': 2, 'mvs_gr': 3}[integrator]
        _check_status(native.swift_begin(float(stellar_mass), parameters, mode,
                      duration, float(timestep), interval,
                      int(gr_step) if mode == 3 else 1, remaining))
        energy_parts = []
        orbit_parts = [[] for _ in planets]
        done = False
        stopped_nonelliptic = False
        while not done:
            energy, orbits, valid, count, done, status = native.swift_advance(len(planets), 2048, 256)
            stopped_nonelliptic = status == 5 and integrator == 'mvs_gr'
            if not stopped_nonelliptic:
                _check_status(status)
            if time.monotonic() >= deadline:
                raise TimeoutError('SWIFT exceeded its time limit')
            if count:
                energy_parts.append(energy[:count].copy())
                for index in range(len(planets)):
                    mask = valid[:count, index] != 0
                    if np.any(mask):
                        orbit_parts[index].append(orbits[:count, index, :][mask].copy())
            if stopped_nonelliptic:
                break
        result = {'energy': np.concatenate(energy_parts),
                  'orbits': [np.concatenate(parts) if parts else np.empty((0, 11))
                             for parts in orbit_parts],
                  'completed': not stopped_nonelliptic,
                  'status': ('MVS_GR stopped: a planet became non-elliptic; '
                             'the bound-orbit GR correction no longer applies'
                             if stopped_nonelliptic else 'finished')}
    finally:
        _LOCK.release()
    if folder is not None:
        folder.mkdir(parents=True, exist_ok=False)
        np.savetxt(folder / 'energy.out', result['energy'])
        arrays = {'energy': result['energy']}
        for index, orbit in enumerate(result['orbits']):
            np.savetxt(folder / f'pl_{index+1}.out', orbit)
            arrays[f'planet_{index}'] = orbit
        np.savez(folder / 'result.npz', **arrays)
    return result
