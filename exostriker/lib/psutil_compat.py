"""Guard the supervisor against broken Linux psutil PID-check binaries.

Some system abi3 builds return a borrowed None from check_pid_range when
loaded by an older CPython. Process.children repeatedly constructs Process
objects, so even an idle calculation eventually aborts with none_dealloc.
Only replace that validator when its reference-count defect is observed;
process discovery, PID-reuse checks and cancellation remain psutil's own.
"""
import operator
import os
import sys


def _check_linux_pid_range(value):
    """Equivalent range validation for Linux's signed 32-bit pid_t."""
    value = operator.index(value)
    if value < -(2**31):
        raise OverflowError('signed integer is less than minimum')
    if value >= 2**31:
        raise OverflowError('signed integer is greater than maximum')
    if value < 0:
        raise ValueError('pid must be a positive integer')


def protect_pid_checks(psutil):
    """Return True only when the defective native validator was replaced."""
    if (sys.implementation.name != 'cpython' or sys.version_info >= (3, 12)
            or not sys.platform.startswith('linux')):
        return False
    extension = getattr(getattr(psutil, '_psplatform', None), 'cext', None)
    check = getattr(extension, 'check_pid_range', None)
    if check is None or check is _check_linux_pid_range:
        return False
    pid = os.getpid()
    # Warm lazy state, then use a tiny bounded probe. A correct validator
    # leaves None's count unchanged; the faulty binary loses one per call.
    check(pid)
    before = sys.getrefcount(None)
    for _ in range(8):
        check(pid)
    lost = before - sys.getrefcount(None)
    if lost >= 8:
        extension.check_pid_range = _check_linux_pid_range
        return True
    return False
