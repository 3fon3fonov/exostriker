"""Qt supervisor for disposable numerical processes, with transactional results."""
import os
import json
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import traceback

import dill
import job_serialization
import psutil
from PyQt6 import QtCore


class CalculationSignals(QtCore.QObject):
    # Exactly one terminal signal, after process cleanup. Only 'result' may be
    # adopted by the GUI; cancelled/error jobs leave its session untouched.
    finished = QtCore.pyqtSignal(str, object)
    output = QtCore.pyqtSignal(str)


def kill_process_tree(process, known_children=()):
    """Reap only this job's processes, including descendants in new sessions."""
    children = set(known_children)
    failure = None
    try:
        try:
            root = psutil.Process(process.pid)
            # Prevent pools from spawning replacements during enumeration.
            root.suspend()
            children.update(root.children(recursive=True))
        except psutil.NoSuchProcess:
            pass
        pending = list(children)
        suspended = set()
        while pending:
            child = pending.pop()
            if child in suspended:
                continue
            suspended.add(child)
            try:
                child.suspend()
                descendants = child.children()
                children.update(descendants)
                pending.extend(descendants)
            except psutil.NoSuchProcess:
                pass
    except psutil.Error as exc:
        failure = exc
    finally:
        # Even an enumeration/suspension failure must not leave the root or
        # already suspended children alive. Popen owns the root process handle.
        for child in children:
            try:
                child.kill()
            except psutil.NoSuchProcess:
                pass
            except psutil.Error as exc:
                failure = exc
        if os.name != 'nt':
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            except OSError as exc:
                failure = exc
        if process.poll() is None:
            process.kill()
        process.wait(timeout=10)
        _, alive = psutil.wait_procs(list(children), timeout=2)
        for child in alive:
            try:
                if child.is_running() and child.status() != psutil.STATUS_ZOMBIE:
                    failure = RuntimeError(f'Calculation child {child.pid} did not exit')
            except psutil.NoSuchProcess:
                pass
    if failure is not None:
        raise RuntimeError(f'Calculation process cleanup failed: {failure}') from failure


class PreparedCalculation:
    """Own one imported, unused child. Ownership transfers to exactly one job."""
    def __init__(self, runner_path=None):
        self.directory = tempfile.TemporaryDirectory(prefix='exostriker-warm-')
        self.folder = Path(self.directory.name)
        self.reaped = False
        self.log = (self.folder / 'output.log').open('wb')
        try:
            kwargs = ({'creationflags': subprocess.CREATE_NEW_PROCESS_GROUP}
                      if os.name == 'nt' else {'start_new_session': True})
            self.process = subprocess.Popen(
                [sys.executable, '-u', str(runner_path or Path(__file__).with_name('calculation_runner.py')),
                 '--warm', str(self.folder)], stdin=subprocess.PIPE, stdout=self.log,
                stderr=subprocess.STDOUT,
                env=dict(os.environ, MPLBACKEND='Agg', PYTHONUNBUFFERED='1',
                         PYTHONIOENCODING='utf-8'), **kwargs)
        except BaseException:
            self.log.close()
            self.directory.cleanup()
            raise

    @property
    def ready(self):
        return self.process.poll() is None and (self.folder / 'ready').exists()

    def close(self):
        try:
            if not self.reaped:
                kill_process_tree(self.process)
                self.reaped = True
        finally:
            self.process.stdin.close()
            self.log.close()
            self.directory.cleanup()


class _RetirePrepared(QtCore.QRunnable):
    def __init__(self, prepared):
        super().__init__()
        self.prepared = prepared

    def run(self):
        try:
            self.prepared.close()
        except Exception:
            traceback.print_exc()


class InitializationWarmup(QtCore.QObject):
    """Keep one disposable import ready while idle, never a pool of sessions."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.prepared = None
        self.paused = False
        self.closing = False
        self.enabled = os.environ.get('EXOSTRIKER_PREWARM_INIT', '1') != '0'
        QtCore.QTimer.singleShot(0, self.prepare)

    @property
    def ready(self):
        return self.prepared is not None and self.prepared.ready

    def prepare(self):
        if self.enabled and not self.paused and not self.closing and self.prepared is None:
            try:
                self.prepared = PreparedCalculation()
            except Exception:
                # Prewarming is optional: ordinary job launch remains available.
                traceback.print_exc()

    def acquire(self):
        self.paused = True
        prepared, self.prepared = self.prepared, None
        if prepared is not None and prepared.process.poll() is not None:
            QtCore.QThreadPool.globalInstance().start(_RetirePrepared(prepared))
            return None
        return prepared

    def pause(self):
        prepared = self.acquire()
        if prepared is not None:
            QtCore.QThreadPool.globalInstance().start(_RetirePrepared(prepared))

    def resume(self):
        self.paused = False
        QtCore.QTimer.singleShot(0, self.prepare)

    def shutdown(self):
        self.closing = True
        self.pause()


class ProcessWorker(QtCore.QRunnable):
    def __init__(self, fit, operation, options):
        super().__init__()
        # Snapshot on the GUI thread, before starting any background work.
        self.request = job_serialization.dumps((fit, operation, options))
        self.prepared = None
        self.signals = CalculationSignals()
        self.stop_event = threading.Event()
        self.process = None
        self.runner_path = Path(__file__).with_name('calculation_runner.py')

    def stop(self):
        # Never block the GUI waiting for a native call or subprocess to exit.
        self.stop_event.set()

    @QtCore.pyqtSlot()
    def run(self):
        status, payload = 'cancelled', None
        known_children = set()
        try:
            with tempfile.TemporaryDirectory(prefix='exostriker-job-') as directory:
                folder = self.prepared.folder if self.prepared is not None else Path(directory)
                request, response = folder / 'request.dill', folder / 'response.dill'
                request.write_bytes(self.request)
                self.request = None
                if not self.stop_event.is_set():
                    runner = self.runner_path
                    env = dict(os.environ, MPLBACKEND='Agg', PYTHONUNBUFFERED='1',
                               PYTHONIOENCODING='utf-8')
                    log_path = folder / 'output.log'
                    with log_path.open('ab') as log, log_path.open('rb') as reader:
                        kwargs = ({'creationflags': subprocess.CREATE_NEW_PROCESS_GROUP}
                                  if os.name == 'nt' else {'start_new_session': True})
                        if self.prepared is not None:
                            self.process = self.prepared.process
                            self.process.stdin.write((json.dumps({'cwd': os.getcwd()}) + '\n').encode())
                            self.process.stdin.flush()
                        else:
                            self.process = subprocess.Popen(
                                [sys.executable, '-u', str(runner), str(request), str(response)],
                                stdin=subprocess.PIPE, stdout=log, stderr=subprocess.STDOUT,
                                env=env, **kwargs)
                        try:
                            root = psutil.Process(self.process.pid)
                            while True:
                                chunk = reader.read(65536)
                                if chunk:
                                    self.signals.output.emit(chunk.decode('utf-8', errors='replace'))
                                try:
                                    known_children.update(root.children(recursive=True))
                                except psutil.NoSuchProcess:
                                    pass
                                if self.stop_event.is_set():
                                    break
                                if response.exists():
                                    with response.open('rb') as stream:
                                        status, payload = dill.load(stream)
                                    if self.stop_event.is_set():
                                        status, payload = 'cancelled', None
                                    break
                                if self.process.poll() is not None:
                                    raise RuntimeError(
                                        f'Calculation process exited with code {self.process.returncode} '
                                        'before returning a result. See Stdout/Stderr.')
                                self.stop_event.wait(0.01 if self.prepared is not None else 0.1)
                        finally:
                            kill_process_tree(self.process, known_children)
                            if self.prepared is not None:
                                self.prepared.reaped = True
                            self.process.stdin.close()
                            chunk = reader.read()
                            if chunk:
                                self.signals.output.emit(chunk.decode('utf-8', errors='replace'))
        except BaseException:
            status, payload = 'error', traceback.format_exc()
        finally:
            if self.prepared is not None:
                try:
                    self.prepared.close()
                except BaseException:
                    status, payload = 'error', traceback.format_exc()
                self.prepared = None
        self.signals.finished.emit(status, payload)
