"""Non-blocking, optional PyPI update notification for the desktop app."""
import json
import queue
import threading
import time
from urllib.error import URLError
from urllib.request import Request, urlopen

from packaging.version import InvalidVersion, Version
from PyQt6 import QtCore, QtWidgets


def newer_version(payload, current_version):
    """Return the advertised stable release, or None for unusable metadata."""
    try:
        info = payload['info']
        latest = Version(info['version'])
        if (not info.get('yanked', False) and not latest.is_prerelease
                and not latest.is_devrelease and latest > Version(current_version)):
            return str(latest)
    except (KeyError, TypeError, InvalidVersion):
        pass
    return None


def fetch_latest_version(current_version):
    """Use Python's verified HTTPS, independent of Qt's OpenSSL backend."""
    request = Request('https://pypi.org/pypi/exostriker/json',
                      headers={'Accept': 'application/json'})
    try:
        with urlopen(request, timeout=5) as response:
            return newer_version(json.load(response), current_version)
    except (OSError, URLError, ValueError, TypeError, UnicodeError):
        return None


def _fetch_in_background(results, current_version):
    # No Qt objects are accessed from the daemon thread, including at shutdown.
    results.put(fetch_latest_version(current_version))


class UpdateCheck(QtCore.QObject):
    """Keep network activity and the modeless notification owned by the window."""

    def __init__(self, window, current_version, settings=None):
        super().__init__(window)
        self.window = window
        self.current_version = current_version
        self.settings = settings if settings is not None else QtCore.QSettings(
            'Exostriker', 'Exostriker')
        self.results = queue.Queue()
        self.started = False
        self.dialog = None
        self.poll = QtCore.QTimer(self)
        self.poll.setInterval(100)
        self.poll.timeout.connect(self._poll_result)

    def start(self):
        if self.started:
            return
        self.started = True
        self.deadline = time.monotonic() + 6
        threading.Thread(target=_fetch_in_background,
                         args=(self.results, self.current_version), daemon=True).start()
        self.poll.start()

    def _poll_result(self):
        # Bound the UI's wait even if the OS takes longer to resolve DNS.
        if time.monotonic() >= self.deadline:
            self.poll.stop()
            return
        try:
            latest = self.results.get_nowait()
        except queue.Empty:
            return
        self.poll.stop()
        ignored = self.settings.value('updates/ignored_version', '', type=str)
        if latest and latest != ignored and self.window.isVisible():
            self._show_update(latest)

    def _show_update(self, latest):
        dialog = QtWidgets.QMessageBox(self.window)
        dialog.setWindowTitle('Exo-Striker update available')
        dialog.setIcon(QtWidgets.QMessageBox.Icon.Information)
        dialog.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        dialog.setText('A newer version of the Exo-Striker is available on PyPI.')
        dialog.setInformativeText(
            f'Current version: {self.current_version}\n'
            f'Available version: {latest}\n\n'
            'To update, close the Exo-Striker and run this command in a terminal '
            'in the Python environment where the Exo-Striker is installed:\n\n'
            'pip install exostriker -U\n\n'
            'Restart the Exo-Striker after the update.')
        dialog.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        checkbox = QtWidgets.QCheckBox(f"Don't show again for version {latest}", dialog)
        dialog.setCheckBox(checkbox)
        dialog.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        checkbox.toggled.connect(lambda checked: self._set_ignored_version(latest, checked))
        self.dialog = dialog
        dialog.open()

    def _set_ignored_version(self, version, ignored):
        if ignored:
            self.settings.setValue('updates/ignored_version', version)
        else:
            self.settings.remove('updates/ignored_version')
        # Retire the former permanent opt-out; future releases must be checked.
        self.settings.remove('updates/disabled')
        self.settings.sync()
