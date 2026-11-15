"""Optional Linux X11 terminal embedded in the Designer-owned xterm page."""
import os
import shutil
import sys
import ctypes
import ctypes.util

from PyQt6 import QtCore, QtGui, QtWidgets

from xterm_readme_ui import Ui_XTermReadme


class X11Geometry:
    """Size xterm's foreign X11 child to the native Qt host, in pixels."""

    def __init__(self):
        library = ctypes.util.find_library('X11')
        if not library:
            raise OSError('X11 library is unavailable')
        self.library = ctypes.CDLL(library)
        display = ctypes.c_void_p
        window = ctypes.c_ulong
        window_pointer = ctypes.POINTER(window)
        integer_pointer = ctypes.POINTER(ctypes.c_uint)
        self.library.XOpenDisplay.argtypes = [ctypes.c_char_p]
        self.library.XOpenDisplay.restype = display
        self.library.XCloseDisplay.argtypes = [display]
        self.library.XQueryTree.argtypes = [
            display, window, window_pointer, window_pointer,
            ctypes.POINTER(window_pointer), integer_pointer]
        self.library.XGetGeometry.argtypes = [
            display, window, window_pointer,
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
            integer_pointer, integer_pointer, integer_pointer, integer_pointer]
        self.library.XMoveResizeWindow.argtypes = [
            display, window, ctypes.c_int, ctypes.c_int, ctypes.c_uint, ctypes.c_uint]
        self.library.XSetWindowBorderWidth.argtypes = [display, window, ctypes.c_uint]
        self.library.XFree.argtypes = [ctypes.c_void_p]
        self.library.XGrabServer.argtypes = [display]
        self.library.XUngrabServer.argtypes = [display]
        self.library.XSync.argtypes = [display, ctypes.c_int]
        self.library.XFlush.argtypes = [display]
        self.display = self.library.XOpenDisplay(None)
        if not self.display:
            raise OSError('Cannot connect to the X11 display')

    def fit(self, host_id):
        """Fill the host using its physical size, including high-DPI scaling."""
        root = ctypes.c_ulong()
        parent = ctypes.c_ulong()
        children = ctypes.POINTER(ctypes.c_ulong)()
        count = ctypes.c_uint()
        x = ctypes.c_int()
        y = ctypes.c_int()
        width = ctypes.c_uint()
        height = ctypes.c_uint()
        border = ctypes.c_uint()
        depth = ctypes.c_uint()
        # Keep another X client from destroying a child between querying and
        # resizing it. Always release the server, including on Python errors.
        self.library.XGrabServer(self.display)
        try:
            geometry = self.library.XGetGeometry(
                self.display, host_id, ctypes.byref(root), ctypes.byref(x), ctypes.byref(y),
                ctypes.byref(width), ctypes.byref(height), ctypes.byref(border), ctypes.byref(depth))
            if not geometry:
                return False
            found = self.library.XQueryTree(
                self.display, host_id, ctypes.byref(root), ctypes.byref(parent),
                ctypes.byref(children), ctypes.byref(count))
            if not found or not count.value:
                return False
            for index in range(count.value):
                self.library.XSetWindowBorderWidth(self.display, children[index], 0)
                self.library.XMoveResizeWindow(
                    self.display, children[index], 0, 0,
                    max(1, width.value), max(1, height.value))
            self.library.XSync(self.display, 0)
            return True
        finally:
            if children:
                self.library.XFree(children)
            self.library.XUngrabServer(self.display)
            self.library.XFlush(self.display)

    def close(self):
        """Release only the display connection owned by the resize helper."""
        if self.display:
            self.library.XCloseDisplay(self.display)
            self.display = None


def available_xterm():
    """Return xterm only where its X11 embedding can work."""
    if not sys.platform.startswith('linux'):
        return None
    if QtGui.QGuiApplication.platformName() != 'xcb' or not os.environ.get('DISPLAY'):
        return None
    return shutil.which('xterm')


class EmbeddedXTerm(QtCore.QObject):
    """Run a real xterm inside the optional Linux tab, with its own PTY.

    The terminal retains native keyboard handling, shell configuration, job
    control and fullscreen applications. Start only when selected; hide the
    page when xterm/X11 is unavailable or the terminal cannot launch.
    """

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.readme_dialog = None
        window.xterm_readme_button.clicked.connect(self.show_readme)
        self.executable = available_xterm()
        self.process = QtCore.QProcess(self)
        self.process.errorOccurred.connect(self.failed)
        self.process.finished.connect(self.finished)
        self.process.readyReadStandardError.connect(self.read_errors)
        self.process.started.connect(self.schedule_resize)
        self.geometry = None
        self.resize_attempts = 0
        self.resize_timer = QtCore.QTimer(self)
        self.resize_timer.setInterval(50)
        self.resize_timer.timeout.connect(self.resize_terminal)
        window.xterm_container.installEventFilter(self)
        self.stopping = False
        self.errors = ''
        self.tabs = window.terminal_embeded
        self.page = window.xterm_tab
        self.tabs.setTabVisible(self.tabs.indexOf(self.page), bool(self.executable))
        self.tabs.currentChanged.connect(self.activate)
        QtCore.QCoreApplication.instance().aboutToQuit.connect(self.shutdown)

    def show_readme(self):
        """Show the Designer-owned terminal instructions without blocking the shell."""
        if self.readme_dialog is None:
            self.readme_dialog = QtWidgets.QDialog(self.window)
            instructions = Ui_XTermReadme()
            instructions.setupUi(self.readme_dialog)
        self.readme_dialog.show()
        self.readme_dialog.raise_()
        self.readme_dialog.activateWindow()

    def activate(self, index):
        """Embed xterm into the native host window on the first tab selection."""
        if not self.executable or self.tabs.widget(index) is not self.page:
            return
        if self.process.state() != QtCore.QProcess.ProcessState.NotRunning:
            self.schedule_resize()
            return
        self.stopping = False
        self.errors = ''
        host_id = str(int(self.window.xterm_container.winId()))
        self.process.start(self.executable, ['-into', host_id])

    def eventFilter(self, watched, event):
        """Follow changes to the existing Designer host, without adding widgets."""
        if event.type() in (QtCore.QEvent.Type.Resize, QtCore.QEvent.Type.Show):
            self.schedule_resize()
        return super().eventFilter(watched, event)

    def schedule_resize(self):
        """Wait for native geometry updates and xterm's asynchronous creation."""
        if self.process.state() == QtCore.QProcess.ProcessState.Running:
            self.resize_attempts = 40
            self.resize_timer.start()

    def resize_terminal(self):
        """Match the embedded terminal to the whole available tab area."""
        self.resize_attempts -= 1
        if self.process.state() != QtCore.QProcess.ProcessState.Running:
            self.resize_timer.stop()
            return
        try:
            if self.geometry is None:
                self.geometry = X11Geometry()
            self.geometry.fit(int(self.window.xterm_container.winId()))
        except OSError as error:
            self.resize_timer.stop()
            self.shutdown()
            self.hide_unavailable(str(error))
            return
        # xterm can create its child before completing font/PTY setup, then
        # restore its default geometry after our first successful resize. Keep
        # fitting briefly while startup and the surrounding layout settle.
        if self.resize_attempts <= 0:
            self.resize_timer.stop()

    def read_errors(self):
        """Keep launch diagnostics without filling an unread process pipe."""
        self.errors += bytes(self.process.readAllStandardError()).decode('utf-8', errors='replace')
        self.errors = self.errors[-4000:]

    def failed(self, error):
        """An unavailable optional terminal must not interrupt normal startup."""
        if not self.stopping:
            self.hide_unavailable(self.process.errorString())

    def hide_unavailable(self, reason):
        """Return to Jupyter when embedding fails, leaving the terminal hidden."""
        self.tabs.setCurrentWidget(self.window.console_widget)
        self.tabs.setTabVisible(self.tabs.indexOf(self.page), False)
        print(f'Optional xterm unavailable: {reason.strip()}')

    def finished(self, code, status):
        """Let a normal shell exit be reopened; hide a failed terminal."""
        self.read_errors()
        self.resize_timer.stop()
        if not self.stopping and code != 0:
            self.hide_unavailable(self.errors or f'xterm exited with status {code}')
        elif not self.stopping:
            self.tabs.setCurrentWidget(self.window.console_widget)

    def shutdown(self):
        """Close xterm so it closes its PTY and hangs up the interactive shell."""
        self.stopping = True
        self.resize_timer.stop()
        if self.geometry is not None:
            self.geometry.close()
            self.geometry = None
        if self.process.state() == QtCore.QProcess.ProcessState.NotRunning:
            return
        self.process.terminate()
        if not self.process.waitForFinished(1500):
            self.process.kill()
            self.process.waitForFinished(1000)
