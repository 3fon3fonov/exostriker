import sys
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
from pyqtgraph.dockarea import DockArea, Dock

class Accel_GUI():
    def __init__(self, window, dock_area):
        self.testing = 0
        self.pen = pg.mkPen(color='y')
        """Setup the UI"""
        self.window = window
        self.dock_area = dock_area
        self.window.setCentralWidget(self.dock_area)
        self.spec_dock = []
        self.spec_dock.append(Dock("Spectro 1",
                                   size=(1200, 600),
                                  autoOrientation=False))
        self.spec_dock.append(Dock("Spectro 2",
                                   size=(1200, 600),
                               autoOrientation=False))
        self.dock_area.addDock(self.spec_dock[0], "top")
        self.dock_area.addDock(self.spec_dock[1], "below", self.spec_dock[0])

if __name__ == "__main__":
    app = QtGui.QApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)
    win = QtGui.QMainWindow()
    area = DockArea()
    pyqtplot = Accel_GUI(win, area)
    win.show()
    app.exec_()