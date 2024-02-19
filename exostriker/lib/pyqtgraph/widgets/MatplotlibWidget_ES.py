from ..Qt import QtWidgets, QtCore, QT_LIB
import matplotlib

 

if QT_LIB != 'PyQt5':
    if QT_LIB == 'PySide':
        matplotlib.rcParams['backend.qt4']='PySide'

    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    try:
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QTAgg as NavigationToolbar
    except ImportError:
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
else:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure


# matplotlib.rcParams.keys()

matplotlib.rcParams['axes.formatter.useoffset'] = False
matplotlib.rcParams['axes.formatter.limits'] = [-9, 9]
#matplotlib.rcParams['axes.grid'] = True

 
class MatplotlibWidget(QtWidgets.QWidget):
    """
    Implements a Matplotlib figure inside a QWidget.
    Use getFigure() and redraw() to interact with matplotlib.
    
    Example::
    
        mw = MatplotlibWidget()
        subplot = mw.getFigure().add_subplot(111)
        subplot.plot(x,y)
        mw.draw()
    """
    
    def __init__(self, size=(8.0, 6.0), dpi=125):
        QtWidgets.QWidget.__init__(self)
        self.fig = Figure(size, dpi=dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        self.vbox = QtWidgets.QVBoxLayout()
        self.vbox.addWidget(self.toolbar)
        self.vbox.addWidget(self.canvas)
        
        self.setLayout(self.vbox)

    def getFigure(self):
        return self.fig
        
    def draw(self):
        self.canvas.draw()
