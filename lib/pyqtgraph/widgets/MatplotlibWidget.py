from ..Qt import QtGui, QtCore, QT_LIB
import matplotlib

########################### For nice plotting ##################################

#matplotlib.rcParams['axes.linewidth'] = 2.0 #set the value globally
#matplotlib.rcParams['xtick.major.pad']='1'
#matplotlib.rcParams['ytick.major.pad']='2'


# set tick width
#matplotlib.rcParams['xtick.major.size'] = 8
#matplotlib.rcParams['xtick.major.width'] = 2
#matplotlib.rcParams['xtick.minor.size'] = 5
#matplotlib.rcParams['xtick.minor.width'] = 2

#matplotlib.rcParams['ytick.major.size'] = 8
#matplotlib.rcParams['ytick.major.width'] = 2
#matplotlib.rcParams['ytick.minor.size'] = 5
#matplotlib.rcParams['ytick.minor.width'] = 2

matplotlib.rcParams['axes.formatter.useoffset'] = False


#rc('text',usetex=True)
#font = {'family' : 'normal','weight' : 'black','size': 22,'serif':['Helvetica']}
#rc('font', **font)


if QT_LIB != 'PyQt5':
    if QT_LIB == 'PySide':
        matplotlib.rcParams['backend.qt4']='PySide'

    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    try:
        from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
    except ImportError:
        from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
else:
   # from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
   # from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    from .backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from .backend_qt5agg import NavigationToolbar2QT as NavigationToolbar    
    
    #from matplotlib.backends.qt_editor import figureoptions as figureoptions #https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/backends/qt_editor/figureoptions.py
from matplotlib.figure import Figure 

#print(figureoptions)

#class NavigationToolbar(NavigationToolbar):
#    # only display the buttons we need
#    toolitems = [t for t in NavigationToolbar.toolitems if
#                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]






class MatplotlibWidget(QtGui.QWidget):
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
        QtGui.QWidget.__init__(self)
        self.fig = Figure(size, dpi=dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        self.vbox = QtGui.QVBoxLayout()
        self.vbox.addWidget(self.toolbar)
        self.vbox.addWidget(self.canvas)
        
        self.setLayout(self.vbox)
        

    def getFigure(self):
        return self.fig
        
    def draw(self):
        self.canvas.draw()
