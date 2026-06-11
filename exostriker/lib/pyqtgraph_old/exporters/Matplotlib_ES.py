#from ..Qt import QtGui, QtCore
from .Exporter import Exporter
from .. import PlotItem
from .. import functions as fn
from ..Qt import QtCore, QtWidgets

__all__ = ['MatplotlibExporter']

"""
It is helpful when using the matplotlib Exporter if your
.matplotlib/matplotlibrc file is configured appropriately.
The following are suggested for getting usable PDF output that
can be edited in Illustrator, etc.

backend      : Qt4Agg
text.usetex : True  # Assumes you have a findable LaTeX installation
interactive : False
font.family : sans-serif
font.sans-serif : 'Arial'  # (make first in list)
mathtext.default : sf
figure.facecolor : white  # personal preference
# next setting allows pdf font to be readable in Adobe Illustrator
pdf.fonttype : 42   # set fonts to TrueType (otherwise it will be 3
                    # and the text will be vectorized.
text.dvipnghack : True  # primarily to clean up font appearance on Mac

The advantage is that there is less to do to get an exported file cleaned and ready for
publication. Fonts are not vectorized (outlined), and window colors are white.

"""
    
class MatplotlibExporter(Exporter):
    Name = "Matplotlib Window"
    windows = []
    def __init__(self, item):
        Exporter.__init__(self, item)
        
    def parameters(self):
        return None

    def cleanAxes(self, axl):
        if type(axl) is not list:
            axl = [axl]
        for ax in axl:
            if ax is None:
                continue
            for loc, spine in ax.spines.items():
                if loc in ['left', 'bottom']:
                    pass
                elif loc in ['right', 'top']:
                    spine.set_color('none')
                    # do not draw the spine
                else:
                    raise ValueError('Unknown spine location: %s' % loc)
                # turn off ticks when there is no spine
                ax.xaxis.set_ticks_position('bottom')
    
    def export(self, fileName=None):
        
        if isinstance(self.item, PlotItem):
            mpw = MatplotlibWindow()
            MatplotlibExporter.windows.append(mpw)

            stdFont = 'Arial'
            
            fig = mpw.getFigure()
            

            # get labels from the graphic item
            xlabel = self.item.axes['bottom']['item'].label.toPlainText()
            ylabel = self.item.axes['left']['item'].label.toPlainText()
            title = self.item.titleLabel.text

            ax = fig.add_subplot(111, title=title)
            ax.clear()
            self.cleanAxes(ax)
            #ax.grid(True)
 
            for indx, item in enumerate(self.item.curves):
                x, y = item.getData() 
                
                
                if self.item.curves[indx].__class__.__name__ == "PlotDataItem" and self.item.curves[indx].opts['logMode'][0] == True:
                    gls_log = True
                else:
                    gls_log = False
                    

                if x is None:
                    continue
                else:
                    if gls_log == True:
                        x = 10.0**x
                       
 
                opts = item.opts
                #print(self.item.curves[indx].__class__.__name__)

                pen = fn.mkPen(opts['pen'])
                if pen.style() == QtCore.Qt.NoPen:
                    linestyle = ''
                    zorder = 10
                else:
                    linestyle = '-'
                    if len(y) < 500:
                        zorder = -100
                    else:
                        zorder = 10
 
                    
                color = tuple([c/255. for c in fn.colorTuple(pen.color())])
                symbol = opts['symbol']
                if symbol == 't':
                    symbol = '^'
                symbolPen = fn.mkPen(opts['symbolPen'])
                symbolBrush = fn.mkBrush(opts['symbolBrush'])
                markeredgecolor = tuple([c/255. for c in fn.colorTuple(symbolPen.color())])
                markerfacecolor = tuple([c/255. for c in fn.colorTuple(symbolBrush.color())])
                markersize = opts['symbolSize']
             
                if opts['fillLevel'] is not None and opts['fillBrush'] is not None:
                    fillBrush = fn.mkBrush(opts['fillBrush'])
                    fillcolor = tuple([c/255. for c in fn.colorTuple(fillBrush.color())])
                    ax.fill_between(x=x, y1=y, y2=opts['fillLevel'], facecolor=fillcolor)
                
                if symbol == 't':
                    symbol = 'v'
                elif symbol == 't1':
                    symbol = '^'
                elif symbol == 't2':
                    symbol = '>'
                elif symbol == 't3':
                    symbol = '<'
                elif symbol == 'star':
                    symbol = '*'
                    
 
                ax.plot(x, y, marker=symbol, color=color, linewidth=pen.width(), 
                        linestyle=linestyle, 
                        markeredgecolor=markeredgecolor, 
                        markerfacecolor=markerfacecolor,
                        markersize=markersize, zorder = zorder)

 
                xr, yr = self.item.viewRange()

                if xr is None:
                    continue
                else:
                    if gls_log == True:
                        xr[0] = 10.0**xr[0]
                        xr[1] = 10.0**xr[1]                       
                        ax.semilogx()
                        from matplotlib import ticker
                        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())  # set regular formatting
                        ax.ticklabel_format(style='sci', scilimits=(-6, 9))  # disable scientific notation
                        
                ax.set_xbound(*xr)
                ax.set_ybound(*yr)
                
                
            for indx, item in enumerate(self.item.items):
                #if indx <=1:
               #     continue
                

                if self.item.items[indx].__class__.__name__ == "InfiniteLine":
                    level = self.item.items[indx].value() 
                    Inf_line_color = self.item.items[indx].pen.color().name()
                    #pen = fn.mkPen(self.item.items[indx].pen)
                    #print(pen)
                    if Inf_line_color == '#000000': 
                        linestyle="--" 
                    else: 
                        linestyle="-"
                    ax.axhline(y=level, linewidth=0.8, linestyle=linestyle, color=Inf_line_color, zorder=-30)
                    
                    continue   
               
                if self.item.items[indx].__class__.__name__ == "FillBetweenItem":
                    x1,y1 = self.item.items[indx].curves[0].getData()
                    x2,y2 = self.item.items[indx].curves[1].getData()
                    #fillBrush = fn.mkBrush(self.item.items[indx].curves[0].opts['fillBrush'])
                   # fillcolor = tuple([c/255. for c in fn.colorTuple(fillBrush.color())])                   
                    ax.fill_between(x=x1, y1=y1, y2=y2, facecolor="#f48c42" )#fn.mkColor(244,140,66,128))
                    #print(self.item.items[indx].curve1)                               
                    continue

                if self.item.items[indx].__class__.__name__ == "TextItem":
                    continue
                

                
                if "top" in self.item.items[indx].opts:
                    yerr=self.item.items[indx].opts["top"]
                    ax.errorbar(x=self.item.items[indx].opts["x"], y=self.item.items[indx].opts["y"],
                                     yerr=yerr, linestyle='', marker='o',markersize=0.5, linewidth=1.0, 
                                     color=self.item.items[indx].opts["pen"], capsize = 0, elinewidth=1,mew=0.0, zorder=-10)
                                     #color="k", capsize = 0, elinewidth=1,mew=0.0, zorder=-10)
                if "left" in self.item.items[indx].opts:
                    xerr=self.item.items[indx].opts["left"]
                    ax.errorbar(x=self.item.items[indx].opts["x"], y=self.item.items[indx].opts["y"],
                                     xerr=xerr, linestyle='', marker='o',markersize=0.5, linewidth=1.0, 
                                     color=self.item.items[indx].opts["pen"], capsize = 0, elinewidth=1,mew=0.0, zorder=-10)  


                                         
            #print(self.item.items[indx].opts)    
            ax.set_xlabel(xlabel)  # place the labels.
            ax.set_ylabel(ylabel)
            #ax.spines['top'].set_visible(True)
            #ax.spines['right'].set_visible(True)

            ax.spines['top'].set_color('k')
            ax.spines['right'].set_color('k')


            mpw.draw()
        else:
          #  raise Exception("Matplotlib export currently only works with plot items")
            print("Matplotlib export currently only works with plot items, and not with scenes")             
   
MatplotlibExporter.register()        
        

class MatplotlibWindow(QtWidgets.QMainWindow):
    def __init__(self):
        from ..widgets import MatplotlibWidget_ES
        QtWidgets.QMainWindow.__init__(self)
        self.mpl = MatplotlibWidget_ES.MatplotlibWidget()
        self.setCentralWidget(self.mpl)
        self.show()
        
    def __getattr__(self, attr):
        return getattr(self.mpl, attr)
        
    def closeEvent(self, ev):
        MatplotlibExporter.windows.remove(self)


