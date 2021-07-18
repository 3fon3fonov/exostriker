# -*- coding: utf-8 -*-
from ..Qt import QtGui, QtCore
from .Exporter import Exporter
from ..parametertree import Parameter
from .. import PlotItem
from ..python2_3 import asUnicode

__all__ = ['CSVExporter']
    
translate = QtCore.QCoreApplication.translate

class CSVExporter(Exporter):
    Name = "CSV from plot data"
    windows = []
    def __init__(self, item):
        Exporter.__init__(self, item)
        self.params = Parameter(name='params', type='group', children=[
            {'name': 'separator', 'title': translate("Exporter", 'separator'), 'type': 'list', 'value': 'comma', 'values': ['comma', 'tab']},
            {'name': 'precision', 'title': translate("Exporter", 'precision'), 'type': 'int', 'value': 10, 'limits': [0, None]},
            {'name': 'columnMode', 'title': translate("Exporter", 'columnMode'), 'type': 'list', 'values': ['(x,y,err_y,err_x) per plot', '(x,y,y,y) for all plots']}
        ])
        
        
    def parameters(self):
        return self.params
    
    def export(self, fileName=None):
        
        if not isinstance(self.item, PlotItem):
            raise Exception("Must have a PlotItem selected for CSV export.")
        
        if fileName is None:
            self.fileSaveDialog(filter=["*.csv", "*.tsv"])
            return

        data = []
        header = []

        appendAllX = self.params['columnMode'] == '(x,y,err_y,err_x) per plot'

        for i, c in enumerate(self.item.curves):
            
            cd = c.getData()
            if cd[0] is None:
                continue
            data.append(cd)

            if hasattr(c, 'implements') and c.implements('plotData') and c.name() is not None:
                name = c.name().replace('"', '""') + '_'
                xName, yName = '"'+name+'x"', '"'+name+'y"'
            else:
                xName = 'x%04d' % i
                yName = 'y%04d' % i
                e_yName = 'e_y%04d' % i
                e_xName = 'e_x%04d' % i

            if appendAllX or i == 0:
                header.extend([xName, yName,e_yName,e_xName ])
            else:
                header.extend([yName,e_yName,e_xName ])
 
    
     
        for indx, cc in enumerate(self.item.items):
 
            if cc.__class__.__name__ == "ErrorBarItem":                               
            
                for option in ['top', 'left']:
                    if option in cc.opts and i > 0:
                        err = cc.opts[option]
                        try:
                            data[1] = list(data[1])
                            data[1].append(err)
                            data[1] = tuple(data[1])
                        except:
                            continue
                    else:
                        continue

        if self.params['separator'] == 'comma':
            sep = ','
        else:
            sep = '\t'

        with open(fileName, 'w') as fd:
            fd.write(sep.join(map(asUnicode, header)) + '\n')
            i = 0
            numFormat = '%%0.%dg' % self.params['precision']
            numRows = max([len(d[0]) for d in data])
            for i in range(numRows):
                for j, d in enumerate(data):
                    # write x value if this is the first column, or if we want
                    # x for all rows
                    if appendAllX or j == 0:
                        if d is not None and i < len(d[0]):
                            fd.write(numFormat % d[0][i] + sep)
                        else:
                            fd.write(' %s' % sep)

                    # write y value
                    if d is not None and i < len(d[1]):
                        fd.write(numFormat % d[1][i] + sep)
                    else:
                        fd.write(' %s' % sep)

                    # write data e_y value
                    if d is not None and len(d)>2 and d[2] is not None and i < len(d[2]):
                            fd.write(numFormat % d[2][i] + sep)
                    else:
                        fd.write(' %s' % sep)

                    # write data e_x value
                    if d is not None and len(d)>3 and d[3] is not None and i < len(d[3]):
                            fd.write(numFormat % d[3][i] + sep)
                    else:
                        fd.write(' %s' % sep)

                fd.write('\n')


CSVExporter.register()        
                
        
