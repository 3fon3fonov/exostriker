__all__ = ['PDFExporter']

 
from ..Qt import QtCore,QtWidgets
from ..Qt.QtWidgets import QGraphicsItem, QApplication
from ..Qt.QtGui import QPainter, QPdfWriter, QPageSize
from ..Qt.QtCore import QMarginsF, Qt, QSizeF, QRectF
from .Exporter import Exporter
from ..parametertree import Parameter

translate = QtCore.QCoreApplication.translate

class PDFExporter(Exporter):
    """A pdf exporter for pyqtgraph graphs. Based on pyqtgraph's
     ImageExporter.

     There is a bug in Qt<5.12 that makes Qt wrongly use a cosmetic pen
     (QTBUG-68537). Workaround: do not use completely opaque colors.

     There is also a bug in Qt<5.12 with bold fonts that then remain bold.
     To see it, save the OWNomogram output."""

    Name = "Scalable Vector Graphics (PDF)"
    allowCopy=True

    def __init__(self, item):
        Exporter.__init__(self, item)

        tr = self.getTargetRect()

        scene = item.scene() if isinstance(item, QtWidgets.QGraphicsItem) else item
        bgbrush = scene.views()[0].backgroundBrush()
        bg = bgbrush.color()

        self.params = Parameter(name='params', type='group', children=[
            {'name': 'background', 'title': translate("Exporter", 'background'), 'type': 'color', 'value': bg},
            {'name': 'width', 'title': translate("Exporter", 'width'), 'type': 'float', 'value': tr.width(),
             'limits': (0, None)},
            {'name': 'height', 'title': translate("Exporter", 'height'), 'type': 'float', 'value': tr.height(),
             'limits': (0, None)},
            #{'name': 'viewbox clipping', 'type': 'bool', 'value': True},
            #{'name': 'normalize coordinates', 'type': 'bool', 'value': True},
     #       {
      #          'name': 'scaling stroke',
      #          'title': translate("Exporter", 'scaling stroke'),
      #          'type': 'bool',
      #          'value': False,
      #          'tip': "If False, strokes are non-scaling, which means that "
      #                 "they appear the same width on screen regardless of "
      #                 "how they are scaled or how the view is zoomed."
      #      },
        ])
        #self.params.param('width').sigValueChanged.connect(self.widthChanged)
        #self.params.param('height').sigValueChanged.connect(self.heightChanged)


        if isinstance(item, QGraphicsItem):
            scene = item.scene()
        else:
            scene = item
        bgbrush = scene.views()[0].backgroundBrush()
        bg = bgbrush.color()
        if bgbrush.style() == Qt.NoBrush:
            bg.setAlpha(0)
        self.background = bg

        # The following code is a workaround for a bug in pyqtgraph 1.1. The suggested
        # fix upstream was pyqtgraph/pyqtgraph#1458
        try:
            from pyqtgraph.graphicsItems.ViewBox.ViewBox import ChildGroup
            for item in self.getPaintItems():
                if isinstance(item, ChildGroup):
                    if item.flags() & QGraphicsItem.ItemClipsChildrenToShape:
                        item.setFlag(QGraphicsItem.ItemClipsChildrenToShape, False)
        except:  # pylint: disable=bare-except
            pass




    def parameters(self):
        return self.params

    def export(self, filename=None):
        pw = QPdfWriter(filename)
        dpi = int(QApplication.primaryScreen().logicalDotsPerInch())
        pw.setResolution(dpi)
        pw.setPageMargins(QMarginsF(0, 0, 0, 0))
#        pw.setPageSize(
#            QPageSize(QSizeF(self.getTargetRect().size()) / dpi * 25.4,
#                      QPageSize.Millimeter))
        painter = QPainter(pw)
        try:
            self.setExportMode(True, {'antialias': True,
                                      'background': self.background,
                                      'painter': painter})
            painter.setRenderHint(QPainter.Antialiasing, True)
            if QtCore.QT_VERSION >= 0x050D00:
                painter.setRenderHint(QPainter.LosslessImageRendering, True)
            self.getScene().render(painter,
                                   QRectF(self.getTargetRect()),
                                   QRectF(self.getSourceRect()))
        finally:
            self.setExportMode(False)

        painter.end()

PDFExporter.register()        


