"""Sampler parameter selectors with Qt-rendered mathematical labels."""
from PyQt6 import QtCore, QtGui, QtWidgets
from sample_plot_labels import plot_label

RICH_LABEL_ROLE = int(QtCore.Qt.ItemDataRole.UserRole) + 1


def document(text, font, color):
    doc = QtGui.QTextDocument()
    doc.setDocumentMargin(0)
    doc.setDefaultFont(font)
    doc.setDefaultStyleSheet('body { color: %s; }' % color.name())
    doc.setHtml('<body>' + text + '</body>')
    return doc


class SampleLabelDelegate(QtWidgets.QStyledItemDelegate):
    def paint(self, painter, option, index):
        text = index.data(RICH_LABEL_ROLE)
        if text is None:
            return super().paint(painter, option, index)
        opt = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        selected = bool(opt.state & QtWidgets.QStyle.StateFlag.State_Selected)
        role = QtGui.QPalette.ColorRole.HighlightedText if selected else QtGui.QPalette.ColorRole.Text
        doc = document(text, opt.font, opt.palette.color(role))
        opt.text = ''
        style = opt.widget.style() if opt.widget else QtWidgets.QApplication.style()
        style.drawControl(QtWidgets.QStyle.ControlElement.CE_ItemViewItem, opt, painter, opt.widget)
        rect = style.subElementRect(QtWidgets.QStyle.SubElement.SE_ItemViewItemText, opt, opt.widget)
        painter.save()
        painter.setClipRect(rect)
        painter.translate(rect.left(), rect.top() + max(0, (rect.height()-doc.size().height())/2))
        doc.drawContents(painter)
        painter.restore()

    def sizeHint(self, option, index):
        text = index.data(RICH_LABEL_ROLE)
        if text is None:
            return super().sizeHint(option, index)
        doc = document(text, option.font, option.palette.color(QtGui.QPalette.ColorRole.Text))
        return QtCore.QSize(int(doc.idealWidth())+12, int(doc.size().height())+6)


class SampleLabelComboBox(QtWidgets.QComboBox):
    """Render both the selected item and popup with real subscripts/Greek."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setItemDelegate(SampleLabelDelegate(self))

    def addMathItem(self, label, userData=None):
        rich = plot_label(label)
        doc = QtGui.QTextDocument()
        doc.setHtml(rich)
        # Keep DisplayRole plain for accessibility and Qt's item searching.
        self.addItem(doc.toPlainText(), userData)
        self.setItemData(self.count()-1, rich, RICH_LABEL_ROLE)
        self.setItemData(self.count()-1, '<span>'+rich+'</span>', QtCore.Qt.ItemDataRole.ToolTipRole)

    def paintEvent(self, event):
        text = self.currentData(RICH_LABEL_ROLE)
        if text is None:
            return super().paintEvent(event)
        painter = QtWidgets.QStylePainter(self)
        opt = QtWidgets.QStyleOptionComboBox()
        self.initStyleOption(opt)
        painter.drawComplexControl(QtWidgets.QStyle.ComplexControl.CC_ComboBox, opt)
        rect = self.style().subControlRect(QtWidgets.QStyle.ComplexControl.CC_ComboBox, opt,
                                           QtWidgets.QStyle.SubControl.SC_ComboBoxEditField, self)
        doc = document(text, self.font(), opt.palette.color(QtGui.QPalette.ColorRole.ButtonText))
        painter.setClipRect(rect)
        painter.translate(rect.left(), rect.top()+max(0,(rect.height()-doc.size().height())/2))
        doc.drawContents(painter)
