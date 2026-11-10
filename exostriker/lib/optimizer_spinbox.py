"""Scientific numeric input promoted from QDoubleSpinBox in Designer."""
from PyQt6 import QtCore, QtGui, QtWidgets


class OptimizerDoubleSpinBox(QtWidgets.QDoubleSpinBox):
    """Compact decimal/scientific input for small convergence tolerances."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setLocale(QtCore.QLocale.c())
        self.number_validator = QtGui.QDoubleValidator(0., 1e8, 12, self)
        self.number_validator.setLocale(QtCore.QLocale.c())
        self.number_validator.setNotation(QtGui.QDoubleValidator.Notation.ScientificNotation)

    def textFromValue(self, value):
        return format(value, '.12g')

    def valueFromText(self, text):
        value, valid = self.locale().toDouble(text)
        return value if valid else self.minimum()

    def validate(self, text, position):
        if not hasattr(self, 'number_validator'):
            return super().validate(text, position)
        if self.specialValueText() and text == self.specialValueText():
            return QtGui.QValidator.State.Acceptable, text, position
        self.number_validator.setRange(self.minimum(), self.maximum(), self.decimals())
        return self.number_validator.validate(text, position)
