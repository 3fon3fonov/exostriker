"""Options for the OpenAI API chatbot."""
from PyQt6 import QtWidgets
from exochat_backend import MODELS, read_settings, save_settings


class show_chat_opt(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Exo-Striker Chat Settings')
        self.setMinimumWidth(520)
        settings = read_settings()
        layout = QtWidgets.QFormLayout(self)
        self.key_panel = QtWidgets.QLineEdit(settings['api_key'])
        self.key_panel.setEchoMode(QtWidgets.QLineEdit.EchoMode.Password)
        self.key_panel.setPlaceholderText('API key, or leave blank to use OPENAI_API_KEY')
        self.comboBox_model = QtWidgets.QComboBox()
        self.comboBox_model.setEditable(True)
        self.comboBox_model.addItems(MODELS)
        self.comboBox_model.setCurrentText(settings['model'])
        self.comboBox_model.setToolTip('Luna: fast; Sol: balanced; Astra: complex work. '
                                       'You can enter another API model ID available to your account.')
        self.reasoning_box = QtWidgets.QComboBox()
        self.reasoning_box.addItems(['auto', 'low', 'medium', 'high'])
        self.reasoning_box.setCurrentText(settings['reasoning_effort'])
        self.reasoning_box.setToolTip('More reasoning can increase response time and token usage. '
                                     'Use auto for a model that does not accept this option.')
        self.verbosity_box = QtWidgets.QComboBox()
        self.verbosity_box.addItems(['auto', 'low', 'medium', 'high'])
        self.verbosity_box.setCurrentText(settings['verbosity'])
        self.verbosity_box.setToolTip('Default response detail. Use auto for models without this option.')
        self.max_tokens_box = QtWidgets.QSpinBox()
        self.max_tokens_box.setRange(128, 128000)
        self.max_tokens_box.setSingleStep(1024)
        self.max_tokens_box.setValue(settings['max_output_tokens'])
        self.max_tokens_box.setToolTip('Includes visible text and reasoning tokens. '
                                     'The model may have a lower limit.')
        layout.addRow('OpenAI API key', self.key_panel)
        layout.addRow('Model', self.comboBox_model)
        layout.addRow('Reasoning effort', self.reasoning_box)
        layout.addRow('Response detail', self.verbosity_box)
        layout.addRow('Output token budget', self.max_tokens_box)
        note = QtWidgets.QLabel('Uses your OpenAI API account and billing. '
                                'Only the messages in this chat are sent; fit data is not attached automatically.')
        note.setWordWrap(True)
        layout.addRow(note)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Save |
                                             QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.update_OpenAI_settings)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def update_OpenAI_settings(self):
        try:
            save_settings(dict(model=self.comboBox_model.currentText(),
                               api_key=self.key_panel.text(),
                               reasoning_effort=self.reasoning_box.currentText(),
                               verbosity=self.verbosity_box.currentText(),
                               max_output_tokens=self.max_tokens_box.value()))
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.warning(self, 'Cannot save settings', str(exc))
            return
        self.accept()

    @staticmethod
    def get_window(parent=None):
        return show_chat_opt(parent).exec()
