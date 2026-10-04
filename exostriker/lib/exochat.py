"""OpenAI chat widget with worker-thread requests and local conversation state."""
import sys
import os
from PyQt6.QtCore import QObject, QRunnable, QThreadPool, pyqtSignal, pyqtSlot
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QLineEdit, QPushButton, QLabel
from PyQt6.QtGui import QTextCursor
from exochat_opt import show_chat_opt
from exochat_backend import read_settings, request_reply


class ReplySignals(QObject):
    ready = pyqtSignal(object)
    failed = pyqtSignal(str)


class ReplyWorker(QRunnable):
    def __init__(self, settings, history, message):
        super().__init__()
        self.settings, self.history, self.message = settings, list(history), message
        self.signals = ReplySignals()

    @pyqtSlot()
    def run(self):
        try:
            result = request_reply(self.settings, self.history, self.message)
        except Exception as exc:
            # SDK errors can contain request details; never echo the API key.
            detail = str(exc)
            key = self.settings.get('api_key', '').strip() or os.environ.get('OPENAI_API_KEY', '').strip()
            if key:
                detail = detail.replace(key, '[redacted]')
            self.signals.failed.emit(detail)
        else:
            self.signals.ready.emit(result)


class ChatWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.history = []
        self.worker = None
        self.conversation = QTextEdit()
        self.conversation.setReadOnly(True)
        self.message = QLineEdit()
        self.message.setPlaceholderText('Ask about Exo-Striker or exoplanet analysis...')
        self.send_button = QPushButton('Send')
        self.send_button.clicked.connect(self.send_message)
        self.message.returnPressed.connect(self.send_message)
        self.opt_button = QPushButton('Settings')
        self.opt_button.clicked.connect(self.get_symbol)
        self.clear_button = QPushButton('New chat')
        self.clear_button.clicked.connect(self.clear_chat)
        self.status = QLabel('')
        layout = QVBoxLayout(self)
        layout.addWidget(self.conversation)
        layout.addWidget(self.status)
        inputs = QHBoxLayout()
        for widget in (self.message, self.opt_button, self.clear_button, self.send_button):
            inputs.addWidget(widget)
        layout.addLayout(inputs)

    def get_symbol(self):
        show_chat_opt.get_window(self)

    def clear_chat(self):
        self.history = []
        self.conversation.clear()
        self.status.clear()

    def append_text(self, speaker, text):
        # Plain-text insertion prevents pasted user/model HTML from rendering.
        self.conversation.moveCursor(QTextCursor.MoveOperation.End)
        self.conversation.insertPlainText(f'{speaker}: {text}\n\n')
        self.conversation.ensureCursorVisible()

    def set_busy(self, busy):
        for widget in (self.send_button, self.opt_button, self.clear_button):
            widget.setEnabled(not busy)
        self.status.setText('Thinking...' if busy else '')

    def send_message(self):
        message = self.message.text().strip()
        if not message or self.worker is not None:
            return
        try:
            settings = read_settings()
        except (OSError, ValueError) as exc:
            self.append_text('Error', str(exc))
            return
        self.message.clear()
        self.append_text('You', message)
        self.set_busy(True)
        self.worker = ReplyWorker(settings, self.history, message)
        self.worker.signals.ready.connect(self.receive_reply)
        self.worker.signals.failed.connect(self.receive_error)
        QThreadPool.globalInstance().start(self.worker)

    def receive_reply(self, result):
        text, self.history, notice = result
        self.worker = None
        self.set_busy(False)
        self.append_text('Assistant', text)
        if notice:
            self.status.setText(notice)

    def receive_error(self, error):
        self.worker = None
        self.set_busy(False)
        self.append_text('Error', error)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    chat = ChatWidget()
    chat.show()
    sys.exit(app.exec())
