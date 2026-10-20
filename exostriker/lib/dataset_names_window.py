"""Edit display names of imported datasets without changing their paths."""
from PyQt6 import QtWidgets


# title, collection, display-name field, row prefix, GUI button group
DATASET_GROUPS = (
    ('RV', 'rv_data_sets', 4, '', 'buttonGroup_add_RV_data'),
    ('Transits', 'tra_data_sets', -1, '', 'buttonGroup_transit_data'),
    ('TTV', 'ttv_data_sets', 5, '', 'buttonGroup_ttv_data'),
    ('Astrometry', 'ast_data_sets', 7, 'Imaging', 'buttonGroup_ast_data'),
    ('Astrometry', 'ast_data_sets_hipp', 4, 'Hipparcos', 'buttonGroup_ast_data_2'),
    ('Astrometry', 'ast_data_sets_gaia', -1, 'Gaia', 'buttonGroup_ast_data_3'),
)


def imported_dataset_names(session):
    """Yield occupied slots, preserving actual indices and each record layout."""
    for section, attribute, name_index, prefix, button_group in DATASET_GROUPS:
        collection = getattr(session, attribute, {})
        items = collection.items() if hasattr(collection, 'items') else enumerate(collection)
        for slot, record in sorted(items):
            if len(record):
                yield section, attribute, slot, name_index, prefix, button_group, record


class DatasetNamesDialog(QtWidgets.QDialog):
    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self.setWindowTitle('Edit file names')
        self.resize(700, 620)
        layout = QtWidgets.QVBoxLayout(self)
        description = QtWidgets.QLabel('Change the dataset names shown in Exo-Striker.', self)
        description.setWordWrap(True)
        layout.addWidget(description)
        scroll = QtWidgets.QScrollArea(self)
        scroll.setWidgetResizable(True)
        contents = QtWidgets.QWidget(scroll)
        sections_layout = QtWidgets.QVBoxLayout(contents)
        self.sections = {}
        for title in ('RV', 'Transits', 'TTV', 'Astrometry'):
            group = QtWidgets.QGroupBox(title, contents)
            self.sections[title] = QtWidgets.QFormLayout(group)
            sections_layout.addWidget(group)
        self.editors = {}
        self.records = {}
        self.changed_datasets = set()
        for section, attribute, slot, name_index, prefix, _, record in imported_dataset_names(session):
            key = (attribute, slot)
            field = QtWidgets.QLineEdit(str(record[name_index]), contents)
            label = f'{prefix} #{slot + 1}' if prefix else f'#{slot + 1}'
            field.setAccessibleName(f'{section} {label} dataset name')
            field.setMinimumWidth(320)
            self.sections[section].addRow(label, field)
            self.editors[key] = field
            self.records[key] = (record, name_index)
        for form in self.sections.values():
            if form.rowCount() == 0:
                form.addRow(QtWidgets.QLabel('No imported datasets', contents))
        sections_layout.addStretch(1)
        scroll.setWidget(contents)
        layout.addWidget(scroll, 1)
        self.message = QtWidgets.QLabel(self)
        self.message.setWordWrap(True)
        layout.addWidget(self.message)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Apply |
            QtWidgets.QDialogButtonBox.StandardButton.Cancel, parent=self)
        self.apply_button = buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Apply)
        self.apply_button.setEnabled(bool(self.editors))
        self.apply_button.clicked.connect(self.apply_names)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def apply_names(self):
        # Validate every row before changing anything. Cancel/window-close never
        # touch the session, and a dataset replaced during editing is protected.
        changes = []
        for key, field in self.editors.items():
            attribute, slot = key
            original, name_index = self.records[key]
            current = getattr(self.session, attribute)[slot]
            if current is not original:
                self.message.setText('The imported datasets changed. Close and reopen this window.')
                return
            name = field.text().strip()
            if not name:
                self.message.setText('Enter a non-empty name for each imported dataset.')
                field.setFocus()
                return
            if name != str(original[name_index]):
                changes.append((attribute, slot, original, name_index, name))
        for attribute, slot, record, name_index, name in changes:
            record[name_index] = name
            initial = getattr(self.session, attribute + '_init', {})
            try:
                initial_record = initial[slot]
            except (KeyError, IndexError):
                initial_record = []
            if len(initial_record):
                initial_record[name_index] = name
            self.changed_datasets.add(attribute)
        self.accept()
