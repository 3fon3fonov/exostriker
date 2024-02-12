from PyQt6 import QtWidgets,QtGui 
import sys, os
import json
 
class show_chat_opt(QtWidgets.QDialog):


    def __init__(self, parent = None):
       # super(show_chat_opt, self).__init__(parent)
        super(show_chat_opt, self).__init__()

        self.layout = QtWidgets.QGridLayout(self)
        #self.layout = QtWidgets.QVBoxLayout(self)
       # self.parent = parent       
 
        self.title = 'This is a test window'

        self.widget=QtWidgets.QWidget(self)  # central widget

        self.comboBox_model = QtWidgets.QComboBox(self)
        self.key_panel = QtWidgets.QTextEdit('', self)
        self.key_panel.setFixedHeight(20)
        self.key_panel.setFixedWidth(500)

        self.text_label1 = QtWidgets.QLabel("Your OpenAI key:",self)
        self.text_label2 = QtWidgets.QLabel("model",self)
        self.text_label3 = QtWidgets.QLabel("temperature",self)
        self.text_label4 = QtWidgets.QLabel("max_tokens",self)
        self.text_label5 = QtWidgets.QLabel("top_p",self)
        self.text_label6 = QtWidgets.QLabel("frequency_penalty",self)
        self.text_label7 = QtWidgets.QLabel("presence_penalty",self)

        sett_init = self.read_OpenAI_settings()
        self.sett_init = sett_init

        self.temp_box  = QtWidgets.QDoubleSpinBox(self)
        self.temp_box.setMaximum(1.0)
        self.temp_box.setMinimum(0.0)
        self.temp_box.setValue(float(sett_init["openAI"]["temperature"])) 
        self.temp_box.setSingleStep(0.05) 

        self.max_tokens_box  = QtWidgets.QSpinBox(self)
        self.max_tokens_box.setMaximum(8096)
        self.max_tokens_box.setMinimum(10)
        self.max_tokens_box.setValue(int(sett_init["openAI"]["max_tokens"])) 
        self.max_tokens_box.setSingleStep(10) 

        self.top_p_box  = QtWidgets.QDoubleSpinBox(self)
        self.top_p_box.setMaximum(1.0)
        self.top_p_box.setMinimum(0.0)
        self.top_p_box.setValue(float(sett_init["openAI"]["top_p"])) 
        self.top_p_box.setSingleStep(0.05) 

        self.fr_pen_box  = QtWidgets.QDoubleSpinBox(self)
        self.fr_pen_box.setMaximum(2.0)
        self.fr_pen_box.setMinimum(-2.0)
        self.fr_pen_box.setValue(float(sett_init["openAI"]["frequency_penalty"])) 
        self.fr_pen_box.setSingleStep(0.1) 

        self.pr_pen_box = QtWidgets.QDoubleSpinBox(self)
        self.pr_pen_box.setMaximum(2.0)
        self.pr_pen_box.setMinimum(-2.0)
        self.pr_pen_box.setValue(float(sett_init["openAI"]["presence_penalty"])) 
        self.pr_pen_box.setSingleStep(0.1) 

        self.key_panel.setPlainText(str(sett_init["openAI"]["api_key"]))

        self.model = self.comboBox_model.currentText()

        self.comboBox_model.activated.connect(self.update_model)

        self.cancel_button = QtWidgets.QPushButton('Close', self)

        self.save_button = QtWidgets.QPushButton('Save', self)

        self.cancel_button.clicked.connect(self.close)
        self.save_button.clicked.connect(self.update_OpenAI_settings)


        self.init_buttons()     
        self.init_comboBox_model()

    def read_OpenAI_settings(self):
        if os.path.isfile('./lib/ES_settings_dev.json'):
            AIsettings = './lib/ES_settings_dev.json'
        else:
            AIsettings = './lib/ES_settings.json'

        with open(AIsettings, "r") as jsonFile:
            sett = json.load(jsonFile)

        return sett


    def update_OpenAI_settings(self):

        if os.path.isfile('./lib/ES_settings_dev.json'):
            AIsettings = './lib/ES_settings_dev.json'
        else:
            AIsettings = './lib/ES_settings.json'

        with open(AIsettings, "r") as jsonFile:
            sett = json.load(jsonFile)

        self.update_model()
        sett["openAI"]["model"] = str(self.model)
        sett["openAI"]["temperature"] = float(self.temp_box.value())
        sett["openAI"]["max_tokens"] = int(self.max_tokens_box.value())
        sett["openAI"]["frequency_penalty"] = float(self.fr_pen_box.value())
        sett["openAI"]["presence_penalty"] = float(self.pr_pen_box.value())
        sett["openAI"]["top_p"] = float(self.top_p_box.value())

        sett["openAI"]["api_key"] = str(self.key_panel.toPlainText())


        with open(AIsettings, "w") as jsonFile:
            json.dump(sett, jsonFile)



    def init_comboBox_model(self):

        self.comboBox_model.clear()

        model = ["text-davinci-001","text-davinci-002","text-davinci-003"]
        for i in range(len(model)):
            self.comboBox_model.addItem(model[i],i)    
        #self.comboBox_model.setCurrentIndex(0)
        self.comboBox_model.setCurrentIndex(int(int(str(self.sett_init["openAI"]["model"][-1]))-1))


    def update_model(self):
        self.model = self.comboBox_model.currentText()


    def init_buttons(self):
        
        self.layout.addWidget(self.text_label1, 0, 0)     
        self.layout.addWidget(self.text_label2, 1, 0)     
        self.layout.addWidget(self.text_label3, 2, 0)     
        self.layout.addWidget(self.text_label4, 3, 0)     
        self.layout.addWidget(self.text_label5, 4, 0)     
        self.layout.addWidget(self.text_label6, 5, 0)     
        self.layout.addWidget(self.text_label7, 6, 0)     

        self.layout.addWidget(self.key_panel, 0, 1)     
        self.layout.addWidget(self.comboBox_model, 1, 1)        
        self.layout.addWidget(self.temp_box, 2, 1)        
        self.layout.addWidget(self.max_tokens_box, 3, 1)        
        self.layout.addWidget(self.top_p_box, 4, 1)        
        self.layout.addWidget(self.fr_pen_box, 5, 1)        
        self.layout.addWidget(self.pr_pen_box, 6, 1)        

        self.layout.addWidget(self.cancel_button ,8, 3)
        self.layout.addWidget(self.save_button,7, 3)

    def return_but_N(self):
        #   Return list of values. It need map with str (self.lineedit.text() will return QString)
        return self.radio_group.checkedId() 

    # static method to create the dialog and return
    @staticmethod
    def get_window(parent = None):
        dialog = show_chat_opt(parent)
        result = dialog.exec_()
        return result



if __name__ == '__main__':
   
    app = QtWidgets.QApplication([])
    but = show_chat_opt.get_window()
    app.exec_()        
