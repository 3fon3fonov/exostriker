import openai
import sys, os
#import numpy as np
import json
 

from PyQt6.QtCore import Qt, QTextStream, QIODevice
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, 
                             QTextEdit, QLineEdit, QPushButton)
from PyQt6.QtGui import QColor,QIcon

from exochat_opt import show_chat_opt

#try:
#    with open('./openai_api_key.txt', 'r') as fileai: 
#         api_key= str(fileai.readline().rstrip())
#    openai.api_key = api_key
#    key_found = True
#except:
#    openai.api_key = ''
#    key_found = False         




class ChatWidget(QWidget):
    def __init__(self):
        super().__init__()
 
        self.setGeometry(1,1, 495, 325) 
       
        self.conversation = QTextEdit()
        self.conversation.setReadOnly(True)
        self.message = QLineEdit()
        self.send_button = QPushButton("Send")
        self.send_button.clicked.connect(self.send_message)
        
        self.opt_button = QPushButton("")
        self.opt_button.setIcon(QIcon('./lib/UI/opt_icon.png'))
        self.opt_button.clicked.connect(self.get_symbol)

        self.dialog_symbols = show_chat_opt(self)

        # Create the layout
        layout = QVBoxLayout()
        input_layout = QHBoxLayout()
        input_layout.addWidget(self.message)
        input_layout.addWidget(self.opt_button)
        input_layout.addWidget(self.send_button)
        layout.addWidget(self.conversation)
        layout.addLayout(input_layout)
        self.setLayout(layout)


        #self.redColor = QColor(255, 0, 0)
        self.blackColor = QColor(0, 0, 0)
        self.blueColor = QColor(0, 0, 255)    

        self.read_settings()

 

    def read_settings(self):
 
        if os.path.isfile('./lib/ES_settings_dev.json'):
            AIsettings = './lib/ES_settings_dev.json'
        else:
            AIsettings = './lib/ES_settings.json'


        f = open(AIsettings)
        sett = json.load(f)
        f.close()

        self.model = str(sett["openAI"]["model"])
        self.temperature = float(sett["openAI"]["temperature"])
        self.top_p = float(sett["openAI"]["top_p"])
        self.max_tokens = int(sett["openAI"]["max_tokens"])
        self.frequency_penalty = float(sett["openAI"]["frequency_penalty"])
        self.presence_penalty = float(sett["openAI"]["presence_penalty"])
        self.api_key=str(sett["openAI"]["api_key"])

        openai.api_key = self.api_key



        return 


    def get_symbol(self):

        but_n = self.dialog_symbols.get_window()
            
       # if but_n != None:
       #     print("Test") 
       #     model=str(self.dialog_symbols.model)
      #      print(model)
       # else:
      #      return  

    def send_message(self):
        # Get the message from the QLineEdit
        message = self.message.text()
  
        # Clear the QLineEdit
        self.message.clear()

        #if key_found == False:
       #     self.conversation.append(f"You do not have a valid openai api key\n") 
        
        # Append the message to the conversation
        self.conversation.setTextColor(self.blackColor)
        self.conversation.append(f"You: {message}\n")

        
        # Use ChatGPT to generate a response

        try:
            self.read_settings()
            response = openai.Completion.create(
              model=self.model,
              prompt=f"You: {message}\n",
              temperature=self.temperature,
              max_tokens=self.max_tokens,
              top_p=self.top_p,
              frequency_penalty=self.frequency_penalty,
              presence_penalty=self.presence_penalty,
              stop=["You:"]
            )         
     
            #print(response)
            # Append the response to the conversation
            self.conversation.setTextColor(self.blueColor)
            self.conversation.append(f"{response['choices'][0]['text']}\n")
            self.conversation.setTextColor(self.blackColor)

#        except openai.error.RateLimitError: 
#            self.conversation.append(f"openai.error.RateLimitError: You exceeded your current quota, please check your plan and billing details.")
#        except openai.error.AuthenticationError: 
#            self.conversation.append(f"openai.error.AuthenticationError: Incorrect API key provided.")
 
        except Exception as e:
            self.conversation.append(f"Your Exo-Striker ChatBoot does not work! Mostlikely you: \n \n  * Have not provided, or you do not have a valid openai API key. You can obtain an API key from https://beta.openai.com. Then add your API key as a single pressing the 'option' button, and try again! \n \n * You do not have an internet! \n The Exo-Striker ChatBoot is a clone of Chat GPT, which requres an internet connection. \n \n * Or... maybe OpenAI is busy... \n" )    

            self.conversation.append(f"The exact error was: \n %s\n"%e)            
        
        # Scroll to the bottom of the conversation
        self.conversation.verticalScrollBar().setValue(self.conversation.verticalScrollBar().maximum())

if __name__ == "__main__":
    app = QApplication(sys.argv)
    chat = ChatWidget()
    chat.show()
    sys.exit(app.exec_())

