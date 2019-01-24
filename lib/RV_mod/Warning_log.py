 #!/usr/bin/python
__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
#sys.path.insert(0, '../lib')
 

  
class Warning_log(object): # here we can store all warnings raised when calling a certain function, so we can print them out nicely 
 
    warning_count=0
    warning_messages=[]
    procedure_name='Your query' # what function's warnings are we storing, will be needed for printing warning log
 
    def __init__(self,warning_messages, procedure_name):
        self.warning_count=len(warning_messages)
        self.warning_messages=warning_messages
        self.procedure_name=procedure_name
 
    def update_warning_list(self,warning_message): # add a new warning message to the warning list 
        self.warning_messages.append(warning_message)
        self.warning_count = self.warning_count+1
        return
 
    def print_warning_log(self):
        if(self.warning_count==0): # do not print a warning log if there were no warnings
            pass
        else: # printing a custom warning log with a list of all warnings encountered when running the procedure
            print('')
            print('%s resulted in %d warnings:'%(self.procedure_name,self.warning_count))
            print('')
            for i in range(self.warning_count):
                print('%d %s'%(i+1,self.warning_messages[i]))
            print('') 
        return
