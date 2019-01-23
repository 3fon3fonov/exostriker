 #!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'

class Error(Exception):
    pass

class InputError(Error):
 
    def __init__(self, message):
        self.message = message

class FittingError(Error):
 
    def __init__(self, message):
        self.message = message
      