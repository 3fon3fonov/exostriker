import sys #,os
from PyQt5 import QtWidgets,QtGui 
import numpy as np
 
class show_param_boxes(QtWidgets.QDialog):


    #def __init__(self, parent = None):
       # super(show_param_boxes, self).__init__(parent)
    #    super(show_param_boxes, self).__init__()
    def __init__(self,parent):
        super(show_param_boxes, self).__init__()
        
        #self.layout1 = QtWidgets.QVBoxLayout(self)   
        self.layout = QtWidgets.QGridLayout(self)
        self.title = 'This is a test window'
        
        self.parent = parent

        
        

       
        self.checked_boxes = []
        self.text_boxes = []       
       # self.setFixedSize(550, 800)        
       # self.widget = QtWidgets.QWidget(self)
       # self.widget.setLayout(QtWidgets.QGridLayout())
       # self.widget.layout().addWidget(self.text)

       #self.layout=QtWidgets.QGridLayout() # layout for the central widget
        self.widget=QtGui.QWidget(self)  # central widget
      #  self.widget.setLayout(layout)    
        
        self.radio_group =QtGui.QButtonGroup(self.widget) # Number group
        self.radio_group2=QtGui.QButtonGroup(self.widget) # Number group
        
        
        self.radio_median     = QtWidgets.QRadioButton('median', self)
        self.radio_mean       = QtWidgets.QRadioButton('mean', self)
        self.radio_mode       = QtWidgets.QRadioButton('mode', self)
        self.radio_best_samp  = QtWidgets.QRadioButton('max lnL sample', self)
        self.radio_best_gui   = QtWidgets.QRadioButton('max lnL GUI', self) 
        self.radio_no_cross   = QtWidgets.QRadioButton('do not show', self) 
        
        self.label_bestfit    = QtWidgets.QLabel("Truths as:")
        
        self.radio_group.addButton(self.radio_median,0)
        self.radio_group.addButton(self.radio_mean,1)
        self.radio_group.addButton(self.radio_mode,2)
        self.radio_group.addButton(self.radio_best_samp,3)
        self.radio_group.addButton(self.radio_best_gui,4)
        self.radio_group.addButton(self.radio_no_cross,5)

     #   self.radio_mp_So      = QtWidgets.QRadioButton('Solar', self)
     #   self.radio_mp_Mj      = QtWidgets.QRadioButton('Jupiter', self)
     #   self.radio_mp_Me      = QtWidgets.QRadioButton('Earth', self)        
 
        self.mass_check  = QtWidgets.QCheckBox('incl. mass (needed: K,P,e )', self)
        self.semi_check  = QtWidgets.QCheckBox('incl. semi-major axis (needed: P)', self)
        self.radi_check  = QtWidgets.QCheckBox('incl. radius (needed: Rp/Rs)', self)
        
        #self.radi_check.setEnabled(False)
        
       # self.mass_check.setChecked(bool(self.parent.lables_cornerplot['mass']))

        self.radio_bestfit_labels = ["median","mean","mode","best_samp","best_gui","none","mass","semimajor","radius"]
        self.radio_group_list = [self.radio_median,self.radio_mean,self.radio_mode,
                                 self.radio_best_samp,self.radio_best_gui,self.radio_no_cross,
                                 self.mass_check,self.semi_check,self.radi_check]
        
        #self.radio_group_mass_list = [self.radio_mp_So,self.radio_mp_Mj,self.radio_mp_Me,
       #                          self.radio_best_samp,self.radio_best_gui,self.radio_no_cross,
       #                          self.mass_check,self.semi_check,self.radi_check]        
#

        self.label_cornerplot_opt    = QtWidgets.QLabel("Cornerplot options:")


        self.cancel_button = QtGui.QPushButton('Close', self)
        
   

        self.truth_color_button = QtGui.QPushButton('Truths color', self)
        self.truth_color_button.clicked.connect(self.get_color_truth)


        self.samp_color_button = QtGui.QPushButton('Samples color', self)
        self.samp_color_button.clicked.connect(self.get_color_samp)
        
        self.check_plot_datapoints  = QtWidgets.QCheckBox('plot datapoints', self)
        self.check_plot_contours    = QtWidgets.QCheckBox('plot contours', self)
        self.check_plot_show_titles = QtWidgets.QCheckBox('show titles', self)
        self.check_plot_scale_hist  = QtWidgets.QCheckBox('scale hist', self)
        self.check_no_fill_contours = QtWidgets.QCheckBox('no fill contours', self)
        
        
        self.spin_bins  = QtWidgets.QSpinBox(self)
        self.spin_bins.setSuffix(' bins')
        

        if len(self.parent.lables_cornerplot)!=0:
            self.radio_median.setChecked(bool(self.parent.lables_cornerplot["median"]))
            self.radio_mean.setChecked(self.parent.lables_cornerplot["mean"])
            self.radio_mode.setChecked(self.parent.lables_cornerplot["mode"])
            self.radio_best_samp.setChecked(self.parent.lables_cornerplot["best_samp"])
            self.radio_best_gui.setChecked(self.parent.lables_cornerplot["best_gui"])
            self.radio_no_cross.setChecked(self.parent.lables_cornerplot["none"])
            self.mass_check.setChecked(self.parent.lables_cornerplot["mass"])
            self.semi_check.setChecked(self.parent.lables_cornerplot["semimajor"])
            self.radi_check.setChecked(self.parent.lables_cornerplot["radius"])
            
            self.cornerplot_opt2 = {"bins":25,
                              "reverse":True,
                              "upper":True,
                              "quantiles":68.3,
                              "levels":(0.6827, 0.9545,0.9973), 
                              "smooth":1.0, 
                              "smooth1d":1.0, 
                              "dpi":300, 
                              "pad":15, 
                              "labelpad":50,
                              "title_kwargs":{"fontsize": 12},  }      
            
            self.cornerplot_opt = self.parent.lables_cornerplot["cornerplot"]
            self.OrigLabels     = self.parent.lables_cornerplot["OrigLabels"]
            self.truth_color_button.setStyleSheet("color: %s;"%self.cornerplot_opt["truth_color"])
            self.samp_color_button.setStyleSheet("color: %s;"%self.cornerplot_opt["color"])
            
            self.check_plot_datapoints.setChecked(self.cornerplot_opt["plot_datapoints"])       
            self.check_plot_contours.setChecked(self.cornerplot_opt["plot_contours"])
            self.check_plot_show_titles.setChecked(self.cornerplot_opt["show_titles"])
            self.check_plot_scale_hist.setChecked(self.cornerplot_opt["scale_hist"])
            self.check_no_fill_contours.setChecked(self.cornerplot_opt["no_fill_contours"])
           
            self.spin_bins.setValue(self.cornerplot_opt["bins"])
#            print(self.OrigLabels)            
        else:
            self.truth_color_button.setStyleSheet("color: #ff0000;")
            self.samp_color_button.setStyleSheet("color: #ff0000;")

        


        self.initialize_color_dialog()
        self.init_buttons()  
        

    
    def init_buttons(self):


        
        k = 0
        for g in range(len(self.parent.lables_cornerplot)-11):

            
            
            self.button  = QtWidgets.QCheckBox('', self)
            self.button.setChecked(self.parent.lables_cornerplot[g][1]=="True")


            self.text_panel = QtWidgets.QTextEdit('%s'%self.parent.lables_cornerplot[g][0], self)
 

           
            self.button.setFixedHeight(20)
            self.text_panel.setFixedHeight(20)
            self.text_panel.setFixedWidth(200)
           
            self.layout.addWidget(self.button, g, 0)
            self.layout.addWidget(self.text_panel,g, 1)

            self.checked_boxes.append(self.button)
            self.text_boxes.append(self.text_panel) 
            k = g
            #for 
            
            
        #for g in range(4):
        self.layout.addWidget(self.mass_check, k+1, 0, 1, 2)
        self.layout.addWidget(self.semi_check, k+2, 0, 1, 2)
        self.layout.addWidget(self.radi_check, k+3, 0, 1, 2)

         
        self.layout.addWidget(self.label_bestfit, 0,2)              
        self.layout.addWidget(self.radio_median, 1,2)
        self.layout.addWidget(self.radio_mean,2,2)
        self.layout.addWidget(self.radio_mode,3,2)
        self.layout.addWidget(self.radio_best_samp,4,2)
        self.layout.addWidget(self.radio_best_gui, 5,2)  
        self.layout.addWidget(self.radio_no_cross, 6,2)  
        self.layout.addWidget(self.truth_color_button, 7,2)  
        
        
 
        self.layout.addWidget(self.label_cornerplot_opt, 0,3)
        self.layout.addWidget(self.spin_bins, 1,3)
        
        
        self.layout.addWidget(self.check_plot_datapoints, 2,3) 
        self.layout.addWidget(self.check_plot_contours, 3,3) 
        self.layout.addWidget(self.check_plot_show_titles, 4,3) 
        self.layout.addWidget(self.check_plot_scale_hist, 5,3) 
        self.layout.addWidget(self.check_no_fill_contours, 6,3) 

        
             
        self.layout.addWidget(self.samp_color_button, 7,3)              
        
        
        
        #self.setCentralWidget(self.widget)
        #self.Ok_button = QtGui.QPushButton('OK', self)
        #self.layout.addWidget(self.Ok_button)        
        

        self.layout.addWidget(self.cancel_button, k+1+6,2, 2, 3)

        self.cancel_button.clicked.connect(self.close)
        #self.Ok_button.clicked.connect(self.get_radio)
        #self.layout.setRowStretch(k, 1)
        #self.layout.setColumnStretch(k, 1)
        
        


    def return_but_N(self):

        #for g in range(len(self.parent.lables_cornerplot)):
            
        results = {k: np.array([self.text_boxes[k].toPlainText(), self.checked_boxes[k].isChecked()]) for k in range(len(self.checked_boxes))}
        
        for k in range(9):
            #print(self.radio_bestfit_labels[k],self.radio_group_list[k].isChecked() )
            results[self.radio_bestfit_labels[k]] = self.radio_group_list[k].isChecked() 
        
        #print(np.array([[bestfit_labels[k],bestfit_labels_bool[k]] for k in range(6)]) )  
        
        self.cornerplot_opt["plot_datapoints"] = self.check_plot_datapoints.isChecked()
        self.cornerplot_opt["plot_contours"]   = self.check_plot_contours.isChecked()
        self.cornerplot_opt["show_titles"]     = self.check_plot_show_titles.isChecked()
        self.cornerplot_opt["scale_hist"]      = self.check_plot_scale_hist.isChecked()
        self.cornerplot_opt["no_fill_contours"]= self.check_no_fill_contours.isChecked()
          
        self.cornerplot_opt["bins"]            = self.spin_bins.value()
           
      
        
        results["cornerplot"] = self.cornerplot_opt
        results["OrigLabels"] = self.OrigLabels
        
        return results
 
    def get_color_truth(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtGui.QColorDialog.DontUseNativeDialog|QtGui.QColorDialog.ShowAlphaChannel,)
 
        if colorz.isValid(): 
            self.cornerplot_opt["truth_color"]=colorz.name() 
            self.truth_color_button.setStyleSheet("color: %s;"%colorz.name())

        else:
            return
        
    def get_color_samp(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtGui.QColorDialog.DontUseNativeDialog|QtGui.QColorDialog.ShowAlphaChannel,)
 
        if colorz.isValid(): 
            self.cornerplot_opt["color"]=colorz.name() 
            self.samp_color_button.setStyleSheet("color: %s;"%colorz.name())

        else:
            return        
        

    def initialize_color_dialog(self):

        self.colorDialog = QtGui.QColorDialog()
        self.colorDialog.setOption(QtGui.QColorDialog.ShowAlphaChannel, True)
        self.colorDialog.setOption(QtGui.QColorDialog.DontUseNativeDialog, True)

 
    # static method to create the dialog and return button pressed
    @staticmethod
    def get_labels(self):
        dialog = show_param_boxes(self)
        result = dialog.exec_()
        rad_but = dialog.return_but_N()
        return (rad_but)



if __name__ == '__main__':
   # app = QtWidgets.QApplication(sys.argv)
    #w = show_symbols()
   # w.show()
    #sys.exit(app.exec_())    
    app = QtGui.QApplication([])
    but = show_param_boxes.DateDialog.get_radio()
#    print("{} {} {}".format(date, time, ok))
    app.exec_()        