import sys #,os
from PyQt6 import QtWidgets,QtGui 
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
        self.widget=QtWidgets.QWidget(self)  # central widget
      #  self.widget.setLayout(layout)    
        
        self.radio_group =QtWidgets.QButtonGroup(self.widget) # Number group
        self.radio_group2=QtWidgets.QButtonGroup(self.widget) # Number group
        self.radio_group3=QtWidgets.QButtonGroup(self.widget) # Number group
        self.radio_group4=QtWidgets.QButtonGroup(self.widget) # Number group
        
        
        self.radio_median     = QtWidgets.QRadioButton('median', self)
        self.radio_mean       = QtWidgets.QRadioButton('mean', self)
        self.radio_mode       = QtWidgets.QRadioButton('mode', self)
        self.radio_best_samp  = QtWidgets.QRadioButton('max lnL sample', self)
        self.radio_best_gui   = QtWidgets.QRadioButton('max lnL GUI', self) 
        self.radio_no_cross   = QtWidgets.QRadioButton('do not show', self) 
        
        #self.radio_mass_j     = QtWidgets.QRadioButton('Jupiter mass', self)
       # self.radio_mass_e       = QtWidgets.QRadioButton('Earth mass', self)        
        
        self.label_bestfit    = QtWidgets.QLabel("Truths as:")
        
        self.radio_group.addButton(self.radio_median,0)
        self.radio_group.addButton(self.radio_mean,1)
        self.radio_group.addButton(self.radio_mode,2)
        self.radio_group.addButton(self.radio_best_samp,3)
        self.radio_group.addButton(self.radio_best_gui,4)
        self.radio_group.addButton(self.radio_no_cross,5)

        self.radio_mp_So      = QtWidgets.QRadioButton('Solar', self)
        self.radio_mp_Mj      = QtWidgets.QRadioButton('Jupiter', self)
        self.radio_mp_Me      = QtWidgets.QRadioButton('Earth', self)        
 
    
        self.radio_group3.addButton(self.radio_mp_So,0)
        self.radio_group3.addButton(self.radio_mp_Mj,1)
        self.radio_group3.addButton(self.radio_mp_Me,2)
        
        
        self.radio_Rp_So      = QtWidgets.QRadioButton('Solar', self)
        self.radio_Rp_Mj      = QtWidgets.QRadioButton('Jupiter', self)
        self.radio_Rp_Me      = QtWidgets.QRadioButton('Earth', self)        
 
    
        self.radio_group4.addButton(self.radio_Rp_So,0)
        self.radio_group4.addButton(self.radio_Rp_Mj,1)
        self.radio_group4.addButton(self.radio_Rp_Me,2)        
        
    
        self.mass_check  = QtWidgets.QCheckBox('incl. mass (needed: K,P,e )', self)
        self.semi_check  = QtWidgets.QCheckBox('incl. semi-major axis (needed: P)', self)
        self.radi_check  = QtWidgets.QCheckBox('incl. radius (needed: Rp/Rs)', self)
        self.ppm_check   = QtWidgets.QCheckBox('rel. flux --> ppm', self)
        self.lambda_check   = QtWidgets.QCheckBox('incl. lambda (needed: e,omega )', self)


        
        #self.radi_check.setEnabled(False)
        
       # self.mass_check.setChecked(bool(self.parent.lables_cornerplot['mass']))

        self.radio_bestfit_labels = ["median","mean","mode","best_samp","best_gui",
                                     "none","mass","use_Me","use_Mj","use_Ms","semimajor","radius",
                                     "use_Re","use_Rj","use_Rs","use_lambda","use_ppm"]
        self.radio_group_list = [self.radio_median,self.radio_mean,self.radio_mode,
                                 self.radio_best_samp,self.radio_best_gui,self.radio_no_cross,
                                 self.mass_check,self.radio_mp_Me, self.radio_mp_Mj, self.radio_mp_So, 
                                 self.semi_check,self.radi_check,self.radio_Rp_Me, 
                                 self.radio_Rp_Mj, self.radio_Rp_So,self.lambda_check,self.ppm_check]
        
        #self.radio_group_mass_list = [self.radio_mp_So,self.radio_mp_Mj,self.radio_mp_Me,
       #                          self.radio_best_samp,self.radio_best_gui,self.radio_no_cross,
       #                          self.mass_check,self.semi_check,self.radi_check]        
#

        self.label_cornerplot_opt    = QtWidgets.QLabel("Cornerplot options:")


        self.cancel_button = QtWidgets.QPushButton('Close', self)
        
   
        self.truth_color_button = QtWidgets.QPushButton('Truths color', self)
        self.truth_color_button.clicked.connect(self.get_color_truth)


        self.samp_color_button = QtWidgets.QPushButton('Samples color', self)
        self.samp_color_button.clicked.connect(self.get_color_samp)

################################### Stab opt. ###########################
        self.stab_samp_color_button = QtWidgets.QPushButton('Stab. Samples color', self)
        self.stab_samp_color_button.clicked.connect(self.get_color_stab_samp)

        self.spin_max_stab_time  = QtWidgets.QDoubleSpinBox(self)
        self.spin_max_stab_time.setPrefix('stab. threshold = ')
        self.spin_max_stab_time.setSuffix(' [yr]')
        self.spin_max_stab_time.setRange(1, 10000000000)


        self.stab_samp_color_button.setHidden(True)
        self.spin_max_stab_time.setHidden(True)

#########################################################################

        
        self.check_plot_datapoints  = QtWidgets.QCheckBox('plot datapoints', self)
        self.check_plot_contours    = QtWidgets.QCheckBox('plot contours', self)
        self.check_plot_show_titles = QtWidgets.QCheckBox('show titles', self)
        self.check_plot_scale_hist  = QtWidgets.QCheckBox('scale hist', self)
        self.check_no_fill_contours = QtWidgets.QCheckBox('no fill contours', self)
        self.check_fill_contours    = QtWidgets.QCheckBox('fill contours', self)
        
        self.check_reverse          = QtWidgets.QCheckBox('reverse', self)
        
        
        self.spin_bins  = QtWidgets.QSpinBox(self)
        self.spin_bins.setSuffix(' bins')
        
        self.spin_label_pad  = QtWidgets.QDoubleSpinBox(self)
        self.spin_label_pad.setPrefix('label pad = ')        
        self.spin_label_pad.setMaximum(1.0)
        self.spin_label_pad.setMinimum(0.0)
        self.spin_label_pad.setSingleStep(0.05) 

        self.quantiles  = QtWidgets.QDoubleSpinBox(self)
        self.quantiles.setPrefix('quantiles = ')  
        self.quantiles.setMaximum(100.1)
        self.quantiles.setMinimum(0.0)
        self.quantiles.setSingleStep(0.1) 

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
            
            try:
                self.radio_mp_So.setChecked(self.parent.lables_cornerplot["use_Ms"])
                self.radio_mp_Mj.setChecked(self.parent.lables_cornerplot["use_Mj"])
                self.radio_mp_Me.setChecked(self.parent.lables_cornerplot["use_Me"])       
            except:
                self.parent.lables_cornerplot["use_Ms"] = False
                self.parent.lables_cornerplot["use_Mj"] = False
                self.parent.lables_cornerplot["use_Me"] = True
                self.radio_mp_So.setChecked(self.parent.lables_cornerplot["use_Ms"])
                self.radio_mp_Mj.setChecked(self.parent.lables_cornerplot["use_Mj"])
                self.radio_mp_Me.setChecked(self.parent.lables_cornerplot["use_Me"])      
            
            try:
                self.radio_Rp_So.setChecked(self.parent.lables_cornerplot["use_Rs"])
                self.radio_Rp_Mj.setChecked(self.parent.lables_cornerplot["use_Rj"])
                self.radio_Rp_Me.setChecked(self.parent.lables_cornerplot["use_Re"])       
            except:
                self.parent.lables_cornerplot["use_Rs"] = False
                self.parent.lables_cornerplot["use_Rj"] = False
                self.parent.lables_cornerplot["use_Re"] = True
                self.radio_Rp_So.setChecked(self.parent.lables_cornerplot["use_Rs"])
                self.radio_Rp_Mj.setChecked(self.parent.lables_cornerplot["use_Rj"])
                self.radio_Rp_Me.setChecked(self.parent.lables_cornerplot["use_Re"])                  
            try:    
                self.lambda_check.setChecked(self.parent.lables_cornerplot["use_lambda"])
            except:
                self.parent.lables_cornerplot["use_lambda"] = False
                self.lambda_check.setChecked(self.parent.lables_cornerplot["use_lambda"])  
                
            try:    
                self.ppm_check.setChecked(self.parent.lables_cornerplot["use_ppm"])
            except:
                self.parent.lables_cornerplot["use_ppm"] = False
                self.ppm_check.setChecked(self.parent.lables_cornerplot["use_ppm"])

            
            self.cornerplot_opt2 = {"bins":25,
                              "reverse":True,
                              "upper":True,
                              "quantiles":68.3,
                              "levels":(0.6827, 0.9545,0.9973), 
                              "smooth":1.0, 
                              "smooth1d":1.0, 
                              "dpi":300, 
                              "pad":15, 
                              "title_kwargs":{"fontsize": 12},  }      
            
            self.cornerplot_opt = self.parent.lables_cornerplot["cornerplot"]
            self.OrigLabels     = self.parent.lables_cornerplot["OrigLabels"]

#            try:
            if "stab_color" not in self.cornerplot_opt:
                self.cornerplot_opt["stab_color"] = 'r'
                
            if "stab_threshold" not in self.cornerplot_opt:
                self.cornerplot_opt["stab_threshold"] = 100000.0
 

            self.truth_color_button.setStyleSheet("color: %s;"%self.cornerplot_opt["truth_color"])
            self.samp_color_button.setStyleSheet("color: %s;"%self.cornerplot_opt["color"])

            self.stab_samp_color_button.setStyleSheet("color: %s;"%self.cornerplot_opt["stab_color"])
            
            self.check_plot_datapoints.setChecked(self.cornerplot_opt["plot_datapoints"])       
            self.check_plot_contours.setChecked(self.cornerplot_opt["plot_contours"])
            self.check_plot_show_titles.setChecked(self.cornerplot_opt["show_titles"])
            self.check_plot_scale_hist.setChecked(self.cornerplot_opt["scale_hist"])
            self.check_no_fill_contours.setChecked(self.cornerplot_opt["no_fill_contours"])
            
            try:
                self.check_fill_contours.setChecked(self.cornerplot_opt["fill_contours"])
            except:
                self.cornerplot_opt["fill_contours"] = False
                self.check_fill_contours.setChecked(self.cornerplot_opt["fill_contours"])
                
                
            self.check_reverse.setChecked(self.cornerplot_opt["reverse"])   
            
            self.spin_bins.setValue(self.cornerplot_opt["bins"])
            self.spin_label_pad.setValue(self.cornerplot_opt["labelpad"])

            if self.cornerplot_opt["quantiles"] == None:
                self.quantiles.setValue(100.1)            
            else:
                self.quantiles.setValue(self.cornerplot_opt["quantiles"])

            
            #print(self.OrigLabels)            
        else:
            self.truth_color_button.setStyleSheet("color: #ff0000;")
            self.samp_color_button.setStyleSheet("color: #ff0000;")

        


        self.initialize_color_dialog()
        self.init_buttons()  
        

    
    def init_buttons(self):


        
        k = 0
        l = 0
        for g in range(len(self.parent.lables_cornerplot)-19):

            k = g%20

            if g != 0 and k == 0:
                l = l+2
 
                
            self.button  = QtWidgets.QCheckBox('', self)
            self.button.setChecked(self.parent.lables_cornerplot[g][1]=="True")


            self.text_panel = QtWidgets.QTextEdit('%s'%self.parent.lables_cornerplot[g][0], self)
 

           
            self.button.setFixedHeight(20)
            self.text_panel.setFixedHeight(20)
            self.text_panel.setFixedWidth(200)
           
            self.layout.addWidget(self.button, k, l)
            self.layout.addWidget(self.text_panel,k, l+1)

            self.checked_boxes.append(self.button)
            self.text_boxes.append(self.text_panel) 
#            k = g
            #for 
 
        #for g in range(4):
        self.layout.addWidget(self.mass_check, 0, l+2)
        self.layout.addWidget(self.radio_mp_Me, 1,  l+2)
        self.layout.addWidget(self.radio_mp_Mj, 2,  l+2)         
        self.layout.addWidget(self.radio_mp_So, 3, l+2)        
        self.layout.addWidget(self.semi_check, 4, l+2)
        self.layout.addWidget(self.radi_check, 5, l+2)
        self.layout.addWidget(self.radio_Rp_Me, 6,  l+2)
        self.layout.addWidget(self.radio_Rp_Mj, 7,  l+2)         
        self.layout.addWidget(self.radio_Rp_So, 8, l+2) 
        self.layout.addWidget(self.lambda_check, 9, l+2)         
        self.layout.addWidget(self.ppm_check, 10, l+2) 
         
        self.layout.addWidget(self.label_bestfit, 0,l+3)              
        self.layout.addWidget(self.radio_median, 1,l+3)
        self.layout.addWidget(self.radio_mean,2,l+3)
        self.layout.addWidget(self.radio_mode,3,l+3)
        self.layout.addWidget(self.radio_best_samp,4,l+3)
        self.layout.addWidget(self.radio_best_gui, 5,l+3)  
        self.layout.addWidget(self.radio_no_cross, 6,l+3)  
        self.layout.addWidget(self.truth_color_button, 11,l+3)  
        
        
 
        self.layout.addWidget(self.label_cornerplot_opt, 0,l+4)
        self.layout.addWidget(self.spin_bins, 1,l+4)
        
        
        self.layout.addWidget(self.check_plot_datapoints, 2,l+4) 
        self.layout.addWidget(self.check_plot_contours, 3,l+4) 
        self.layout.addWidget(self.check_plot_show_titles, 4,l+4) 
        self.layout.addWidget(self.check_plot_scale_hist, 5,l+4) 
        self.layout.addWidget(self.check_no_fill_contours, 6,l+4) 
        self.layout.addWidget(self.check_fill_contours, 7,l+4) 
        
        self.layout.addWidget(self.check_reverse, 8,l+4) 
        self.layout.addWidget(self.spin_label_pad, 9,l+4)
        self.layout.addWidget(self.quantiles, 10,l+4)          


        self.layout.addWidget(self.samp_color_button, 11,l+4)       

        try: 
            if "max. time" in self.OrigLabels:

                self.stab_samp_color_button.setHidden(False)
                self.spin_max_stab_time.setHidden(False)

                self.layout.addWidget(self.stab_samp_color_button, 11,l+4)     
                self.layout.addWidget(self.spin_max_stab_time, 12,l+4)     
                self.layout.addWidget(self.cancel_button, 15,l+4)
                self.spin_max_stab_time.setValue(self.cornerplot_opt["stab_threshold"])

            else:
                self.layout.addWidget(self.cancel_button, 13,l+4)           
        except:
            self.layout.addWidget(self.cancel_button, 13,l+4)
        self.cancel_button.clicked.connect(self.close)       
            
        
        #self.setCentralWidget(self.widget)
        #self.Ok_button = QtGui.QPushButton('OK', self)
        #self.layout.addWidget(self.Ok_button)        
        


        #self.Ok_button.clicked.connect(self.get_radio)
        #self.layout.setRowStretch(k, 1)
        #self.layout.setColumnStretch(k, 1)
        
        


    def return_but_N(self):

        #for g in range(len(self.parent.lables_cornerplot)):
            
        results = {k: np.array([self.text_boxes[k].toPlainText(), self.checked_boxes[k].isChecked()]) for k in range(len(self.checked_boxes))}
        
        for k in range(17):
            #print(self.radio_bestfit_labels[k],self.radio_group_list[k].isChecked() )
            results[self.radio_bestfit_labels[k]] = self.radio_group_list[k].isChecked() 
        
        #print(np.array([[bestfit_labels[k],bestfit_labels_bool[k]] for k in range(6)]) )  
        
        self.cornerplot_opt["plot_datapoints"] = self.check_plot_datapoints.isChecked()
        self.cornerplot_opt["plot_contours"]   = self.check_plot_contours.isChecked()
        self.cornerplot_opt["show_titles"]     = self.check_plot_show_titles.isChecked()
        self.cornerplot_opt["scale_hist"]      = self.check_plot_scale_hist.isChecked()
        self.cornerplot_opt["no_fill_contours"]= self.check_no_fill_contours.isChecked()
        self.cornerplot_opt["fill_contours"]   = self.check_fill_contours.isChecked()

        self.cornerplot_opt["reverse"]         = self.check_reverse.isChecked()
          
        self.cornerplot_opt["bins"]            = self.spin_bins.value()
        self.cornerplot_opt["labelpad"]        = self.spin_label_pad.value()
        self.cornerplot_opt["stab_threshold"]   = self.spin_max_stab_time.value()

        if float(self.quantiles.value()) > 100.0:
            self.cornerplot_opt["quantiles"]        = None   
        else:     
            self.cornerplot_opt["quantiles"]        = self.quantiles.value()     

        
        results["cornerplot"] = self.cornerplot_opt
        results["OrigLabels"] = self.OrigLabels
        
        return results
 
    def get_color_truth(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog|QtWidgets.QColorDialog.ShowAlphaChannel,)
 
        if colorz.isValid(): 
            self.cornerplot_opt["truth_color"]=colorz.name() 
            self.truth_color_button.setStyleSheet("color: %s;"%colorz.name())

        else:
            return
        
    def get_color_samp(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog|QtWidgets.QColorDialog.ShowAlphaChannel,)
 
        if colorz.isValid(): 
            self.cornerplot_opt["color"]=colorz.name() 
            self.samp_color_button.setStyleSheet("color: %s;"%colorz.name())

        else:
            return        


    def get_color_stab_samp(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog|QtWidgets.QColorDialog.ShowAlphaChannel,)
 
        if colorz.isValid(): 
            self.cornerplot_opt["stab_color"]=colorz.name() 
            self.stab_samp_color_button.setStyleSheet("color: %s;"%colorz.name())

        else:
            return        

        

    def initialize_color_dialog(self):

        self.colorDialog = QtWidgets.QColorDialog()
        self.colorDialog.setOption(QtWidgets.QColorDialog.ShowAlphaChannel, True)
        self.colorDialog.setOption(QtWidgets.QColorDialog.DontUseNativeDialog, True)

 
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
    app = QtWidgets.QApplication([])
    but = show_param_boxes.DateDialog.get_radio()
#    print("{} {} {}".format(date, time, ok))
    app.exec_()        
