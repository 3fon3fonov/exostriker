import sys
from PyQt5 import QtWidgets,QtGui,QtCore
import io
import numpy as np
from print_info_window import print_info


if sys.version_info[0] == 3:
    from urllib.request import urlopen
elif sys.version_info[0] == 2:
    from urllib2 import urlopen

#font = QtGui.QFont()
#font.setPointSize(8)
#font.setBold(False)


class RVBank_window(QtWidgets.QDialog):


    def __init__(self, parent=None):
        super(RVBank_window, self).__init__(parent=parent)
        vLayout = QtWidgets.QVBoxLayout(self)
        hLayout = QtWidgets.QHBoxLayout()


        self.info_dialog = print_info(self)


        self.lineEdit = QtWidgets.QLineEdit(self)
        hLayout.addWidget(self.lineEdit)

#        self.filter = QtWidgets.QPushButton("Search", self)
#        hLayout.addWidget(self.filter)
#        self.filter.clicked.connect(self.filterClicked)
        self.lineEdit.textChanged.connect(self.filterClicked)

        self.radio_group=QtWidgets.QButtonGroup(hLayout) # Number group
        self.button1  = QtWidgets.QRadioButton('HARPS RVBank', self)
        self.button2  = QtWidgets.QRadioButton('HIRES NZP', self)
        self.button1.setChecked(True)
#        self.button2.setEnabled(False)
#        self.button3  = QtWidgets.QRadioButton('symbol "t1"', self)

        self.radio_group.addButton(self.button1)
        self.radio_group.addButton(self.button2)

        hLayout.addWidget(self.button1)
        hLayout.addWidget(self.button2)

        self.radio_group.buttonClicked.connect(self.init_model)


        self.readme_button = QtWidgets.QPushButton('READ ME', self)
        hLayout.addWidget(self.readme_button)

        self.readme_button.clicked.connect(self.info)

        #self.radio_group.buttonClicked.connect(self.final_output)
        #self.button2.toggled.connect(self.final_output)
########  Targets ##############
        self.list = QtWidgets.QListView(self)
        vLayout.addLayout(hLayout)
        vLayout.addWidget(self.list)
        self.model = QtGui.QStandardItemModel(self.list)


########  Options ##############

        self.list_opt = QtWidgets.QListView(self)
        vLayout.addWidget(self.list_opt)
        self.model_opt = QtGui.QStandardItemModel(self.list_opt)

        self.setGeometry(2,2, 555, 525) 
        
        
        self.list.clicked.connect(self.on_clicked)
        self.list_opt.clicked.connect(self.on_clicked_opt)

        #self.list.setSelectionRectVisible(True)
        #print(self.list.isSelectionRectVisible())

        url = "http://www2.mpia.de/homes/trifonov/%s_RVs/%s_harps_all-data_v1.dat"%(targets_HARPS[0],targets_HARPS[0])
        self.path = url
        self.data_index  = 1
        self.row_opt = 0
        self.row = 0
        self.type_data = "HARPS"
        
        self.try_connection(url)
        if self.url_success == False:
            return

        try:
            self.x_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[0])
            self.y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[1])
            self.e_y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[2])

        except:
            print("Something is wrong.... No connection to the HARPS RVBank!!!")
            self.url_success = False            

        self.init_model()


    def try_connection(self,url):

        try:
            self.resp = urlopen(url).read() 
            self.url_success = True
            return
        except:
            print("No internet connection, or no connection to the HARPS RVBank!!!")
            self.url_success = False

 
    def init_model(self):


        self.model.clear()
        self.model_opt.clear()

        if self.button1.isChecked():
            targets    = targets_HARPS
            data_files = data_files_HARPS
            data_files_ind = data_files_ind_HARPS
        elif self.button2.isChecked():
            targets    = targets_HIRES
            data_files = data_files_HIRES
            data_files_ind = data_files_ind_HIRES

        for code in targets:
            item = QtGui.QStandardItem(code)
            item.setCheckable(False)
            item.setEditable(False)
            self.model.appendRow(item)
        self.list.setModel(self.model)

        for code in data_files:
            item = QtGui.QStandardItem(code)
            item.setCheckable(False)
            item.setEditable(False)

            self.model_opt.appendRow(item)
        self.list_opt.setModel(self.model_opt)

        self.filterClicked()



#####################################
#        self.list.selectionModel().currentChanged.connect(self.on_row_changed)
#
#    def on_row_changed(self, current, previous):
#        self._INDEX = current.row()
#        print('Row %d selected' % current.row())
#####################################



    def filterClicked(self):

        filter_text = str(self.lineEdit.text()).lower()
        valid_row = []
        for row in range(self.model.rowCount()):
            if filter_text in str(self.model.item(row).text()).lower():
                self.list.setRowHidden(row, False)
                valid_row.append(row)
            else:
                self.list.setRowHidden(row, True)

        try:
            ix = self.model.index(valid_row[0], 0)
        except:
            ix = self.model.index(0, 0)   
         
        sm = self.list.selectionModel()
        sm.select(ix, QtCore.QItemSelectionModel.Select)
        
#        index = self.model.indexFromItem(self.model.item(valid_row[-1]))
 
        #self.on_clicked(ix)
        self.row = ix.row()
        self.final_output()

    def on_clicked(self, index):

        self.row = index.row()
        self.final_output()

    def on_clicked_opt(self, index):
        self.row_opt = index.row()
        self.final_output()
        
    def final_output(self):
        
        row =self.row
        row_opt = self.row_opt
        
        #print(row,row_opt)
        
        
        if self.button1.isChecked():

            self.type_data = "HARPS"
            self.data_index = data_files_ind_HARPS[row_opt]
            self.data_name = data_files_HARPS[row_opt]
            
            url = "http://www2.mpia.de/homes/trifonov/%s_RVs/%s_harps_all-data_v1.dat"%(targets_HARPS[row],targets_HARPS[row])
           # resp = urlopen(url).read() 
            
            self.try_connection(url)
            if self.url_success == False:
                return
            
            self.x_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[0])
            if self.data_index <22:
                self.y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[self.data_index])
                self.e_y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[self.data_index+1])
            elif self.data_index == 22 or self.data_index == 24:
                self.y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[self.data_index])* 1000.0
                self.e_y_data = np.array([1.0]*len(self.y_data))
            else:
                self.y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[self.data_index])
                self.e_y_data = np.array([np.mean(self.y_data)*0.01]*len(self.y_data))

            self.path = url
            self.target_name = targets_HARPS[row]


        elif self.button2.isChecked():

            self.type_data = "HIRES"
            
            if row_opt > 3:
                row_opt = 0
            
            self.data_index = data_files_ind_HIRES[row_opt]
            self.data_name = data_files_HIRES[row_opt]
            
            #print(self.data_index,self.data_name)

            url = "http://www2.mpia.de/homes/trifonov/HIRES/%s_RVs/%s.dat"%(targets_HIRES[row],targets_HIRES[row])
            #resp = urlopen(url).read() 
            
            self.try_connection(url)
            if self.url_success == False:
                return
            
            self.x_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[0])
            if self.data_index <5:
                self.y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[self.data_index])
                self.e_y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[self.data_index+1])
            elif self.data_index >= 5:
                self.y_data = np.genfromtxt(io.BytesIO(self.resp),usecols=[self.data_index])
                self.e_y_data = np.array([0.001]*len(self.y_data))

            self.path = url
            self.target_name = targets_HIRES[row]



    def info(self):
        
       # self.info_dialog.setGeometry(300, 300, 150, 150)
        self.info_dialog.setFixedSize(550, 700)

        self.info_dialog.setWindowTitle('RVBank info')
 
    
        text = ''
        self.info_dialog.text.setText(text) 
        
        text = """
<br>A few things you need to know about the RVBank data:
<br>
<br>
<br>The data are retrieved via an internet connection from the RVBank web pages: 
<br>* <a href='http://www.mpia.de/homes/trifonov/HARPS_RVBank.html'>HARPS RVBank</a>, and
<br>* <a href='http://www.mpia.de/homes/trifonov/HIRES_RVBank.html'>HIRES RVBank</a>
<br>
<br>Therefore, you need an internet connection for browsing the RVBank data.
"""
        self.info_dialog.text.append(text)


        text = """
<br>
<br>
If you made the use of the HARPS data for your paper, please cite:
<br>* <a href='https://ui.adsabs.harvard.edu/abs/2020arXiv200105942T/abstract'> Trifonov et al. (2020)</a>
<br>(also, please check for relevant references in the paper).
<br>
<br>If you made the use of the original HIRES data, please cite:
<br>* <a href='https://ui.adsabs.harvard.edu/abs/2017AJ....153..208B/abstract'> Butler et al. (2017)</a>
<br>
<br>If you made the use of the NZP corrected HIRES data, please cite:
<br>* <a href='https://ui.adsabs.harvard.edu/abs/2019MNRAS.484L...8T/abstract'> Tal-Or et al. (2019)</a>
together with  <a href='https://ui.adsabs.harvard.edu/abs/2017AJ....153..208B/abstract'> Butler et al. (2017)</a>  
<br>
<br>
<br>
<br> Remember! 
<br> 
<br> * HARPS-DRS activity data such as FWHM, Contrast, and Bisector-span, originally do not have error bars! Thus, the Exo-Striker adopts "sigma=1 m/s" for CCF and FWHM data, and "sigma=0.01" for the Contrast data. 
<br> 
<br> * HIRES s- and h-index data originally do not have error bars! Thus, the Exo-Striker adopts "sigma"=0.001.
<br> 
<br> * Some HIRES h-index data appear as -1.0. This means that for the given epoch the h-index data was not possible to be computed.
"""

 
        self.info_dialog.text.append(text)

    
        self.info_dialog.text.setReadOnly(True)
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))
        self.info_dialog.show()










data_files_HARPS = ["RVs SERVAL + NZP correction",
              "RVs DRS + NZP correction",
              "RVs SERVAL",
              "RVs DRS",
              "CRX","dLW","Halpha","NaD1","NaD2","FWHM_DRS","CONTRAST_DRS","Bisector"]

data_files_ind_HARPS = [1,3,5,7,11,13,15,17,19,22,23,24]


targets_HARPS = ['BD+053829', 'BD+053835', 'BD+072351', 'BD+072474', 'BD+101802', 'BD+144559', 'BD+20594', 'BD+23527', 'BD-013943', 'BD-053461', 'BD-055432', 'BD-073537', 'BD-07436B', 'BD-096003', 'BD-114672', 'BD-122361', 'BD-122371', 'BD-132130', 'BD-136424', 'BD-156276', 'BD-15705', 'BD-17751', 'BD-184705', 'BD-184710', 'BD-194341', 'BD-202976', 'BD-221092', 'CD-2312854', 'CD-2412030', 'CD-2510553', 'CD-2514757', 'CD-261407', 'CD-281745', 'CD-287411', 'CD-301812', 'CD-3019757', 'CD-306011', 'CD-306015', 'CD-312415', 'CD-329927', 'CD-332771', 'CD-337795B', 'CD-348618', 'CD-3510525B', 'CD-381059', 'CD-383745', 'CD-383754', 'CD-386968', 'CD-395833', 'CD-426963', 'CD-452997', 'CD-461064', 'CD-4911401', 'CD-4911402', 'CD-4911404', 'CD-4911410', 'CD-4911415', 'CD-50714', 'CD-52381', 'CD-535815', 'CD-535820', 'CD-535829', 'CD-536474', 'CD-546002', 'CD-546003', 'CD-551709', 'CD-567708', 'CD-584411', 'CD-60370', 'CD-692101', 'CD-8028', 'CD-8480', 'CoRoT-13', 'CoRoT-1', 'CoRoT-22', 'CoRoT-24', 'CoRoT-25', 'CoRoT-32', 'CoRoT-7', 'CoRoT101332685', 'CoRoT101614469', 'CoRoT102296945', 'CoRoT102352354', 'CoRoT102387834', 'CoRoT102568803', 'CoRoT102615551', 'CoRoT102755837', 'CoRoT104620638', 'CoRoT110676555', 'CoRoT110752102', 'CoRoT110838244', 'CoRoT110858446', 'CoRoT221699621', 'CoRoT604178318', 'CoRoT631423929', 'EPIC210558622', 'EPIC210744674', 'EPIC210968143', 'EPIC212460519', 'EPIC212469831', 'EPIC212572439', 'EPIC212645891', 'EPIC212735333', 'EPIC219409830', 'EPIC219503117', 'EPIC219661601', 'EPIC219722212', 'GJ1001', 'GJ1002', 'GJ1008', 'GJ1009', 'GJ1012', 'GJ1020', 'GJ1021', 'GJ1022', 'GJ1032', 'GJ1043', 'GJ1044', 'GJ1046', 'GJ1048', 'GJ1050', 'GJ1056', 'GJ1057', 'GJ105A', 'GJ105B', 'GJ1061', 'GJ1065', 'GJ1066', 'GJ1068', 'GJ1075', 'GJ1079', 'GJ1084', 'GJ1085', 'GJ1089', 'GJ108', 'GJ1094', 'GJ1099', 'GJ10', 'GJ1123', 'GJ1126', 'GJ1129', 'GJ1132', 'GJ1135', 'GJ1137', 'GJ1144', 'GJ1145', 'GJ115', 'GJ1161', 'GJ1177', 'GJ117', 'GJ118', 'GJ1191', 'GJ1214', 'GJ121', 'GJ1224', 'GJ1232', 'GJ1236', 'GJ1239', 'GJ1247', 'GJ1248', 'GJ1256', 'GJ1265', 'GJ126', 'GJ1279', 'GJ1280', 'GJ1282', 'GJ128', 'GJ12', 'GJ131', 'GJ132', 'GJ135', 'GJ136', 'GJ138', 'GJ143.1', 'GJ143', 'GJ146', 'GJ149', 'GJ150.0', 'GJ159', 'GJ162.2', 'GJ163', 'GJ166C', 'GJ167.2', 'GJ16', 'GJ173', 'GJ177.1', 'GJ177', 'GJ178', 'GJ179', 'GJ17', 'GJ18.0', 'GJ180', 'GJ182', 'GJ183', 'GJ187', 'GJ188', 'GJ189.2', 'GJ191', 'GJ19', 'GJ1', 'GJ2001', 'GJ2018', 'GJ2030', 'GJ2033', 'GJ2037', 'GJ204.2', 'GJ2046', 'GJ2049', 'GJ204', 'GJ2056', 'GJ2058', 'GJ205', 'GJ2060', 'GJ2066', 'GJ2085', 'GJ208', 'GJ2119', 'GJ2126', 'GJ2127', 'GJ2128', 'GJ2137', 'GJ213', 'GJ214', 'GJ218', 'GJ219.0', 'GJ221', 'GJ223.3', 'GJ224', 'GJ225.1', 'GJ225', 'GJ227', 'GJ230', 'GJ233', 'GJ234', 'GJ236', 'GJ238', 'GJ243', 'GJ245.1', 'GJ250', 'GJ253', 'GJ260', 'GJ268.2', 'GJ273', 'GJ275', 'GJ279', 'GJ280', 'GJ281', 'GJ285', 'GJ288', 'GJ291.1', 'GJ293.1', 'GJ294A', 'GJ294B', 'GJ298', 'GJ299', 'GJ29', 'GJ300', 'GJ3018', 'GJ3021', 'GJ302', 'GJ3049', 'GJ3053', 'GJ3071', 'GJ3072', 'GJ3082', 'GJ3091', 'GJ309', 'GJ31.3', 'GJ3102', 'GJ3110', 'GJ312', 'GJ3138', 'GJ313', 'GJ3141', 'GJ3148', 'GJ3161', 'GJ317', 'GJ3187', 'GJ3189', 'GJ3193', 'GJ31', 'GJ3200', 'GJ3202', 'GJ3205', 'GJ3207', 'GJ321.3', 'GJ3212', 'GJ3218', 'GJ3221', 'GJ3222', 'GJ3256', 'GJ3258', 'GJ3260', 'GJ3263', 'GJ3279', 'GJ327', 'GJ3293', 'GJ3300', 'GJ3305', 'GJ3307', 'GJ3313', 'GJ3314', 'GJ3317', 'GJ3323', 'GJ3325', 'GJ3328', 'GJ3332', 'GJ333', 'GJ3340', 'GJ3341', 'GJ334', 'GJ3356', 'GJ3367', 'GJ3376', 'GJ3379', 'GJ337', 'GJ33', 'GJ3404', 'GJ341', 'GJ3436', 'GJ3446', 'GJ3449', 'GJ3455', 'GJ3470', 'GJ3492', 'GJ3493', 'GJ3502', 'GJ3508', 'GJ3519', 'GJ3530', 'GJ3535', 'GJ3543', 'GJ355', 'GJ3563', 'GJ357', 'GJ358', 'GJ3591', 'GJ3598', 'GJ3603', 'GJ3615', 'GJ3618', 'GJ361', 'GJ3620', 'GJ3633', 'GJ3634', 'GJ364.0', 'GJ3643', 'GJ3644', 'GJ3648', 'GJ3663', 'GJ3671', 'GJ3677', 'GJ367', 'GJ369', 'GJ3700', 'GJ3701', 'GJ3707', 'GJ3708', 'GJ3727', 'GJ3737', 'GJ3746', 'GJ3752', 'GJ3759', 'GJ3769', 'GJ3778', 'GJ377', 'GJ379.1', 'GJ3792', 'GJ3796', 'GJ3813', 'GJ3817', 'GJ3822', 'GJ3827', 'GJ382', 'GJ3830', 'GJ3835', 'GJ3838', 'GJ386', 'GJ3871', 'GJ3874', 'GJ3885', 'GJ388', 'GJ3896', 'GJ3900', 'GJ3905', 'GJ390', 'GJ3915', 'GJ3918', 'GJ393', 'GJ3944', 'GJ3952', 'GJ3962', 'GJ3969', 'GJ3970', 'GJ398.1', 'GJ3987', 'GJ3996', 'GJ3998', 'GJ399', 'GJ3', 'GJ4001', 'GJ4008', 'GJ4014', 'GJ401', 'GJ402.1', 'GJ402', 'GJ404', 'GJ4052', 'GJ4056', 'GJ406', 'GJ4071', 'GJ4077', 'GJ4088', 'GJ4092', 'GJ40', 'GJ4100', 'GJ410', 'GJ412.3', 'GJ413.1', 'GJ4130', 'GJ4140', 'GJ4154', 'GJ416.1', 'GJ4160', 'GJ419', 'GJ4206', 'GJ421', 'GJ423.1', 'GJ4244', 'GJ4248', 'GJ4291', 'GJ4303', 'GJ4306', 'GJ4322', 'GJ432', 'GJ433', 'GJ4340', 'GJ4347', 'GJ4353', 'GJ435', 'GJ4364', 'GJ436', 'GJ4377', 'GJ4383', 'GJ438', 'GJ442', 'GJ446', 'GJ447', 'GJ448.0', 'GJ449', 'GJ452.1', 'GJ452', 'GJ45', 'GJ465', 'GJ468', 'GJ472', 'GJ476', 'GJ477', 'GJ478', 'GJ479', 'GJ480', 'GJ488', 'GJ493', 'GJ496.1', 'GJ503.3', 'GJ506.1', 'GJ506', 'GJ508.3', 'GJ510', 'GJ511', 'GJ513', 'GJ514', 'GJ517', 'GJ52.1', 'GJ521.1', 'GJ522', 'GJ524.1', 'GJ526', 'GJ529', 'GJ530', 'GJ534', 'GJ536', 'GJ538.1', 'GJ54.1', 'GJ541.1', 'GJ542', 'GJ544', 'GJ550', 'GJ551', 'GJ552', 'GJ553', 'GJ557.0', 'GJ558.1', 'GJ559A', 'GJ559B', 'GJ55', 'GJ56.1', 'GJ560', 'GJ563.2', 'GJ565', 'GJ569', 'GJ570.1', 'GJ570A', 'GJ570B', 'GJ571', 'GJ574', 'GJ579.2', 'GJ579.3', 'GJ57', 'GJ58.2', 'GJ581', 'GJ583', 'GJ587', 'GJ588', 'GJ58', 'GJ594', 'GJ599', 'GJ603', 'GJ604', 'GJ606', 'GJ611.3', 'GJ613.0', 'GJ615', 'GJ616', 'GJ618.1', 'GJ618.4', 'GJ620', 'GJ622', 'GJ624', 'GJ628', 'GJ629.3', 'GJ62', 'GJ631', 'GJ634.1', 'GJ634', 'GJ637', 'GJ641', 'GJ643', 'GJ645', 'GJ65.2', 'GJ650', 'GJ653', 'GJ654', 'GJ656', 'GJ657', 'GJ660.1', 'GJ660', 'GJ664', 'GJ666', 'GJ667', 'GJ66', 'GJ674', 'GJ676', 'GJ680', 'GJ683', 'GJ686', 'GJ688', 'GJ691', 'GJ692', 'GJ693', 'GJ695.1', 'GJ696', 'GJ698.1', 'GJ699.2', 'GJ699', 'GJ701', 'GJ702', 'GJ707', 'GJ70', 'GJ710', 'GJ711', 'GJ715', 'GJ722', 'GJ724.2', 'GJ724', 'GJ726', 'GJ729', 'GJ739', 'GJ740', 'GJ744', 'GJ747.4', 'GJ752', 'GJ754', 'GJ755', 'GJ762', 'GJ768.1', 'GJ770', 'GJ771', 'GJ776', 'GJ780', 'GJ781.2', 'GJ783', 'GJ784.1', 'GJ784', 'GJ785', 'GJ787', 'GJ788.1', 'GJ792.1', 'GJ796', 'GJ797', 'GJ798', 'GJ79', 'GJ7', 'GJ800', 'GJ801', 'GJ803', 'GJ805', 'GJ808', 'GJ811', 'GJ812.1', 'GJ812', 'GJ817', 'GJ818', 'GJ821', 'GJ825.2', 'GJ825', 'GJ826.1', 'GJ827.1', 'GJ827', 'GJ83.1', 'GJ836.7', 'GJ838', 'GJ83', 'GJ84.1', 'GJ841', 'GJ842', 'GJ845', 'GJ846', 'GJ848', 'GJ849.1', 'GJ849', 'GJ851.3', 'GJ853', 'GJ855.1', 'GJ855', 'GJ85', 'GJ862', 'GJ864', 'GJ86', 'GJ875', 'GJ876', 'GJ877', 'GJ879', 'GJ87', 'GJ880', 'GJ881', 'GJ882', 'GJ887', 'GJ891', 'GJ894', 'GJ9003', 'GJ9007', 'GJ9008', 'GJ900', 'GJ9010', 'GJ9012', 'GJ9015', 'GJ9018', 'GJ9024', 'GJ9027', 'GJ9029', 'GJ902', 'GJ9031', 'GJ904', 'GJ9055', 'GJ9056', 'GJ9067', 'GJ9075', 'GJ9079', 'GJ9080', 'GJ908', 'GJ9103A', 'GJ9107', 'GJ9109', 'GJ9116', 'GJ9118', 'GJ9121', 'GJ912', 'GJ9133', 'GJ9144', 'GJ9156', 'GJ9163', 'GJ9165', 'GJ9177', 'GJ9178', 'GJ9190', 'GJ9191', 'GJ9201', 'GJ9210', 'GJ9219', 'GJ9222', 'GJ9223', 'GJ9241', 'GJ9252', 'GJ9254', 'GJ9255', 'GJ9263', 'GJ9268', 'GJ9269', 'GJ9273', 'GJ9283', 'GJ9296', 'GJ9299', 'GJ9316', 'GJ9320', 'GJ9322', 'GJ9324', 'GJ9328', 'GJ9332', 'GJ9338', 'GJ9356', 'GJ9367', 'GJ9384', 'GJ9396', 'GJ9411', 'GJ9415', 'GJ9423', 'GJ9429', 'GJ9441', 'GJ9467', 'GJ9472', 'GJ9473', 'GJ9476', 'GJ9479', 'GJ9481', 'GJ9482', 'GJ9498', 'GJ9507', 'GJ9527', 'GJ9536', 'GJ9547', 'GJ9563', 'GJ9565', 'GJ9583', 'GJ9590', 'GJ9592', 'GJ95', 'GJ9626', 'GJ9630', 'GJ9646', 'GJ9655', 'GJ9659', 'GJ9663', 'GJ9673', 'GJ9700', 'GJ9714', 'GJ9732', 'GJ9733', 'GJ9762', 'GJ9765', 'GJ9771', 'GJ9774', 'GJ9787', 'GJ9789', 'GJ9790', 'GJ9797', 'GJ97', 'GJ9801', 'GJ9803', 'GJ9824', 'GJ9827', 'GJ9828', 'GJ9829', 'GJ9833', 'GJ9836', 'GJ9847', 'HATS-17', 'HATS-31', 'HATS411-024', 'HATS411-062', 'HATS412-016', 'HATS511-001', 'HATS524-001', 'HATS535-001', 'HATS537-009', 'HATS538-002', 'HATS547-004', 'HATS549-008', 'HATS551-012', 'HATS551-017', 'HATS554-051', 'HATS554-059', 'HATS555-001', 'HATS557-026', 'HATS557-030', 'HATS557-037', 'HATS560-012', 'HATS561-002', 'HATS563-014', 'HATS565-001', 'HATS565-005', 'HATS566-006', 'HATS586-003', 'HATS598-010', 'HATS598-013', 'HATS600-004', 'HATS600-027', 'HATS601-056', 'HATS601-059', 'HATS601-061', 'HATS602-002', 'HATS602-003', 'HATS602-035', 'HATS602-040', 'HATS602-072', 'HATS602-124', 'HATS606-001', 'HATS606-010', 'HATS610-007', 'HATS625-003', 'HATS625-020', 'HATS625-025', 'HATS647-005', 'HATS698-015', 'HATS698-041', 'HATS699-039', 'HATS699-049', 'HATS699-050', 'HATS700-003', 'HATS700-026', 'HATS709-003', 'HATS737-002', 'HATS745-004', 'HATS745-015', 'HATS745-021', 'HATS746-009', 'HATS747-013', 'HATS755-007', 'HATS756-008', 'HATS762-007', 'HATS762-014', 'HATS772-088', 'HATS772-095', 'HATS772-126', 'HATS772-143', 'HATS772-186', 'HATS778-004', 'HATS778-015', 'HATS778-018', 'HATS778-024', 'HD100069', 'HD10008', 'HD100289', 'HD100363', 'HD100407', 'HD10042', 'HD100508', 'HD100655', 'HD10069', 'HD100777', 'HD100889', 'HD101063', 'HD101065', 'HD101197', 'HD101339', 'HD101348', 'HD101367', 'HD101412', 'HD101612', 'HD101615', 'HD101644', 'HD101650', 'HD10180', 'HD101886', 'HD10188', 'HD101930', 'HD102117', 'HD102124', 'HD102188', 'HD102195', 'HD102196', 'HD102200', 'HD102272', 'HD102300', 'HD102346', 'HD102361', 'HD102458', 'HD10278', 'HD102843', 'HD103197', 'HD103234', 'HD10337', 'HD103589', 'HD103720', 'HD103743', 'HD103774', 'HD10380', 'HD103847', 'HD103891', 'HD103949', 'HD104006', 'HD104067', 'HD1040', 'HD104181', 'HD104263', 'HD10472', 'HD104760', 'HD104800', 'HD104819', 'HD104982', 'HD105575', 'HD105665', 'HD105690', 'HD105707', 'HD105779', 'HD105837', 'HD105850', 'HD105938', 'HD105', 'HD10615', 'HD106275', 'HD106290', 'HD106315', 'HD106444', 'HD10647', 'HD106589', 'HD106661', 'HD106906', 'HD106937', 'HD10700', 'HD107148', 'HD107181', 'HD107383', 'HD10761', 'HD108063', 'HD108147', 'HD108177', 'HD108202', 'HD108309', 'HD10830', 'HD108341', 'HD108564', 'HD108570', 'HD108767B', 'HD108767', 'HD108768', 'HD108857', 'HD108935', 'HD108953', 'HD10895', 'HD109098', 'HD109138', 'HD109271', 'HD109310', 'HD109379', 'HD109408', 'HD109409', 'HD109536', 'HD109573', 'HD109684', 'HD109787', 'HD109988', 'HD110014', 'HD110058', 'HD11025', 'HD110291', 'HD110411', 'HD110537', 'HD110557', 'HD110619', 'HD110621', 'HD110668', 'HD111133', 'HD111232', 'HD111261B', 'HD111519', 'HD111564', 'HD11168', 'HD111777', 'HD111814', 'HD11195', 'HD111998', 'HD112100', 'HD112109', 'HD11226', 'HD11231', 'HD112540', 'HD11262', 'HD112646', 'HD11271', 'HD112934', 'HD113226', 'HD11332', 'HD113449', 'HD113569', 'HD113679', 'HD113904', 'HD11397', 'HD113998', 'HD114076', 'HD114330', 'HD114386', 'HD114561', 'HD114613', 'HD114642', 'HD114729', 'HD114747', 'HD114853', 'HD115031', 'HD11505', 'HD115169', 'HD115231', 'HD115341', 'HD115496', 'HD115585', 'HD115659', 'HD115674', 'HD115773', 'HD115902', 'HD116087', 'HD11608', 'HD116114', 'HD116160', 'HD11616', 'HD116235', 'HD116259', 'HD116284', 'HD116410', 'HD116434', 'HD116458', 'HD116568', 'HD116883', 'HD116941', 'HD117105', 'HD117126', 'HD117190', 'HD1171', 'HD117207', 'HD117359', 'HD117618', 'HD118022', 'HD118466', 'HD118551', 'HD118563', 'HD118679', 'HD119027', 'HD119173', 'HD119288', 'HD119329', 'HD11938', 'HD119503', 'HD119537', 'HD119629', 'HD119638', 'HD11964', 'HD11977', 'HD119782', 'HD119802', 'HD119949', 'HD120250', 'HD120344', 'HD120362', 'HD12039', 'HD12055', 'HD120734', 'HD121004', 'HD121416', 'HD121504', 'HD122173', 'HD122194', 'HD122308', 'HD122430', 'HD122474', 'HD122603', 'HD122640', 'HD12264', 'HD12270', 'HD122970', 'HD123265', 'HD123319', 'HD12345', 'HD123517', 'HD123619', 'HD123651', 'HD123732', 'HD123767', 'HD12387', 'HD123999', 'HD124244', 'HD12438', 'HD124523', 'HD124785', 'HD124900', 'HD124931', 'HD12524', 'HD125271', 'HD125522', 'HD125612B', 'HD125823', 'HD12583', 'HD125881', 'HD12617', 'HD126248', 'HD126515', 'HD126535', 'HD126583', 'HD126661', 'HD126793', 'HD126803', 'HD126829', 'HD126999', 'HD127124', 'HD127423', 'HD128207', 'HD128356', 'HD128429', 'HD128431', 'HD128571', 'HD128674', 'HD128760', 'HD129191', 'HD129229', 'HD12932', 'HD129422', 'HD129642', 'HD129685', 'HD129704', 'HD129750', 'HD129814', 'HD129899', 'HD129926', 'HD129956', 'HD130109', 'HD130322', 'HD13060', 'HD130694', 'HD130930', 'HD130952', 'HD130989', 'HD131218', 'HD131399', 'HD13147', 'HD131653', 'HD131664', 'HD13183', 'HD132052', 'HD132307', 'HD13252', 'HD132569', 'HD13263', 'HD132648', 'HD132944', 'HD133112', 'HD133469', 'HD13350', 'HD133574', 'HD13357', 'HD133600', 'HD133633', 'HD133813', 'HD134060', 'HD134088', 'HD134113', 'HD134214', 'HD13423', 'HD134305', 'HD134330', 'HD134353', 'HD134481', 'HD134606', 'HD134664', 'HD134702', 'HD134888', 'HD134929', 'HD134987', 'HD13513A', 'HD135468', 'HD135559', 'HD135592', 'HD135625', 'HD135778', 'HD13578', 'HD135953', 'HD136352', 'HD136894', 'HD13692', 'HD137010', 'HD137057', 'HD13724', 'HD137388', 'HD13789', 'HD137909', 'HD13808', 'HD138403', 'HD138549', 'HD138573', 'HD138716', 'HD138763', 'HD138799', 'HD138914', 'HD13904', 'HD139084', 'HD139129', 'HD139189', 'HD139211', 'HD139271', 'HD139332', 'HD13940', 'HD139536', 'HD139590', 'HD139710', 'HD139763', 'HD139879', 'HD13997', 'HD140637', 'HD140785', 'HD141011', 'HD141024', 'HD141128', 'HD14127', 'HD141344', 'HD141513', 'HD141597', 'HD141598', 'HD141624', 'HD141637', 'HD141851', 'HD141937', 'HD141943', 'HD141969', 'HD142070', 'HD142072', 'HD142139', 'HD142331', 'HD142527', 'HD142629', 'HD142630', 'HD142879', 'HD142921', 'HD142A', 'HD143114', 'HD143361', 'HD143436', 'HD143638', 'HD14374', 'HD143790', 'HD144059', 'HD144411', 'HD144497', 'HD14452', 'HD144550', 'HD144585', 'HD144589', 'HD144608', 'HD144846', 'HD144848', 'HD144880', 'HD144897', 'HD144899', 'HD145005', 'HD145344', 'HD145377', 'HD145598', 'HD145666', 'HD145689', 'HD145927', 'HD1461', 'HD146296', 'HD146389', 'HD146464', 'HD146514', 'HD14651', 'HD146546', 'HD146624', 'HD1466', 'HD14703', 'HD147135', 'HD147195', 'HD147386', 'HD14744', 'HD14745', 'HD14747', 'HD147512', 'HD147518', 'HD147642', 'HD147644', 'HD147787', 'HD147873', 'HD147935', 'HD148156', 'HD148184', 'HD148211', 'HD14830', 'HD14832', 'HD148427', 'HD148577', 'HD14868', 'HD148816', 'HD148998', 'HD149143', 'HD149189', 'HD149200', 'HD149277', 'HD149396', 'HD149438', 'HD14943', 'HD149724', 'HD149747', 'HD149782', 'HD149933', 'HD150139', 'HD150177', 'HD150248', 'HD150437', 'HD150562', 'HD150633', 'HD150798', 'HD150936', 'HD15115', 'HD151504', 'HD151692', 'HD151772', 'HD152079', 'HD1522', 'HD152433', 'HD152533', 'HD152555', 'HD153053', 'HD15318', 'HD153276', 'HD153363', 'HD15337', 'HD153851', 'HD153950', 'HD154088', 'HD154160', 'HD154195', 'HD154221', 'HD154387', 'HD154672', 'HD154708', 'HD154962', 'HD15507', 'HD155245', 'HD155712', 'HD155806', 'HD155968', 'HD156079', 'HD156098', 'HD15612', 'HD156411', 'HD156423', 'HD156751', 'HD156846', 'HD156854', 'HD156991', 'HD157060', 'HD157172', 'HD157338', 'HD157347', 'HD157527', 'HD157751', 'HD157830', 'HD157950', 'HD157999', 'HD158148', 'HD158352', 'HD158469', 'HD158809', 'HD15906', 'HD159170', 'HD159194', 'HD159376', 'HD15942', 'HD160013', 'HD160089', 'HD16008', 'HD16031', 'HD160613', 'HD160617', 'HD160635', 'HD160720', 'HD160836', 'HD160839', 'HD161098', 'HD161256', 'HD16141', 'HD161459', 'HD161555', 'HD161566', 'HD161592', 'HD161868', 'HD16189', 'HD16193', 'HD162049', 'HD162236', 'HD162396', 'HD162598', 'HD16280', 'HD16297', 'HD163441', 'HD16382', 'HD16400', 'HD16417', 'HD164233', 'HD164249A', 'HD164509', 'HD165131', 'HD165135', 'HD165155', 'HD165204', 'HD16522', 'HD165474', 'HD16548', 'HD165760', 'HD165920', 'HD166473', 'HD166599', 'HD166724', 'HD166745', 'HD167060', 'HD167096', 'HD16714', 'HD167300', 'HD167359', 'HD167468', 'HD167554', 'HD16765', 'HD167677', 'HD167858', 'HD16815', 'HD168432', 'HD168525', 'HD168746', 'HD168769', 'HD168863', 'HD168871', 'HD168960', 'HD16905', 'HD1690', 'HD169142', 'HD169178', 'HD169370', 'HD169690', 'HD16975', 'HD169822', 'HD169830', 'HD169836', 'HD170031', 'HD170053', 'HD170174', 'HD170200', 'HD170231', 'HD170580', 'HD170699', 'HD170706', 'HD170958', 'HD170973', 'HD171028', 'HD171488', 'HD171587', 'HD171759', 'HD171834', 'HD171942', 'HD172046', 'HD172103', 'HD172211', 'HD172513', 'HD172555', 'HD172568', 'HD172643', 'HD173282', 'HD173427', 'HD173540', 'HD17374', 'HD173787', 'HD1737', 'HD173885', 'HD174153', 'HD174295', 'HD17439', 'HD174429', 'HD174512', 'HD174545', 'HD175128', 'HD175219', 'HD17548', 'HD175541', 'HD175607', 'HD175639', 'HD17599', 'HD176354', 'HD176387', 'HD176535', 'HD176638', 'HD176986', 'HD177033', 'HD177122', 'HD177178', 'HD177409', 'HD177474', 'HD177756', 'HD177758', 'HD17793', 'HD17824', 'HD178340', 'HD178484', 'HD17848', 'HD178787', 'HD178904', 'HD179079', 'HD179346', 'HD179433', 'HD179949', 'HD18001', 'HD180134', 'HD180263', 'HD180409', 'HD18083', 'HD181249', 'HD181296', 'HD181327', 'HD18134', 'HD181387', 'HD181428', 'HD181433', 'HD181517', 'HD181544', 'HD181720', 'HD181907', 'HD182893', 'HD18293', 'HD18322', 'HD183263', 'HD18330', 'HD183414', 'HD183579', 'HD183826', 'HD18386', 'HD183877', 'HD183993', 'HD184552', 'HD184588', 'HD184711', 'HD184768', 'HD185256', 'HD185283', 'HD185615', 'HD185679', 'HD18599', 'HD186194', 'HD186265', 'HD186302', 'HD18632', 'HD186704', 'HD186803', 'HD186843', 'HD18692', 'HD18702', 'HD187085', 'HD18708', 'HD18719', 'HD187532', 'HD18754', 'HD187669', 'HD18777', 'HD188042', 'HD188091', 'HD188114', 'HD188228', 'HD188298', 'HD188345', 'HD188748', 'HD188815', 'HD18885', 'HD189005', 'HD189242', 'HD189625', 'HD189627', 'HD189631', 'HD189987', 'HD190125', 'HD190204', 'HD190290', 'HD190524', 'HD190590', 'HD190594', 'HD190613', 'HD190647', 'HD190954', 'HD190984', 'HD1909', 'HD191033', 'HD19107', 'HD191122', 'HD191760B', 'HD192031', 'HD192263', 'HD19230', 'HD192425', 'HD193334', 'HD193406', 'HD193690', 'HD193721', 'HD193756', 'HD193826', 'HD193844', 'HD193901', 'HD193995', 'HD19423', 'HD194490', 'HD194717', 'HD19493', 'HD195019', 'HD195145', 'HD19518', 'HD195200', 'HD195302', 'HD19545', 'HD195633', 'HD195961', 'HD196050', 'HD196286', 'HD196378', 'HD196384', 'HD196385', 'HD196390', 'HD196397', 'HD19641', 'HD196470', 'HD19668', 'HD196724', 'HD196800', 'HD196892', 'HD197027', 'HD197051', 'HD197083', 'HD197197', 'HD197210', 'HD197286', 'HD197300', 'HD197536', 'HD197557', 'HD197623', 'HD19773', 'HD197746', 'HD197777', 'HD197818', 'HD197823', 'HD197890', 'HD197921', 'HD197952', 'HD1979', 'HD198075', 'HD198227', 'HD198232', 'HD19827', 'HD198390', 'HD198748', 'HD199086', 'HD19918', 'HD199190', 'HD199254', 'HD199289', 'HD199604', 'HD199665', 'HD199847', 'HD199868', 'HD199960', 'HD20003', 'HD200078', 'HD200083', 'HD20010', 'HD200143', 'HD200349', 'HD20037', 'HD200538', 'HD200565', 'HD200633', 'HD200761', 'HD200763', 'HD200869', 'HD200962', 'HD201161', 'HD201219', 'HD201484', 'HD201496', 'HD201601A', 'HD201757', 'HD201852', 'HD202206', 'HD202389', 'HD202605', 'HD20278', 'HD202819', 'HD202871', 'HD203335', 'HD203384', 'HD203387', 'HD203413', 'HD203432', 'HD203473', 'HD203771', 'HD203897', 'HD203932', 'HD204287', 'HD204313', 'HD204385', 'HD20439', 'HD20492', 'HD204941', 'HD204943', 'HD204961', 'HD205289', 'HD205294', 'HD205471', 'HD205478', 'HD205536', 'HD205591', 'HD20559B', 'HD205739', 'HD206116', 'HD206172', 'HD206561', 'HD206630', 'HD206642', 'HD206658', 'HD206683', 'HD206837', 'HD206893', 'HD206998', 'HD207190', 'HD2071', 'HD207575', 'HD207699', 'HD207700', 'HD20781', 'HD20782', 'HD207832', 'HD207869', 'HD20794', 'HD208215', 'HD208272', 'HD208487', 'HD20852', 'HD208552', 'HD208573', 'HD20868', 'HD208704', 'HD208880', 'HD20894', 'HD208', 'HD209458', 'HD209625', 'HD209779', 'HD21011', 'HD210139', 'HD210272', 'HD210277', 'HD210329', 'HD210573', 'HD210591', 'HD210739', 'HD210752', 'HD210797', 'HD210975', 'HD211188', 'HD21120', 'HD211317', 'HD21132', 'HD211420', 'HD2114', 'HD211532', 'HD21161', 'HD211723', 'HD211976', 'HD212036', 'HD21209A', 'HD21209B', 'HD212231', 'HD212291', 'HD212301B', 'HD212571', 'HD212580', 'HD212658', 'HD212708', 'HD212728', 'HD212918', 'HD212943', 'HD213240', 'HD213366', 'HD2134', 'HD213575', 'HD213637', 'HD214094', 'HD21430', 'HD214385', 'HD214853', 'HD214867', 'HD214910', 'HD214953', 'HD215257', 'HD215497', 'HD215625', 'HD215641', 'HD215789', 'HD215864', 'HD215906', 'HD2160', 'HD216435', 'HD216437', 'HD216627', 'HD216770', 'HD21693', 'HD217065', 'HD217343', 'HD217395', 'HD217522', 'HD21759', 'HD217786', 'HD218249', 'HD218379', 'HD218396', 'HD218495', 'HD218504', 'HD218544', 'HD218566', 'HD218642', 'HD218750', 'HD21882', 'HD218860', 'HD218960', 'HD218994', 'HD219011', 'HD219057', 'HD219077', 'HD219249', 'HD21938', 'HD219449', 'HD219556', 'HD219688', 'HD21977', 'HD219781', 'HD219784', 'HD219828', 'HD220096', 'HD220308', 'HD220367', 'HD220456', 'HD22049', 'HD220507', 'HD22054', 'HD220572', 'HD220657', 'HD220689', 'HD220729', 'HD221146', 'HD221287', 'HD221503', 'HD221550', 'HD221575', 'HD221580', 'HD221627', 'HD221760', 'HD22177', 'HD221950', 'HD221954', 'HD222095', 'HD222259', 'HD22231', 'HD222422', 'HD222480', 'HD22249', 'HD222582', 'HD222595', 'HD222603', 'HD222669', 'HD222910', 'HD223121', 'HD223171', 'HD223238', 'HD223272', 'HD223282', 'HD223340', 'HD223352', 'HD223681', 'HD223781', 'HD223854', 'HD223889', 'HD224047', 'HD224063', 'HD22409', 'HD224230', 'HD224347', 'HD224383', 'HD224392', 'HD224393', 'HD224433', 'HD224538', 'HD224619', 'HD224685', 'HD224693', 'HD224789', 'HD224817', 'HD225197', 'HD225200', 'HD225297', 'HD22610', 'HD22663', 'HD22676', 'HD230017', 'HD23079', 'HD23127', 'HD231630', 'HD231701', 'HD23207', 'HD23319', 'HD23356', 'HD23456', 'HD23472', 'HD23719', 'HD23901', 'HD23940', 'HD24015', 'HD24062', 'HD240648', 'HD24085', 'HD24112', 'HD24160', 'HD24194', 'HD242780', 'HD24302', 'HD24331', 'HD244138', 'HD24465', 'HD24558', 'HD24633', 'HD24706', 'HD24712', 'HD24745', 'HD24916', 'HD25061', 'HD25102', 'HD25105', 'HD25120', 'HD25171', 'HD2529', 'HD25490', 'HD25535A', 'HD25565', 'HD25587', 'HD25673', 'HD2567', 'HD256', 'HD25825', 'HD25874', 'HD25912', 'HD26071', 'HD26297', 'HD2638', 'HD26430', 'HD26729', 'HD26864', 'HD26887', 'HD26911', 'HD26913', 'HD26923', 'HD269556', 'HD26956', 'HD26965', 'HD26967', 'HD2696', 'HD26990', 'HD27063', 'HD27371', 'HD273832', 'HD274255', 'HD27442A', 'HD27471', 'HD27483', 'HD27536', 'HD27631', 'HD2768', 'HD27697', 'HD27894', 'HD27935', 'HD27947', 'HD28028', 'HD28093', 'HD28185', 'HD28192', 'HD28254', 'HD28305', 'HD28307', 'HD2834', 'HD284113', 'HD285968', 'HD286123', 'HD28635', 'HD28701', 'HD28807', 'HD28821', 'HD28969', 'HD290038', 'HD290327', 'HD29137', 'HD29263', 'HD29291', 'HD29391', 'HD29399', 'HD29428', 'HD29488', 'HD29737', 'HD297396', 'HD29751', 'HD298541', 'HD298542', 'HD29875', 'HD298936', 'HD299043', 'HD29985', 'HD299896', 'HD30053', 'HD30306', 'HD304043', 'HD30447', 'HD304864', 'HD30669', 'HD30703', 'HD3074A', 'HD308353', 'HD308357', 'HD31103', 'HD312256', 'HD31527', 'HD31532', 'HD3167', 'HD318107', 'HD31822', 'HD320345', 'HD3221', 'HD3229', 'HD323684', 'HD325367', 'HD32564', 'HD32646', 'HD32724', 'HD327692', 'HD33142', 'HD33473', 'HD34642', 'HD347929', 'HD35650', 'HD358155', 'HD36051', 'HD36152', 'HD36379', 'HD37160', 'HD37226', 'HD37572', 'HD37990', 'HD38078', 'HD3821', 'HD38283', 'HD38355', 'HD38385', 'HD38459', 'HD38467', 'HD38677', 'HD39091', 'HD39427', 'HD3964', 'HD4021', 'HD40293', 'HD40483', 'HD40657', 'HD40808', 'HD41047', 'HD41071', 'HD41076', 'HD41087', 'HD41245', 'HD41248', 'HD41291', 'HD41323', 'HD41511', 'HD41695', 'HD41742', 'HD41842', 'HD4211', 'HD42269', 'HD42270', 'HD42538', 'HD42545', 'HD42556', 'HD42581', 'HD42606', 'HD42618', 'HD42659', 'HD42719', 'HD42729', 'HD42813', 'HD42936', 'HD4293', 'HD4308', 'HD4313', 'HD43180', 'HD43197', 'HD43717', 'HD43785', 'HD43834B', 'HD43940', 'HD43989', 'HD4398', 'HD44219', 'HD44420', 'HD44447', 'HD44573', 'HD4457', 'HD44594', 'HD44627', 'HD44665', 'HD44762', 'HD44804', 'HD45021', 'HD45081', 'HD45133', 'HD45184', 'HD45270', 'HD45289', 'HD45346', 'HD45364', 'HD45398', 'HD45665', 'HD45669', 'HD457', 'HD45977', 'HD46089', 'HD46116', 'HD46202', 'HD46375A', 'HD46415', 'HD46487', 'HD4649', 'HD46568', 'HD46769', 'HD47001', 'HD47186', 'HD4737', 'HD47536', 'HD47863', 'HD48115', 'HD48152', 'HD48265', 'HD48381', 'HD4838', 'HD483', 'HD48611', 'HD49035', 'HD49050', 'HD49091', 'HD49105', 'HD4915', 'HD49212', 'HD49334', 'HD49798', 'HD49855', 'HD49866', 'HD49877', 'HD49933A', 'HD49947', 'HD50060', 'HD50062', 'HD50169', 'HD50304', 'HD50310', 'HD50337', 'HD50445', 'HD50499', 'HD50590', 'HD50652', 'HD50778', 'HD50896', 'HD51754', 'HD51797', 'HD52265', 'HD52479', 'HD5248', 'HD52973', 'HD53116', 'HD5349', 'HD5388', 'HD54046', 'HD54100', 'HD54351', 'HD54521', 'HD5470', 'HD54810', 'HD55279', 'HD55524', 'HD55693', 'HD56096', 'HD56160', 'HD56259', 'HD56341', 'HD56351', 'HD56380', 'HD56413', 'HD564', 'HD56537', 'HD565', 'HD56618', 'HD56957', 'HD57061', 'HD57150', 'HD5722', 'HD57240', 'HD5737', 'HD57568', 'HD58489', 'HD58676', 'HD5891', 'HD58972', 'HD59088', 'HD59219', 'HD5938', 'HD59435', 'HD59711A', 'HD60060', 'HD60435', 'HD60574', 'HD60584', 'HD60629', 'HD61005', 'HD61051', 'HD6107', 'HD61273', 'HD61383', 'HD61447', 'HD61475', 'HD61772', 'HD61902', 'HD6192', 'HD61935', 'HD61986', 'HD6204', 'HD62128', 'HD62364', 'HD62412', 'HD6245', 'HD62713', 'HD62816', 'HD62847', 'HD62849', 'HD62943', 'HD63454', 'HD63487', 'HD6348', 'HD63563', 'HD63608', 'HD63685', 'HD63754', 'HD63765', 'HD63847', 'HD64152', 'HD64181', 'HD6434', 'HD64640', 'HD64942', 'HD65114', 'HD6512', 'HD6532', 'HD65562', 'HD6559', 'HD65695', 'HD6569', 'HD65949', 'HD65982', 'HD66039', 'HD66040', 'HD66168', 'HD66221', 'HD66255', 'HD66340', 'HD66428', 'HD66552', 'HD66653', 'HD66664', 'HD66740', 'HD66838', 'HD6718', 'HD67200', 'HD6735', 'HD67523', 'HD67578', 'HD6790', 'HD6793', 'HD67959', 'HD67990', 'HD67', 'HD68089', 'HD68161', 'HD68168', 'HD68284', 'HD68287', 'HD68402', 'HD68607', 'HD68694', 'HD68860', 'HD68978', 'HD69013', 'HD69247', 'HD69611', 'HD69655', 'HD69721', 'HD6980', 'HD70569', 'HD70642', 'HD7082', 'HD70889', 'HD70903', 'HD70982', 'HD71043', 'HD71155', 'HD7134', 'HD71479', 'HD71685', 'HD71835', 'HD7199', 'HD72359', 'HD72374', 'HD72579', 'HD72659', 'HD72660', 'HD7279', 'HD72892', 'HD72968', 'HD73121', 'HD73143', 'HD73256', 'HD73262', 'HD73267', 'HD73316', 'HD7345', 'HD73502', 'HD7355', 'HD73887', 'HD7399', 'HD74000', 'HD74006', 'HD74014', 'HD74455', 'HD7449', 'HD74521', 'HD74591', 'HD74698', 'HD74873', 'HD74874', 'HD74957', 'HD750', 'HD7515', 'HD75168', 'HD75171', 'HD75289', 'HD75302', 'HD75328', 'HD75413', 'HD75416', 'HD75445', 'HD75469', 'HD75530', 'HD75732', 'HD75745', 'HD75881', 'HD76188', 'HD76440', 'HD7661', 'HD76700', 'HD7672', 'HD76849', 'HD77058', 'HD77087', 'HD770', 'HD77110', 'HD77191', 'HD77278', 'HD77338', 'HD77462', 'HD77825', 'HD78004', 'HD78130', 'HD78286', 'HD78429', 'HD78534', 'HD78538', 'HD78558', 'HD78647', 'HD78660', 'HD78747', 'HD78964', 'HD7950', 'HD79601', 'HD80316', 'HD8038', 'HD80447', 'HD80883', 'HD80950', 'HD81009', 'HD81101', 'HD81169', 'HD81639', 'HD81700', 'HD81767', 'HD81797', 'HD82114', 'HD82165', 'HD82446', 'HD82516', 'HD82783', 'HD82895', 'HD8291', 'HD82943', 'HD83368', 'HD83373', 'HD83443', 'HD83446', 'HD83509', 'HD83529', 'HD83683', 'HD83953', 'HD84041', 'HD8406', 'HD84075', 'HD84273', 'HD84305', 'HD8446', 'HD84627', 'HD84687', 'HD84810', 'HD84937', 'HD85119', 'HD85249', 'HD8535', 'HD85390', 'HD85444', 'HD85512', 'HD8558', 'HD85725', 'HD85859', 'HD86006', 'HD86021', 'HD86065', 'HD86081', 'HD86133', 'HD86140', 'HD86171', 'HD86181', 'HD8651', 'HD86652', 'HD86997', 'HD8705', 'HD87109', 'HD87320', 'HD87479', 'HD87488', 'HD87521', 'HD87526', 'HD87566', 'HD8763', 'HD87810', 'HD87833', 'HD87838', 'HD88072', 'HD88084', 'HD8828', 'HD88474', 'HD8859', 'HD88656', 'HD88885', 'HD88955', 'HD8912', 'HD89147', 'HD89328', 'HD89454', 'HD89668', 'HD89749', 'HD89839', 'HD8985', 'HD89920', 'HD89965', 'HD90028', 'HD90081', 'HD90132', 'HD90133', 'HD90156', 'HD90177', 'HD90422', 'HD90520', 'HD90569', 'HD9065', 'HD90722', 'HD90774', 'HD90812', 'HD9081', 'HD90905', 'HD90926', 'HD90994', 'HD91121', 'HD91267', 'HD91316', 'HD91345', 'HD91375', 'HD91379', 'HD91682', 'HD9174', 'HD91816', 'HD91988', 'HD92003', 'HD92106', 'HD92139', 'HD92245', 'HD9246', 'HD92536', 'HD92547', 'HD92588', 'HD92719', 'HD92788', 'HD9289', 'HD93153', 'HD93351', 'HD93372', 'HD93385', 'HD93396', 'HD93489', 'HD93507', 'HD9362', 'HD93745', 'HD93849', 'HD93932', 'HD94126', 'HD94151', 'HD94270', 'HD94771', 'HD94964', 'HD95086', 'HD95136', 'HD95456', 'HD95521', 'HD95533', 'HD95542', 'HD9578', 'HD95799', 'HD95860', 'HD95879', 'HD96064', 'HD9608', 'HD96116', 'HD96118', 'HD96175', 'HD96273', 'HD96276', 'HD96290', 'HD96423', 'HD96445', 'HD96494', 'HD96544', 'HD96568', 'HD965', 'HD96673', 'HD96700', 'HD96723', 'HD9672', 'HD96789', 'HD967', 'HD96819', 'HD97037', 'HD97244', 'HD97320', 'HD9742', 'HD97783', 'HD9782', 'HD9796', 'HD97998', 'HD98248', 'HD98280', 'HD98284', 'HD98356', 'HD98430', 'HD98457', 'HD984', 'HD98586', 'HD987', 'HD9920', 'HD99211', 'HD99492', 'HD9986', 'HIP100233', 'HIP100396', 'HIP101846', 'HIP102025', 'HIP102301', 'HIP102495', 'HIP102964', 'HIP103019', 'HIP103039', 'HIP103150', 'HIP10395', 'HIP103972', 'HIP104083', 'HIP1044', 'HIP104644', 'HIP104856', 'HIP105368', 'HIP105388', 'HIP105506', 'HIP10689', 'HIP106931', 'HIP107133', 'HIP107345', 'HIP10741', 'HIP107772', 'HIP108020', 'HIP108065', 'HIP108095', 'HIP108340', 'HIP108930', 'HIP108950', 'HIP109012', 'HIP109149', 'HIP109638', 'HIP110156', 'HIP110443', 'HIP110692', 'HIP110714', 'HIP111162', 'HIP111978', 'HIP1119', 'HIP112414', 'HIP113596', 'HIP113850', 'HIP114044', 'HIP115126', 'HIP115332', 'HIP115381', 'HIP115803', 'HIP115861', 'HIP116003', 'HIP116317', 'HIP116350', 'HIP116374', 'HIP116492', 'HIP117219', 'HIP11739', 'HIP11989', 'HIP12065', 'HIP12147', 'HIP12961', 'HIP14475', 'HIP15587', 'HIP1608', 'HIP16445', 'HIP16536', 'HIP17044', 'HIP17157', 'HIP17316', 'HIP17365', 'HIP17695', 'HIP17956', 'HIP18180', 'HIP18918', 'HIP1993', 'HIP20232', 'HIP20267', 'HIP20444', 'HIP20699', 'HIP21539', 'HIP21991', 'HIP22059', 'HIP22424', 'HIP2247', 'HIP23309', 'HIP23344', 'HIP25578', 'HIP25612', 'HIP25654', 'HIP26013', 'HIP26542', 'HIP27397', 'HIP27928', 'HIP29322', 'HIP29788', 'HIP31293', 'HIP31639', 'HIP31862', 'HIP31878', 'HIP32127', 'HIP32393', 'HIP3249', 'HIP32812', 'HIP34285', 'HIP35305', 'HIP35366', 'HIP35937', 'HIP35965', 'HIP36338', 'HIP36985', 'HIP37217', 'HIP38594', 'HIP39470', 'HIP39611', 'HIP40038', 'HIP40459', 'HIP40774', 'HIP41659', 'HIP41662', 'HIP41795', 'HIP4328', 'HIP43310', 'HIP43313', 'HIP4394', 'HIP43973', 'HIP44376', 'HIP4468', 'HIP44899', 'HIP45301', 'HIP46634', 'HIP46655', 'HIP46933', 'HIP47', 'HIP48762', 'HIP4877', 'HIP48953', 'HIP49067', 'HIP51073', 'HIP5114', 'HIP51344', 'HIP51443', 'HIP51588', 'HIP5215', 'HIP523', 'HIP52776', 'HIP52997', 'HIP53559', 'HIP53639', 'HIP53911', 'HIP54365', 'HIP54373', 'HIP54446', 'HIP54469', 'HIP54597', 'HIP55787', 'HIP55789', 'HIP55955', 'HIP56489', 'HIP56838', 'HIP57459', 'HIP58285', 'HIP58289', 'HIP58348', 'HIP59109', 'HIP59341', 'HIP59925', 'HIP60051', 'HIP60718', 'HIP60910', 'HIP61406', 'HIP61505', 'HIP61872', 'HIP62452', 'HIP6276', 'HIP63157', 'HIP63719', 'HIP64428', 'HIP65924', 'HIP65956', 'HIP6639', 'HIP66678', 'HIP66918', 'HIP67126', 'HIP67164', 'HIP68900', 'HIP69224', 'HIP69570', 'HIP69929', 'HIP70317', 'HIP70475', 'HIP7058', 'HIP70681', 'HIP70975', 'HIP71253', 'HIP7228', 'HIP72378', 'HIP72779', 'HIP73194', 'HIP74190', 'HIP7588', 'HIP7646', 'HIP76550', 'HIP76668', 'HIP78242', 'HIP7829', 'HIP78345', 'HIP78395', 'HIP78408', 'HIP79578', 'HIP80018', 'HIP80083', 'HIP80683', 'HIP80817', 'HIP81018', 'HIP81084', 'HIP8119', 'HIP81465', 'HIP82283', 'HIP82323', 'HIP82357', 'HIP8361', 'HIP8438', 'HIP84405', 'HIP85360', 'HIP86214', 'HIP86509', 'HIP8691', 'HIP87607', 'HIP88064', 'HIP88316', 'HIP88371', 'HIP89892', 'HIP90194', 'HIP90246', 'HIP92444', 'HIP92451', 'HIP92592', 'HIP93287', 'HIP93745', 'HIP9398', 'HIP948', 'HIP95903', 'HIP9603', 'HIP96101', 'HIP96240', 'HIP9786', 'HIP98105', 'HIP98106', 'HIP98470', 'HIP99273', 'HIP99385', 'HIP99695', 'K2-132', 'K2-138', 'K2-140', 'K2-27', 'K2-32', 'TYC1405-1421-1', 'TYC1472-1436-2', 'TYC15-882-1', 'TYC1800-2170-1', 'TYC2-1155-1', 'TYC211-706-1', 'TYC22-1354-1', 'TYC281-875-1', 'TYC30-116-1', 'TYC339-329-1', 'TYC459-1288-1', 'TYC4682-800-1', 'TYC4851-299-1', 'TYC4858-1263-1', 'TYC4873-651-1', 'TYC4874-1474-1', 'TYC4890-625-1', 'TYC4916-897-1', 'TYC4927-1063-1', 'TYC4940-797-1', 'TYC4943-1016-1', 'TYC4947-834-1', 'TYC4948-327-1', 'TYC5126-2775-1', 'TYC5204-993-1', 'TYC5290-462-1', 'TYC534-2183-1', 'TYC5434-3430-1', 'TYC5839-876-1', 'TYC5889-271-1', 'TYC5965-749-1', 'TYC6051-445-1', 'TYC6121-11-1', 'TYC6170-95-1', 'TYC6193-663-1', 'TYC6224-870-1', 'TYC6244-604-1', 'TYC6283-1022-1', 'TYC6296-1176-1', 'TYC6296-1341-1', 'TYC6296-1556-1', 'TYC6296-181-1', 'TYC6296-1857-1', 'TYC6296-303-1', 'TYC6296-478-1', 'TYC6296-483-1', 'TYC6296-654-1', 'TYC6300-257-1', 'TYC6300-98-1', 'TYC6360-482-1', 'TYC6430-1245-1', 'TYC6484-222-1', 'TYC6604-547-1', 'TYC6609-79-1', 'TYC6609-880-1', 'TYC6622-794-1', 'TYC6627-950-1', 'TYC6646-1539-1', 'TYC6716-330-1', 'TYC6764-1356-1', 'TYC6822-2535-1', 'TYC6871-1581-1', 'TYC6887-2075-1', 'TYC6901-2187-1', 'TYC6914-1673-1', 'TYC7047-1324-1', 'TYC7087-2536-1', 'TYC7104-505-1', 'TYC7193-571-1', 'TYC72-1041-1', 'TYC7282-1298-1', 'TYC7442-558-1', 'TYC7501-987-1', 'TYC7528-890-1', 'TYC7588-318-1', 'TYC7646-3253-1', 'TYC7752-648-1', 'TYC7971-310-1', 'TYC8015-1020-1', 'TYC8015-1177-1', 'TYC8096-1049-1', 'TYC813-243-1', 'TYC8224-859-1', 'TYC8229-2228-1', 'TYC8233-68-1', 'TYC8252-1071-1', 'TYC8378-64-1', 'TYC8444-1168-1', 'TYC8527-329-1', 'TYC8545-158-1', 'TYC8633-394-1', 'TYC8647-2057-1', 'TYC8688-1180-1', 'TYC8770-239-1', 'TYC8778-850-1', 'TYC8870-372-1', 'TYC8963-1429-1', 'TYC8963-1543-1', 'TYC9034-968-1', 'TYC9055-1987-1', 'TYC9111-1423-1', 'TYC9114-1267-1', 'TYC9379-1149-1', 'WASP-101', 'WASP-105', 'WASP-107', 'WASP-121', 'WASP-127', 'WASP-129', 'WASP-130', 'WASP-132', 'WASP-157', 'WASP-15', 'WASP-16', 'WASP-17', 'WASP-19', 'WASP-20', 'WASP-22', 'WASP-25', 'WASP-2A', 'WASP-30', 'WASP-31', 'WASP-34', 'WASP-36', 'WASP-40', 'WASP-41', 'WASP-42', 'WASP-43', 'WASP-4', 'WASP-52', 'WASP-54', 'WASP-55', 'WASP-61', 'WASP-62', 'WASP-6', 'WASP-76', 'WASP-77A', 'WASP-80', 'WASP-89', 'WASP-8', 'WASP-99']


data_files_HIRES = ["RVs HIRES + NZP correction",
              "RVs HIRES (original)",
              "s-index","h-index"]

data_files_ind_HIRES = [1,3,5,6]

targets_HIRES = ['0429-01631-1', '0438-00260-1', '0748-01711-1', 'BD+423146', 'BD-103166', 'G097-054', 'G161-29', 'G192-13', 'G195-59', 'G205-028', 'G207-019', 'G244-047', 'G60-06', 'GJ3470', 'GL105B', 'GL107B', 'GL109', 'GL1245B', 'GL2066', 'GL226', 'GL239', 'GL250B', 'GL26', 'GL273', 'GL317', 'GL357', 'GL382', 'GL388', 'GL393', 'GL397', 'GL4063', 'GL406', 'GL408', 'GL412A', 'GL433', 'GL445', 'GL450', 'GL47', 'GL486', 'GL48', 'GL49', 'GL514', 'GL528B', 'GL569A', 'GL625', 'GL667C', 'GL686', 'GL687', 'GL694', 'GL699', 'GL745A', 'GL745B', 'GL793', 'GL803', 'GL806', 'GL83.1', 'GL876', 'GL87', 'GL905', 'GL908', 'GMAUR', 'HD10002', 'HD10008', 'HD10013', 'HD10015', 'HD100180', 'HD100337', 'HD100623', 'HD101206', 'HD101259', 'HD101348', 'HD10145', 'HD101472', 'HD101501', 'HD101675', 'HD101847', 'HD101904', 'HD101959', 'HD102071', 'HD102158', 'HD102195', 'HD102361', 'HD102365', 'HD102444', 'HD102956', 'HD103047', 'HD103095', 'HD103417', 'HD103432', 'HD103459', 'HD103616', 'HD103813', 'HD103828', 'HD103829', 'HD103847', 'HD103890', 'HD103932', 'HD104017', 'HD104067', 'HD104263', 'HD104304', 'HD10436', 'HD104389', 'HD10442', 'HD104437', 'HD104556', 'HD104588', 'HD10476', 'HD104800', 'HD104860', 'HD105113', 'HD105279', 'HD105304', 'HD105618', 'HD105631', 'HD105811', 'HD105', 'HD106088', 'HD106116', 'HD106156', 'HD106279', 'HD106314', 'HD106421', 'HD106949', 'HD10697', 'HD10700', 'HD107087', 'HD107146', 'HD107148', 'HD107211', 'HD10780', 'HD10790', 'HD108189', 'HD108351', 'HD108863', 'HD108874', 'HD108942', 'HD109011', 'HD109159', 'HD109218', 'HD109331', 'HD109358', 'HD109409', 'HD109542', 'HD109718', 'HD109749', 'HD109929', 'HD10995', 'HD110044', 'HD11020', 'HD110315', 'HD110463', 'HD110537', 'HD110743', 'HD110897', 'HD111031', 'HD111153', 'HD11131', 'HD111395', 'HD111431', 'HD111484A', 'HD111484B', 'HD111515', 'HD111563', 'HD111631', 'HD111814', 'HD112019', 'HD112115', 'HD112208', 'HD112257', 'HD112337', 'HD112415', 'HD11271', 'HD112914', 'HD112973', 'HD112988', 'HD113194', 'HD113414', 'HD113490', 'HD113595', 'HD11373', 'HD113983', 'HD114161', 'HD114174', 'HD114375', 'HD11437', 'HD114613', 'HD114659', 'HD114729', 'HD114783', 'HD114946', 'HD11506', 'HD115404A', 'HD115589', 'HD115617', 'HD116029', 'HD116258', 'HD116321', 'HD116442', 'HD116443', 'HD116956', 'HD117122', 'HD117176', 'HD117207', 'HD117434', 'HD117497', 'HD117623', 'HD117936', 'HD11850', 'HD118670', 'HD118744', 'HD118914', 'HD119058', 'HD11964A', 'HD119802', 'HD119850', 'HD120066', 'HD120219', 'HD12039', 'HD120467', 'HD120476A', 'HD12051', 'HD120528', 'HD1205', 'HD120636', 'HD121320', 'HD121550', 'HD121579', 'HD122064', 'HD122120', 'HD122253', 'HD122255', 'HD122303', 'HD12235', 'HD122517', 'HD122652', 'HD122973', 'HD123239', 'HD123265', 'HD1234-00069-1', 'HD12380', 'HD123812', 'HD124102', 'HD124106', 'HD124257A', 'HD124257B', 'HD124292', 'HD124641', 'HD124642', 'HD12484', 'HD125184', 'HD125455', 'HD125607', 'HD125612', 'HD126053', 'HD126203', 'HD126583', 'HD126614', 'HD12661', 'HD126631', 'HD126831', 'HD126990', 'HD126991', 'HD127334', 'HD127374', 'HD127506', 'HD127741', 'HD128095', 'HD128165', 'HD128311', 'HD128428', 'HD12846', 'HD128642', 'HD129010', 'HD12911', 'HD129191', 'HD129333', 'HD1293', 'HD129471', 'HD12974', 'HD129814', 'HD130004', 'HD130048', 'HD130087', 'HD130307', 'HD130322', 'HD13043', 'HD130666', 'HD130672', 'HD130871', 'HD130992', 'HD131117', 'HD131183', 'HD131496', 'HD131509', 'HD132133', 'HD132142', 'HD132173', 'HD132307', 'HD132505', 'HD1326B', 'HD1326', 'HD133125', 'HD133233', 'HD133295', 'HD13345', 'HD13357', 'HD13361', 'HD13382', 'HD134319', 'HD134353', 'HD134439', 'HD134440', 'HD13483', 'HD134987', 'HD13507', 'HD135101A', 'HD1352-01149-1', 'HD135446', 'HD135599', 'HD13579', 'HD13584', 'HD135872', 'HD13612B', 'HD136159', 'HD136274', 'HD136352', 'HD136418', 'HD136442', 'HD136513', 'HD136618', 'HD136713', 'HD136834', 'HD136925', 'HD137368', 'HD13747', 'HD137631', 'HD13773', 'HD137778', 'HD138004', 'HD13836', 'HD138549', 'HD138600', 'HD13871', 'HD138776', 'HD1388', 'HD13931', 'HD139323', 'HD139457', 'HD139477', 'HD139813', 'HD139879', 'HD139907', 'HD13997', 'HD140025', 'HD140538A', 'HD140721', 'HD140913', 'HD141004', 'HD141085', 'HD141272', 'HD141399', 'HD141937', 'HD142091', 'HD14214', 'HD142229', 'HD14223', 'HD142245', 'HD142267', 'HD142626', 'HD142943', 'HD143006', 'HD143174', 'HD143291', 'HD143332', 'HD143761', 'HD14412', 'HD144287', 'HD144363', 'HD144579', 'HD144585', 'HD144872', 'HD144988', 'HD1449', 'HD14519', 'HD145229', 'HD145331', 'HD145428', 'HD145675', 'HD145809', 'HD145934', 'HD145958A', 'HD145958B', 'HD146050', 'HD1461', 'HD146233', 'HD146278', 'HD146362B', 'HD14651', 'HD14655', 'HD146775', 'HD147062', 'HD147231', 'HD147379A', 'HD147379B', 'HD147506', 'HD147750', 'HD147752', 'HD147776', 'HD147887', 'HD148164', 'HD148284', 'HD148319', 'HD148428', 'HD148467', 'HD14855', 'HD149026', 'HD149143', 'HD149661', 'HD149724', 'HD149750', 'HD149760', 'HD1497', 'HD149806', 'HD150122', 'HD1502', 'HD150331', 'HD150433', 'HD150437', 'HD150554', 'HD150698', 'HD150706', 'HD150936', 'HD151288', 'HD151329', 'HD151504', 'HD151522', 'HD151541', 'HD151852', 'HD151877', 'HD151995', 'HD152125', 'HD152391', 'HD152555', 'HD152581', 'HD152733', 'HD152792', 'HD152878', 'HD15335', 'HD153378', 'HD153458', 'HD153525', 'HD153557', 'HD15367', 'HD154088', 'HD154144', 'HD154325', 'HD154345', 'HD154363', 'HD154656', 'HD154697', 'HD154994', 'HD155358', 'HD155415', 'HD155456', 'HD155712', 'HD155817', 'HD155968', 'HD156026', 'HD156079', 'HD156146', 'HD156279', 'HD156342', 'HD156365', 'HD156549', 'HD156668', 'HD156826', 'HD156846', 'HD156985', 'HD157172', 'HD157214', 'HD157299', 'HD157338', 'HD157347', 'HD157881', 'HD158038', 'HD158173', 'HD158449', 'HD158614', 'HD158633', 'HD159062', 'HD159063', 'HD159222', 'HD159868', 'HD1605', 'HD160693', 'HD161198', 'HD161284', 'HD16141', 'HD161424', 'HD161479', 'HD16160', 'HD16175', 'HD161797', 'HD161848', 'HD161897', 'HD162231', 'HD162232', 'HD16275', 'HD16287', 'HD16297', 'HD163153', 'HD163489', 'HD163589', 'HD163607', 'HD16397', 'HD16417', 'HD164330', 'HD164507', 'HD164509', 'HD164595', 'HD164651', 'HD164922', 'HD165109', 'HD165173', 'HD165222', 'HD165341A', 'HD165401', 'HD16559', 'HD165672', 'HD16623', 'HD166620', 'HD1666', 'HD166', 'HD16702', 'HD167215', 'HD167216', 'HD167389', 'HD16760', 'HD167665', 'HD168009', 'HD168443', 'HD168603', 'HD168746', 'HD169830', 'HD169889', 'HD170174', 'HD17037', 'HD170469', 'HD170493', 'HD170512', 'HD170657', 'HD171067', 'HD171238', 'HD17152', 'HD17156', 'HD171665', 'HD17190', 'HD171918', 'HD171999', 'HD172051', 'HD17230', 'HD172310', 'HD172513', 'HD17354', 'HD173701', 'HD173739', 'HD173740', 'HD173818', 'HD17382', 'HD174080', 'HD175425', 'HD175541', 'HD175742', 'HD176377', 'HD176414', 'HD17660', 'HD176733', 'HD176845', 'HD176982', 'HD177033', 'HD1770', 'HD177153', 'HD177274', 'HD177830', 'HD178251', 'HD178911B', 'HD179079', 'HD179306', 'HD179596', 'HD179949', 'HD179957', 'HD179958', 'HD180053', 'HD180161', 'HD180617', 'HD180684', 'HD181234', 'HD181253', 'HD18131', 'HD18143', 'HD182189', 'HD182407', 'HD182488', 'HD182572', 'HD182619', 'HD183162', 'HD183263', 'HD183298', 'HD1832', 'HD183473', 'HD183650', 'HD183658', 'HD183870', 'HD18436A', 'HD18436B', 'HD18445', 'HD185144', 'HD185269', 'HD185295', 'HD185414', 'HD185501', 'HD186104', 'HD18632', 'HD186408', 'HD186427', 'HD186932', 'HD18702', 'HD187123', 'HD187237', 'HD18747', 'HD187897', 'HD187923', 'HD187944', 'HD188015', 'HD18803', 'HD188268', 'HD188298', 'HD188311', 'HD188345', 'HD188376', 'HD188510', 'HD188512', 'HD18907', 'HD18916', 'HD189625', 'HD189627', 'HD189733', 'HD18975', 'HD18993', 'HD190007', 'HD190067', 'HD19019', 'HD190228', 'HD19034', 'HD190360', 'HD190404', 'HD190406', 'HD190412', 'HD190571', 'HD190594', 'HD190821', 'HD190931', 'HD191067', 'HD191089', 'HD191408', 'HD191785', 'HD192020', 'HD192263', 'HD192310', 'HD192343', 'HD192344', 'HD192699', 'HD19308', 'HD193202', 'HD193342', 'HD193690', 'HD19373', 'HD193795', 'HD193901', 'HD194080', 'HD194110', 'HD194541', 'HD19467', 'HD195019B', 'HD195019', 'HD19502', 'HD195034', 'HD19522', 'HD1955-00658-1', 'HD195564', 'HD195787', 'HD195987', 'HD196124', 'HD19617', 'HD19618', 'HD196201', 'HD19632', 'HD19638', 'HD19659', 'HD196676', 'HD19668', 'HD196761', 'HD196850', 'HD197076', 'HD197162', 'HD197623', 'HD19773', 'HD198387', 'HD198425', 'HD198483', 'HD198550', 'HD198599', 'HD198683', 'HD198802', 'HD199019', 'HD199100', 'HD199178', 'HD199255', 'HD199305', 'HD199381', 'HD199476', 'HD199580', 'HD199683', 'HD199960', 'HD200078', 'HD200156', 'HD200538', 'HD200565', 'HD200964', 'HD200968', 'HD201091', 'HD201092', 'HD201203', 'HD201219', 'HD201378', 'HD20155', 'HD201651', 'HD20165', 'HD201924', 'HD201989', 'HD202560', 'HD202575', 'HD2025', 'HD202696', 'HD202751', 'HD202917', 'HD203030', 'HD203471', 'HD203473', 'HD204277', 'HD20439', 'HD204587', 'HD204814', 'HD205163', 'HD205351', 'HD205353', 'HD205739', 'HD205855', 'HD205905', 'HD206116', 'HD20618', 'HD20619', 'HD206374', 'HD206387', 'HD206610', 'HD206658', 'HD20675', 'HD207077', 'HD207485', 'HD207583', 'HD207832', 'HD207839', 'HD207874', 'HD207897', 'HD207992', 'HD207994', 'HD208038', 'HD208202', 'HD208313', 'HD2085', 'HD208801', 'HD208880', 'HD208897', 'HD209203', 'HD209253', 'HD209290', 'HD209340', 'HD209393', 'HD209458', 'HD209599', 'HD209706', 'HD209875', 'HD210011', 'HD210144', 'HD21019A', 'HD210277', 'HD210302', 'HD210312', 'HD210320', 'HD210323', 'HD210373', 'HD210391', 'HD210392', 'HD210702', 'HD211038', 'HD211080', 'HD211681', 'HD211810', 'HD21197', 'HD212291', 'HD212315', 'HD212585', 'HD212733', 'HD212801', 'HD212924', 'HD213042', 'HD213066', 'HD21313', 'HD213329', 'HD21340', 'HD213472', 'HD213519', 'HD213628', 'HD21449', 'HD214683', 'HD214749', 'HD214823', 'HD215032', 'HD215152', 'HD215274', 'HD215500', 'HD215578', 'HD215625', 'HD215704', 'HD216175', 'HD216191', 'HD216259', 'HD216275', 'HD216520', 'HD216722', 'HD216772', 'HD216783', 'HD216803', 'HD216899', 'HD217004', 'HD217014', 'HD217107', 'HD217165', 'HD217357', 'HD21742', 'HD217523', 'HD217591', 'HD21774', 'HD217850', 'HD217877', 'HD217987', 'HD218133', 'HD218168', 'HD218209', 'HD218445', 'HD21847', 'HD218566', 'HD218868', 'HD218935', 'HD219134', 'HD219428', 'HD219538', 'HD219542', 'HD219623', 'HD219770', 'HD219781', 'HD219828', 'HD219834B', 'HD220221', 'HD220339', 'HD22049', 'HD220554', 'HD22072', 'HD220952', 'HD221149', 'HD221354', 'HD221356', 'HD221561', 'HD221822', 'HD221851', 'HD221974', 'HD222038', 'HD222089', 'HD222582', 'HD222697', 'HD22282', 'HD222986', 'HD223238', 'HD223315', 'HD223498', 'HD223691', 'HD223869', 'HD224032', 'HD224040', 'HD224383', 'HD224601', 'HD224619', 'HD224693', 'HD22484', 'HD224983', 'HD225118', 'HD225213', 'HD225261', 'HD22879', 'HD230409', 'HD230999', 'HD231157', 'HD231701', 'HD23249', 'HD232979', 'HD233153', 'HD2331', 'HD23356', 'HD233641', 'HD234314', 'HD23439', 'HD237903', 'HD238069', 'HD238433', 'HD23952', 'HD239960', 'HD24040', 'HD24213', 'HD24238', 'HD24341', 'HD24365', 'HD24451', 'HD24496', 'HD24505', 'HD245409', 'HD24727', 'HD24892', 'HD24916', 'HD249', 'HD25069', 'HD25311', 'HD25329', 'HD25457', 'HD2564', 'HD25665', 'HD25682', 'HD25825', 'HD2589', 'HD25998', 'HD26151', 'HD26161', 'HD26257', 'HD265866', 'HD26633', 'HD26736', 'HD26794', 'HD26965', 'HD26990', 'HD27282', 'HD27496', 'HD27530', 'HD27732', 'HD28005', 'HD28099', 'HD281540', 'HD2816', 'HD28185', 'HD28187', 'HD281934', 'HD28237', 'HD28343', 'HD28388', 'HD283', 'HD285968', 'HD28946', 'HD29419', 'HD29461', 'HD2946', 'HD29528', 'HD29818', 'HD29883', 'HD2992', 'HD30339', 'HD30495', 'HD30649', 'HD30708', 'HD30712', 'HD3074A', 'HD30882', 'HD31018', 'HD31253', 'HD31392', 'HD31412', 'HD3141', 'HD31423', 'HD31451', 'HD31560', 'HD31664', 'HD31675', 'HD31693', 'HD31966', 'HD32147', 'HD32483', 'HD32850', 'HD32923', 'HD32963', 'HD33021', 'HD33142', 'HD33283', 'HD33298', 'HD33334', 'HD33636', 'HD33793', 'HD33822', 'HD3388-01009-1', 'HD3404', 'HD3411-02491-1', 'HD34411', 'HD34445', 'HD34575', 'HD34721', 'HD34745', 'HD34887', 'HD34957', 'HD3545', 'HD355183', 'HD35627', 'HD35850', 'HD35974', 'HD36003', 'HD36130', 'HD36308', 'HD36395', 'HD3651', 'HD3684', 'HD36974', 'HD37006', 'HD37008', 'HD37124', 'HD37213', 'HD37216', 'HD37250', 'HD37394', 'HD37445', 'HD37484', 'HD37605', 'HD3765', 'HD377', 'HD3795', 'HD37962', 'HD37986', 'HD38207', 'HD38230', 'HD38392', 'HD38393', 'HD38467', 'HD38505', 'HD38529', 'HD38801', 'HD38858', 'HD38949', 'HD38A', 'HD38B', 'HD39094', 'HD39715', 'HD39828', 'HD39881', 'HD3995-01436-1', 'HD40397', 'HD40537', 'HD40647', 'HD4075', 'HD40979', 'HD4113', 'HD41484', 'HD41593', 'HD41700', 'HD4203', 'HD4208', 'HD42250', 'HD4256', 'HD42581', 'HD42618', 'HD4307', 'HD4313', 'HD43162', 'HD4356-01014-1', 'HD43691', 'HD43745', 'HD43947', 'HD44420', 'HD44663', 'HD44985', 'HD45067', 'HD45184', 'HD45210', 'HD45350', 'HD45588', 'HD45652', 'HD457', 'HD46090', 'HD46122', 'HD4614B', 'HD4614', 'HD4628', 'HD4635', 'HD46375', 'HD47157', 'HD47186', 'HD4741', 'HD4747', 'HD47562', 'HD47752', 'HD48122', 'HD48345', 'HD4835-00774-1', 'HD48682', 'HD48938', 'HD4915', 'HD49197', 'HD49674', 'HD5015', 'HD50275', 'HD50281', 'HD50499', 'HD50554', 'HD50639', 'HD50692', 'HD50806', 'HD51046', 'HD51219', 'HD51272', 'HD5133', 'HD51419', 'HD51866', 'HD52265', 'HD52456', 'HD52711', 'HD52919', 'HD5319', 'HD531A', 'HD531B', 'HD533', 'HD53532', 'HD53665', 'HD5372', 'HD5470', 'HD55575', 'HD55647', 'HD55696', 'HD56083', 'HD5608', 'HD56122', 'HD56274', 'HD56303', 'HD56957', 'HD57204', 'HD57813', 'HD58727', 'HD5873', 'HD58781', 'HD5946', 'HD60234', 'HD60491', 'HD60737', 'HD61005', 'HD6101', 'HD61606', 'HD61994', 'HD61995', 'HD62613', 'HD63754', 'HD64090', 'HD64324', 'HD64413', 'HD64730', 'HD6512', 'HD65277', 'HD65430', 'HD65486', 'HD65583', 'HD6558', 'HD66171', 'HD66428', 'HD66751', 'HD6697', 'HD6715', 'HD6734', 'HD67458', 'HD67767', 'HD68017', 'HD68165', 'HD68168', 'HD6872A', 'HD6872B', 'HD68978', 'HD68988', 'HD69076', 'HD691', 'HD6963', 'HD69809', 'HD69830', 'HD69960', 'HD70516', 'HD71334', 'HD71479', 'HD71881', 'HD72003', 'HD72429', 'HD72440', 'HD72490', 'HD72659', 'HD72673', 'HD72905', 'HD73226', 'HD73256', 'HD73534', 'HD73667', 'HD74156', 'HD74390', 'HD745', 'HD74669', 'HD7510', 'HD75393', 'HD75407', 'HD75732B', 'HD75732', 'HD75784', 'HD75898', 'HD76218', 'HD76445', 'HD76909', 'HD77172', 'HD7727', 'HD77818', 'HD78255', 'HD79210', 'HD79211', 'HD7924', 'HD79498', 'HD79555', 'HD80355', 'HD80367', 'HD8038', 'HD804', 'HD80606', 'HD80811', 'HD8250', 'HD82886', 'HD82943', 'HD83024', 'HD8328', 'HD83394', 'HD83443', 'HD834', 'HD8375', 'HD8389', 'HD84035', 'HD84117', 'HD84453', 'HD8446', 'HD8467', 'HD84737', 'HD85301', 'HD85440', 'HD85472', 'HD85512', 'HD8553', 'HD85689', 'HD85725', 'HD8574', 'HD8594', 'HD86081', 'HD8648', 'HD86728', 'HD87001', 'HD87230', 'HD87359', 'HD87424', 'HD8765', 'HD87669', 'HD87836', 'HD87883', 'HD88072', 'HD88133', 'HD88218', 'HD88230', 'HD8828', 'HD88371', 'HD88654', 'HD88656', 'HD88725', 'HD88986', 'HD8907', 'HD8912', 'HD89269', 'HD89391', 'HD8939', 'HD8956', 'HD89793', 'HD90043', 'HD90054', 'HD90125', 'HD90156', 'HD9070', 'HD90711', 'HD90722', 'HD90792', 'HD9081', 'HD90875', 'HD90905', 'HD9113', 'HD91204', 'HD91876', 'HD9218', 'HD92222A', 'HD92222B', 'HD92266', 'HD92719', 'HD92788', 'HD9280', 'HD92855', 'HD92945', 'HD9331', 'HD93396', 'HD93461', 'HD93745', 'HD93864', 'HD9407', 'HD94151', 'HD94178', 'HD9446', 'HD9472', 'HD94834', 'HD95088', 'HD95089', 'HD95128', 'HD95188', 'HD9518A', 'HD9540A', 'HD95526', 'HD95622', 'HD9562', 'HD95650', 'HD95735', 'HD95900', 'HD96167', 'HD96361', 'HD96529', 'HD96612', 'HD96683', 'HD96700', 'HD96937', 'HD97101B', 'HD97101', 'HD97343', 'HD97601', 'HD97658', 'HD97854', 'HD98219', 'HD98281', 'HD98553', 'HD98618', 'HD98630', 'HD98736', 'HD98744', 'HD99109', 'HD99491', 'HD99492', 'HD99706', 'HD9986', 'HD99934', 'HIP10072', 'HIP101262', 'HIP102870', 'HIP103039', 'HIP103256', 'HIP103269', 'HIP10337', 'HIP103650', 'HIP104092', 'HIP10416', 'HIP104432', 'HIP10449', 'HIP105341', 'HIP1055', 'HIP105904', 'HIP106924', 'HIP1078', 'HIP108056', 'HIP108940', 'HIP109388', 'HIP109555', 'HIP109980', 'HIP11000', 'HIP110245', 'HIP11048', 'HIP110655', 'HIP110750', 'HIP110774', 'HIP111854', 'HIP112460', 'HIP112496', 'HIP112918', 'HIP113026', 'HIP113207', 'HIP113409', 'HIP114156', 'HIP114411', 'HIP114914', 'HIP115004', 'HIP115332', 'HIP115562', 'HIP116215', 'HIP116838', 'HIP117197', 'HIP117492', 'HIP117559', 'HIP117886', 'HIP117946', 'HIP117972', 'HIP118261', 'HIP118310', 'HIP12493', 'HIP12709', 'HIP12929', 'HIP1294', 'HIP13258', 'HIP13342', 'HIP13375', 'HIP13398', 'HIP13460', 'HIP1368', 'HIP1386', 'HIP14113', 'HIP14729', 'HIP14810', 'HIP15095', 'HIP15312', 'HIP1532', 'HIP15366', 'HIP15563', 'HIP15673', 'HIP15904', 'HIP16134', 'HIP16404', 'HIP17346', 'HIP1734', 'HIP17496', 'HIP18774', 'HIP19165', 'HIP19981', 'HIP20218', 'HIP20359', 'HIP21276', 'HIP21556', 'HIP22288', 'HIP2247', 'HIP22627', 'HIP22762', 'HIP22907', 'HIP23512', 'HIP23516', 'HIP24121', 'HIP24284', 'HIP25220', 'HIP26196', 'HIP26857', 'HIP29052', 'HIP29067', 'HIP29548', 'HIP30112', 'HIP30979', 'HIP3125', 'HIP3143', 'HIP32769', 'HIP32892', 'HIP32919', 'HIP33241', 'HIP33955', 'HIP3418', 'HIP35093', 'HIP36338', 'HIP36551', 'HIP36635', 'HIP36834', 'HIP37217', 'HIP37766', 'HIP37798', 'HIP38939', 'HIP38969', 'HIP3998', 'HIP40375', 'HIP40910', 'HIP41130', 'HIP41443', 'HIP41689', 'HIP42220', 'HIP42491', 'HIP428', 'HIP43534', 'HIP4353', 'HIP4454', 'HIP45042', 'HIP45839', 'HIP46199', 'HIP46343', 'HIP46417', 'HIP46655', 'HIP46769', 'HIP47098', 'HIP47201', 'HIP47513', 'HIP47650', 'HIP48139', 'HIP48411', 'HIP4845', 'HIP48714', 'HIP48740', 'HIP48855', 'HIP49091', 'HIP49197', 'HIP5004', 'HIP50341', 'HIP50960', 'HIP51007', 'HIP51443', 'HIP5215', 'HIP5247', 'HIP52942A', 'HIP53020', 'HIP53327', 'HIP53541', 'HIP54459', 'HIP54498', 'HIP54532', 'HIP54651', 'HIP54810', 'HIP55360', 'HIP55507', 'HIP55915', 'HIP5643', 'HIP56630', 'HIP5663', 'HIP57050', 'HIP57058', 'HIP57087', 'HIP57274', 'HIP5741', 'HIP57450', 'HIP57493', 'HIP57548', 'HIP5938', 'HIP59406B', 'HIP59406', 'HIP59496', 'HIP59748', 'HIP60093', 'HIP60357', 'HIP60559', 'HIP60633', 'HIP61205', 'HIP61706', 'HIP62406', 'HIP6276', 'HIP62794', 'HIP62847', 'HIP63257', 'HIP6339', 'HIP6344', 'HIP63510', 'HIP63759', 'HIP63762', 'HIP63894', 'HIP64048', 'HIP64262', 'HIP65016', 'HIP66074', 'HIP66193', 'HIP66222', 'HIP66283', 'HIP66459', 'HIP66840', 'HIP67164', 'HIP67691', 'HIP67842', 'HIP70865', 'HIP70975', 'HIP71253', 'HIP71898', 'HIP73427', 'HIP74346', 'HIP74981', 'HIP74995', 'HIP75672', 'HIP77908', 'HIP78184', 'HIP7830', 'HIP78423', 'HIP78999', 'HIP79308', 'HIP79431', 'HIP79698', 'HIP80096', 'HIP80295', 'HIP8051', 'HIP8070', 'HIP80824', 'HIP83043', 'HIP8361', 'HIP83762', 'HIP84099', 'HIP84790', 'HIP85582', 'HIP85977', 'HIP86961', 'HIP87062', 'HIP87123', 'HIP87464', 'HIP87579', 'HIP87803', 'HIP89087', 'HIP89215', 'HIP90376', 'HIP91605', 'HIP91699', 'HIP92403', 'HIP92922', 'HIP93119', 'HIP93703', 'HIP93871', 'HIP94931', 'HIP96561', 'HIP97051', 'HIP9788', 'HIP98381', 'HIP99205', 'HIP99332', 'HIP99385', 'HTR125-001', 'HTR127-008', 'HTR133-004', 'HTR136-001', 'HTR145-001', 'HTR145-002', 'HTR152-001', 'HTR152-004', 'HTR153-004', 'HTR154-011', 'HTR155-001', 'HTR161-003', 'HTR161-009', 'HTR169-024', 'HTR170-004', 'HTR176-002', 'HTR182-001', 'HTR185-002', 'HTR188-002', 'HTR191-001', 'HTR194-006', 'HTR195-003', 'HTR196-004', 'HTR198-002', 'HTR204-007', 'HTR204-010', 'HTR204-011', 'HTR204-014', 'HTR205-024', 'HTR205-22', 'HTR205-23', 'HTR213-001', 'HTR239-001', 'HTR239-004', 'HTR248-002', 'LHS462', 'RXJ0348.9', 'RXJ0434.3', 'S101438B', 'S11844', 'S122446', 'S130811', 'S92823', 'T-CYG0-04140', 'T-LAC0-01033', 'T-LAC0-14888', 'T-ORI0-01816', 'T-ORI0-02163', 'T-PER2-04344', 'TRES1', 'TRES2', 'TRES3', 'TRES4', 'V383LAC', 'VULCAN8842', 'WASP-1', 'WASP14', 'WASP234318', 'XO-2']


 
