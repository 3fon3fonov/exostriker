#!/usr/bin/python


def param_gui(self):
    
    param_gui = [
            self.K1, self.P1, self.e1, self.om1, self.ma1, self.incl1, self.Omega1,
            self.K2, self.P2, self.e2, self.om2, self.ma2, self.incl2, self.Omega2,
            self.K3, self.P3, self.e3, self.om3, self.ma3, self.incl3, self.Omega3,
            self.K4, self.P4, self.e4, self.om4, self.ma4, self.incl4, self.Omega4, 
            self.K5, self.P5, self.e5, self.om5, self.ma5, self.incl5, self.Omega5,
            self.K6, self.P6, self.e6, self.om6, self.ma6, self.incl6, self.Omega6,
            self.K7, self.P7, self.e7, self.om7, self.ma7, self.incl7, self.Omega7, 
            self.K8, self.P8, self.e8, self.om8, self.ma8, self.incl8, self.Omega8,
            self.K9, self.P9, self.e9, self.om9, self.ma9, self.incl9, self.Omega9,
            ]

    return param_gui

def param_gui_wd(self):
    
    param_gui_wd = [
            self.om_dot_1, self.om_dot_2, self.om_dot_3, 
            self.om_dot_4, self.om_dot_5, self.om_dot_6, 
            self.om_dot_7, self.om_dot_8, self.om_dot_9
            ]
    return param_gui_wd

def param_gui_tr(self):
    
    param_gui_tr = [
             self.t0_1, self.pl_rad_1, self.a_sol_1,
             self.t0_2, self.pl_rad_2, self.a_sol_2,
             self.t0_3, self.pl_rad_3, self.a_sol_3,
             self.t0_4, self.pl_rad_4, self.a_sol_4, 
             self.t0_5, self.pl_rad_5, self.a_sol_5,
             self.t0_6, self.pl_rad_6, self.a_sol_6,
             self.t0_7, self.pl_rad_7, self.a_sol_7, 
             self.t0_8, self.pl_rad_8, self.a_sol_8,
             self.t0_9, self.pl_rad_9, self.a_sol_9,
             ]
    return param_gui_tr






def rvs_data_gui(self):
        rvs_data_gui = [
                self.Data1,self.Data2,self.Data3,self.Data4,self.Data5,
                self.Data6,self.Data7,self.Data8,self.Data9,self.Data10
                ]
        return rvs_data_gui

def rvs_data_jitter_gui(self):        
        rvs_data_jitter_gui = [
                self.jitter_Data1,self.jitter_Data2,self.jitter_Data3,self.jitter_Data4,self.jitter_Data5,
                self.jitter_Data6,self.jitter_Data7,self.jitter_Data8,self.jitter_Data9,self.jitter_Data10
                ]
        return rvs_data_jitter_gui


def tra_data_gui(self):        
        tra_data_gui = [
                self.trans_Data1,self.trans_Data2,self.trans_Data3,self.trans_Data4,self.trans_Data5,
                self.trans_Data6,self.trans_Data7,self.trans_Data8,self.trans_Data9,self.trans_Data10
                ]
        return tra_data_gui
               
def tra_data_jitter_gui(self):        
        tra_data_jitter_gui = [
                self.jitter_trans_Data1,self.jitter_trans_Data2,self.jitter_trans_Data3,self.jitter_trans_Data4,self.jitter_trans_Data5,
                self.jitter_trans_Data6,self.jitter_trans_Data7,self.jitter_trans_Data8,self.jitter_trans_Data9,self.jitter_trans_Data10
                ]
        return tra_data_jitter_gui
        
            



def param_bounds_gui(self):
  
    param_bounds_gui = [
    [self.K_min_1,self.K_max_1],[self.P_min_1,self.P_max_1], [self.e_min_1,self.e_max_1],[self.om_min_1,self.om_max_1], [self.ma_min_1,self.ma_max_1],[self.incl_min_1,self.incl_max_1], [self.Omega_min_1,self.Omega_max_1],[self.t0_min_1,self.t0_max_1],[self.pl_rad_min_1,self.pl_rad_max_1],[self.a_sol_min_1,self.a_sol_max_1],
    [self.K_min_2,self.K_max_2],[self.P_min_2,self.P_max_2], [self.e_min_2,self.e_max_2],[self.om_min_2,self.om_max_2], [self.ma_min_2,self.ma_max_2],[self.incl_min_2,self.incl_max_2], [self.Omega_min_2,self.Omega_max_2],[self.t0_min_2,self.t0_max_2],[self.pl_rad_min_2,self.pl_rad_max_2],[self.a_sol_min_2,self.a_sol_max_2],
    [self.K_min_3,self.K_max_3],[self.P_min_3,self.P_max_3], [self.e_min_3,self.e_max_3],[self.om_min_3,self.om_max_3], [self.ma_min_3,self.ma_max_3],[self.incl_min_3,self.incl_max_3], [self.Omega_min_3,self.Omega_max_3],[self.t0_min_3,self.t0_max_3],[self.pl_rad_min_3,self.pl_rad_max_3],[self.a_sol_min_3,self.a_sol_max_3],
    [self.K_min_4,self.K_max_4],[self.P_min_4,self.P_max_4], [self.e_min_4,self.e_max_4],[self.om_min_4,self.om_max_4], [self.ma_min_4,self.ma_max_4],[self.incl_min_4,self.incl_max_4], [self.Omega_min_4,self.Omega_max_4],[self.t0_min_4,self.t0_max_4],[self.pl_rad_min_4,self.pl_rad_max_4],[self.a_sol_min_4,self.a_sol_max_4],
    [self.K_min_5,self.K_max_5],[self.P_min_5,self.P_max_5], [self.e_min_5,self.e_max_5],[self.om_min_5,self.om_max_5], [self.ma_min_5,self.ma_max_5],[self.incl_min_5,self.incl_max_5], [self.Omega_min_5,self.Omega_max_5],[self.t0_min_5,self.t0_max_5],[self.pl_rad_min_5,self.pl_rad_max_5],[self.a_sol_min_5,self.a_sol_max_5],
    [self.K_min_6,self.K_max_6],[self.P_min_6,self.P_max_6], [self.e_min_6,self.e_max_6],[self.om_min_6,self.om_max_6], [self.ma_min_6,self.ma_max_6],[self.incl_min_6,self.incl_max_6], [self.Omega_min_6,self.Omega_max_6],[self.t0_min_6,self.t0_max_6],[self.pl_rad_min_6,self.pl_rad_max_6],[self.a_sol_min_6,self.a_sol_max_6],
    [self.K_min_7,self.K_max_7],[self.P_min_7,self.P_max_7], [self.e_min_7,self.e_max_7],[self.om_min_7,self.om_max_7], [self.ma_min_7,self.ma_max_7],[self.incl_min_7,self.incl_max_7], [self.Omega_min_7,self.Omega_max_7],[self.t0_min_7,self.t0_max_7],[self.pl_rad_min_7,self.pl_rad_max_7],[self.a_sol_min_7,self.a_sol_max_7],
    [self.K_min_8,self.K_max_8],[self.P_min_8,self.P_max_8], [self.e_min_8,self.e_max_8],[self.om_min_8,self.om_max_8], [self.ma_min_8,self.ma_max_8],[self.incl_min_8,self.incl_max_8], [self.Omega_min_8,self.Omega_max_8],[self.t0_min_8,self.t0_max_8],[self.pl_rad_min_8,self.pl_rad_max_8],[self.a_sol_min_8,self.a_sol_max_8],
    [self.K_min_9,self.K_max_9],[self.P_min_9,self.P_max_9], [self.e_min_9,self.e_max_9],[self.om_min_9,self.om_max_9], [self.ma_min_9,self.ma_max_9],[self.incl_min_9,self.incl_max_9], [self.Omega_min_9,self.Omega_max_9],[self.t0_min_9,self.t0_max_9],[self.pl_rad_min_9,self.pl_rad_max_9],[self.a_sol_min_9,self.a_sol_max_9]               
    ]
 
    return param_bounds_gui




def offset_bounds_gui(self):
    
    offset_bounds_gui = [
    [self.Data1_min,self.Data1_max], [self.Data2_min,self.Data2_max], [self.Data3_min,self.Data3_max], [self.Data4_min,self.Data4_max], [self.Data5_min,self.Data5_max],   
    [self.Data6_min,self.Data6_max], [self.Data7_min,self.Data7_max], [self.Data8_min,self.Data8_max], [self.Data9_min,self.Data9_max], [self.Data10_min,self.Data10_max]
    ]
    
    return offset_bounds_gui        
        


def jitter_bounds_gui(self):  
      
    jitter_bounds_gui = [
    [self.jitter1_min,self.jitter1_max], [self.jitter2_min,self.jitter2_max], [self.jitter3_min,self.jitter3_max], [self.jitter4_min,self.jitter4_max], [self.jitter5_min,self.jitter5_max],   
    [self.jitter6_min,self.jitter6_max], [self.jitter7_min,self.jitter7_max], [self.jitter8_min,self.jitter8_max], [self.jitter9_min,self.jitter9_max], [self.jitter10_min,self.Data10_max]   
    ]  
    
    return jitter_bounds_gui    















    
        