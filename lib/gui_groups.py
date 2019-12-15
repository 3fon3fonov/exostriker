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

def param_errors_gui(self):

    param_errors_gui = [self.err_K1,self.err_P1,self.err_e1,self.err_om1,self.err_ma1, self.err_i1, self.err_Om1,
                        self.err_K2,self.err_P2,self.err_e2,self.err_om2,self.err_ma2, self.err_i2, self.err_Om2,
                        self.err_K3,self.err_P3,self.err_e3,self.err_om3,self.err_ma3, self.err_i3, self.err_Om3,
                        self.err_K4,self.err_P4,self.err_e4,self.err_om4,self.err_ma4, self.err_i4, self.err_Om4,  
                        self.err_K5,self.err_P5,self.err_e5,self.err_om5,self.err_ma5, self.err_i5, self.err_Om5,
                        self.err_K6,self.err_P6,self.err_e6,self.err_om6,self.err_ma6, self.err_i6, self.err_Om6,
                        self.err_K7,self.err_P7,self.err_e7,self.err_om7,self.err_ma7, self.err_i7, self.err_Om7, 
                        self.err_K8,self.err_P8,self.err_e8,self.err_om8,self.err_ma8, self.err_i8, self.err_Om8,
                        self.err_K9,self.err_P9,self.err_e9,self.err_om9,self.err_ma9, self.err_i9, self.err_Om9,                       
                        ]
    return param_errors_gui

def use_param_gui(self):

    use_param_gui  = [self.use_K1, self.use_P1, self.use_e1, self.use_om1, self.use_ma1, self.use_incl1, self.use_Omega1,
                      self.use_K2, self.use_P2, self.use_e2, self.use_om2, self.use_ma2, self.use_incl2, self.use_Omega2,
                      self.use_K3, self.use_P3, self.use_e3, self.use_om3, self.use_ma3, self.use_incl3, self.use_Omega3,                        
                      self.use_K4, self.use_P4, self.use_e4, self.use_om4, self.use_ma4, self.use_incl4, self.use_Omega4,    
                      self.use_K5, self.use_P5, self.use_e5, self.use_om5, self.use_ma5, self.use_incl5, self.use_Omega5,    
                      self.use_K6, self.use_P6, self.use_e6, self.use_om6, self.use_ma6, self.use_incl6, self.use_Omega6, 
                      self.use_K7, self.use_P7, self.use_e7, self.use_om7, self.use_ma7, self.use_incl7, self.use_Omega7,    
                      self.use_K8, self.use_P8, self.use_e8, self.use_om8, self.use_ma8, self.use_incl8, self.use_Omega8,    
                      self.use_K9, self.use_P9, self.use_e9, self.use_om9, self.use_ma9, self.use_incl9, self.use_Omega9,                       
                      ]
    return use_param_gui




###########################################################################



def param_gui_wd(self):
    
    param_gui_wd = [
            self.om_dot_1, self.om_dot_2, self.om_dot_3, 
            self.om_dot_4, self.om_dot_5, self.om_dot_6, 
            self.om_dot_7, self.om_dot_8, self.om_dot_9
            ]
    return param_gui_wd

def use_param_gui_wd(self):

    use_param_gui_wd = [
            self.use_om_dot_1, self.use_om_dot_2, self.use_om_dot_3, 
            self.use_om_dot_4, self.use_om_dot_5, self.use_om_dot_6, 
            self.use_om_dot_7, self.use_om_dot_8, self.use_om_dot_9
            ]
    return use_param_gui_wd

def param_errors_gui_wd(self):

    param_errors_gui_wd = [
                self.err_om_dot_1,self.err_om_dot_2,self.err_om_dot_3,
                self.err_om_dot_4,self.err_om_dot_5,self.err_om_dot_6,
                self.err_om_dot_7,self.err_om_dot_8,self.err_om_dot_9,
            ]
    return param_errors_gui_wd

###########################################################################



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



def use_param_gui_tr(self):

    use_param_gui_tr = [self.use_t0_1, self.use_pl_rad_1, self.use_a_sol_1,
             self.use_t0_2, self.use_pl_rad_2, self.use_a_sol_2,
             self.use_t0_3, self.use_pl_rad_3, self.use_a_sol_3,
             self.use_t0_4, self.use_pl_rad_4, self.use_a_sol_4, 
             self.use_t0_5, self.use_pl_rad_5, self.use_a_sol_5,
             self.use_t0_6, self.use_pl_rad_6, self.use_a_sol_6,
             self.use_t0_7, self.use_pl_rad_7, self.use_a_sol_7, 
             self.use_t0_8, self.use_pl_rad_8, self.use_a_sol_8,
             self.use_t0_9, self.use_pl_rad_9, self.use_a_sol_9,
             ]


    return use_param_gui_tr

###########################################################################


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


def use_data_offset_gui(self):

    use_data_offset_gui = [self.use_offset_Data1,self.use_offset_Data2,self.use_offset_Data3,self.use_offset_Data4,
                           self.use_offset_Data5,self.use_offset_Data6,self.use_offset_Data7,self.use_offset_Data8,
                           self.use_offset_Data9,self.use_offset_Data10]
    
    return use_data_offset_gui


def use_data_jitter_gui(self):

    use_data_jitter_gui = [self.use_jitter_Data1,self.use_jitter_Data2,self.use_jitter_Data3,self.use_jitter_Data4,self.use_jitter_Data5,
                           self.use_jitter_Data6,self.use_jitter_Data7,self.use_jitter_Data8,self.use_jitter_Data9,self.use_jitter_Data10]
    
    return use_data_jitter_gui


def data_errors_gui(self):

    data_errors_gui = [
            self.err_Data1,self.err_Data2,self.err_Data3,self.err_Data4,self.err_Data5,
            self.err_Data6,self.err_Data7,self.err_Data8,self.err_Data9,self.err_Data10
            ]
    
    return data_errors_gui


def data_errors_jitter_gui(self):

    data_errors_jitter_gui = [
            self.err_jitter_Data1,self.err_jitter_Data2,self.err_jitter_Data3,self.err_jitter_Data4,self.err_jitter_Data5,
            self.err_jitter_Data6,self.err_jitter_Data7,self.err_jitter_Data8,self.err_jitter_Data9,self.err_jitter_Data10
            ]

    return data_errors_jitter_gui


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


def use_tra_data_offset_gui(self):

    use_tra_data_offset_gui = [
            self.use_offset_trans_Data1,self.use_offset_trans_Data2,self.use_offset_trans_Data3,self.use_offset_trans_Data4,
            self.use_offset_trans_Data5,self.use_offset_trans_Data6,self.use_offset_trans_Data7,self.use_offset_trans_Data8,
            self.use_offset_trans_Data9,self.use_offset_trans_Data10
            ]

    return use_tra_data_offset_gui


def use_tra_data_jitter_gui(self):
    
    use_tra_data_jitter_gui = [
            self.use_jitter_trans_Data1,self.use_jitter_trans_Data2,self.use_jitter_trans_Data3,self.use_jitter_trans_Data4,
            self.use_jitter_trans_Data5,self.use_jitter_trans_Data6,self.use_jitter_trans_Data7,self.use_jitter_trans_Data8,
            self.use_jitter_trans_Data9,self.use_jitter_trans_Data10
            ]
    
    return use_tra_data_jitter_gui


def tra_data_errors_gui(self):

    tra_data_errors_gui = [
            self.err_trans_Data1,self.err_trans_Data2,self.err_trans_Data3,self.err_trans_Data4,self.err_trans_Data5,
            self.err_trans_Data6,self.err_trans_Data7,self.err_trans_Data8,self.err_trans_Data9,self.err_trans_Data10
            ]
    
    return tra_data_errors_gui

def tra_data_errors_jitter_gui(self):

    tra_data_errors_jitter_gui = [
            self.err_jitter_trans_Data1,self.err_jitter_trans_Data2,self.err_jitter_trans_Data3,self.err_jitter_trans_Data4,
            self.err_jitter_trans_Data5,self.err_jitter_trans_Data6,self.err_jitter_trans_Data7,self.err_jitter_trans_Data8,
            self.err_jitter_trans_Data9,self.err_jitter_trans_Data10
            ]
    
    return tra_data_errors_jitter_gui






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


################### OmDot ########################

def om_dot_bounds_gui(self):

    om_dot_bounds_gui = [
    [self.omega_dot_min_1,self.omega_dot_max_1], [self.omega_dot_min_2,self.omega_dot_max_2], 
    [self.omega_dot_min_3,self.omega_dot_max_3], [self.omega_dot_min_4,self.omega_dot_max_4], 
    [self.omega_dot_min_5,self.omega_dot_max_5], [self.omega_dot_min_6,self.omega_dot_max_6], 
    [self.omega_dot_min_7,self.omega_dot_max_7], [self.omega_dot_min_8,self.omega_dot_max_8], 
    [self.omega_dot_min_9,self.omega_dot_max_9] 
    ]  
    return om_dot_bounds_gui



################### GP ########################

def gp_rot_params(self):

    gp_rot_params = [
            self.GP_rot_kernel_Amp,
            self.GP_rot_kernel_time_sc,
            self.GP_rot_kernel_Per,
            self.GP_rot_kernel_fact
            ]

    return gp_rot_params

def gp_rot_errors_gui(self):

    gp_rot_errors_gui = [
            self.err_rot_kernel_Amp,
            self.err_rot_kernel_time_sc,
            self.err_rot_kernel_Per,
            self.err_rot_kernel_fact
            ]
        
    return gp_rot_errors_gui

def use_gp_rot_params(self):
    use_gp_rot_params = [
            self.use_GP_rot_kernel_Amp,
            self.use_GP_rot_kernel_time_sc,
            self.use_GP_rot_kernel_Per,
            self.use_GP_rot_kernel_fact
            ]
    return use_gp_rot_params


def gp_sho_params(self):

    gp_sho_params = [
            self.GP_sho_kernel_S,
            self.GP_sho_kernel_Q,
            self.GP_sho_kernel_omega
            ]

    return gp_sho_params

def use_gp_sho_params(self):

    use_gp_sho_params = [
            self.use_GP_sho_kernel_S,
            self.use_GP_sho_kernel_Q,
            self.use_GP_sho_kernel_omega
            ]
    return use_gp_sho_params


def gp_sho_errors_gui(self):

    gp_sho_errors_gui = [
            self.err_sho_kernel_S,
            self.err_sho_kernel_Q,
            self.err_sho_kernel_omega
            ]
    
    return gp_sho_errors_gui



def tra_gp_rot_params(self):

    tra_gp_rot_params = [
            self.tra_GP_rot_kernel_Amp,
            self.tra_GP_rot_kernel_time_sc,
            self.tra_GP_rot_kernel_Per,
            self.tra_GP_rot_kernel_fact]

    return tra_gp_rot_params


def use_tra_gp_rot_params(self):

    use_tra_gp_rot_params = [
            self.use_tra_GP_rot_kernel_Amp,
            self.use_tra_GP_rot_kernel_time_sc,
            self.use_tra_GP_rot_kernel_Per,
            self.use_tra_GP_rot_kernel_fact
            ]

    return use_tra_gp_rot_params



def tra_gp_sho_params(self):

    tra_gp_sho_params = [
            self.tra_GP_sho_kernel_S,
            self.tra_GP_sho_kernel_Q,
            self.tra_GP_sho_kernel_omega]

    return tra_gp_sho_params

def use_tra_gp_sho_params(self):

    use_tra_gp_sho_params = [
            self.use_tra_GP_sho_kernel_S,
            self.use_tra_GP_sho_kernel_Q,
            self.use_tra_GP_sho_kernel_omega]

    return use_tra_gp_sho_params


################ labels ##########################

def param_a_gui(self): 

    param_a_gui = [
            self.label_a1, self.label_a2, self.label_a3, 
            self.label_a4, self.label_a5, self.label_a6, 
            self.label_a7, self.label_a8, self.label_a9
            ]

    return param_a_gui

def param_mass_gui(self):

    param_mass_gui = [
            self.label_mass1, self.label_mass2, self.label_mass3, 
            self.label_mass4, self.label_mass5, self.label_mass6, 
            self.label_mass7, self.label_mass8, self.label_mass9
            ]

    return param_mass_gui


def param_t_peri_gui(self): 

    param_t_peri_gui = [
            self.label_t_peri1, self.label_t_peri2, self.label_t_peri3, 
            self.label_t_peri4, self.label_t_peri5, self.label_t_peri6, 
            self.label_t_peri7, self.label_t_peri8, self.label_t_peri9
            ]

    return param_t_peri_gui


def planet_checked_gui(self): 

    planet_checked_gui = [
            self.use_Planet1,self.use_Planet2,self.use_Planet3,
            self.use_Planet4,self.use_Planet5,self.use_Planet6,
            self.use_Planet7,self.use_Planet8,self.use_Planet9
            ]

    return planet_checked_gui



        

def add_rv_error(self): 
    add_rv_error = [
            [self.inflate_RV_sigma_1,self.use_inflate_RV_sigma_1],[self.inflate_RV_sigma_2,self.use_inflate_RV_sigma_2],
            [self.inflate_RV_sigma_3,self.use_inflate_RV_sigma_3],[self.inflate_RV_sigma_4,self.use_inflate_RV_sigma_4],
            [self.inflate_RV_sigma_5,self.use_inflate_RV_sigma_5],[self.inflate_RV_sigma_6,self.use_inflate_RV_sigma_6],
            [self.inflate_RV_sigma_7,self.use_inflate_RV_sigma_7],[self.inflate_RV_sigma_8,self.use_inflate_RV_sigma_8],
            [self.inflate_RV_sigma_9,self.use_inflate_RV_sigma_9],[self.inflate_RV_sigma_10,self.use_inflate_RV_sigma_10]
            ]
 
    return add_rv_error



######################### Arb N-body #################################
    
def arb_param_gui(self): 

    arb_param_gui = [
                 self.arb_K_1, self.arb_P_1, self.arb_e_1, self.arb_om_1, self.arb_ma_1, self.arb_incl_1, self.arb_Om_1,
                 self.arb_K_2, self.arb_P_2, self.arb_e_2, self.arb_om_2, self.arb_ma_2, self.arb_incl_2, self.arb_Om_2,
                 self.arb_K_3, self.arb_P_3, self.arb_e_3, self.arb_om_3, self.arb_ma_3, self.arb_incl_3, self.arb_Om_3,
                 self.arb_K_4, self.arb_P_4, self.arb_e_4, self.arb_om_4, self.arb_ma_4, self.arb_incl_4, self.arb_Om_4, 
                 self.arb_K_5, self.arb_P_5, self.arb_e_5, self.arb_om_5, self.arb_ma_5, self.arb_incl_5, self.arb_Om_5,
                 self.arb_K_6, self.arb_P_6, self.arb_e_6, self.arb_om_6, self.arb_ma_6, self.arb_incl_6, self.arb_Om_6,
                 self.arb_K_7, self.arb_P_7, self.arb_e_7, self.arb_om_7, self.arb_ma_7, self.arb_incl_7, self.arb_Om_7, 
                 self.arb_K_8, self.arb_P_8, self.arb_e_8, self.arb_om_8, self.arb_ma_8, self.arb_incl_8, self.arb_Om_8,
                 self.arb_K_9, self.arb_P_9, self.arb_e_9, self.arb_om_9, self.arb_ma_9, self.arb_incl_9, self.arb_Om_9,
                 ]
    return arb_param_gui


def arb_param_gui_use(self): 
    arb_param_gui_use = [
            self.use_arb_Planet_1,self.use_arb_Planet_2,self.use_arb_Planet_3,
            self.use_arb_Planet_4,self.use_arb_Planet_5,self.use_arb_Planet_6,
            self.use_arb_Planet_7,self.use_arb_Planet_8,self.use_arb_Planet_9
            ]
 
    return arb_param_gui_use




def ttv_data_to_planet(self): 
    
    ttv_data_to_planet = [
            self.ttv_data_planet_1,self.ttv_data_planet_2,self.ttv_data_planet_3,self.ttv_data_planet_4,self.ttv_data_planet_5,
            self.ttv_data_planet_6,self.ttv_data_planet_7,self.ttv_data_planet_8,self.ttv_data_planet_9,self.ttv_data_planet_10,
            ]
   
    return ttv_data_to_planet


def use_ttv_data_to_planet(self): 
    
    use_ttv_data_to_planet = [
            self.use_ttv_data_1,self.use_ttv_data_2,self.use_ttv_data_3,self.use_ttv_data_4,self.use_ttv_data_5,
            self.use_ttv_data_6,self.use_ttv_data_7,self.use_ttv_data_8,self.use_ttv_data_9,self.use_ttv_data_10,
            ]
   
    return use_ttv_data_to_planet
