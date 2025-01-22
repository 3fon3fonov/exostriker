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

    param_errors_gui = [self.err_K1,self.err_P1,self.err_e1,self.err_om1,self.err_ma1, self.err_incl1, self.err_Om1,
                        self.err_K2,self.err_P2,self.err_e2,self.err_om2,self.err_ma2, self.err_incl2, self.err_Om2,
                        self.err_K3,self.err_P3,self.err_e3,self.err_om3,self.err_ma3, self.err_incl3, self.err_Om3,
                        self.err_K4,self.err_P4,self.err_e4,self.err_om4,self.err_ma4, self.err_incl4, self.err_Om4,  
                        self.err_K5,self.err_P5,self.err_e5,self.err_om5,self.err_ma5, self.err_incl5, self.err_Om5,
                        self.err_K6,self.err_P6,self.err_e6,self.err_om6,self.err_ma6, self.err_incl6, self.err_Om6,
                        self.err_K7,self.err_P7,self.err_e7,self.err_om7,self.err_ma7, self.err_incl7, self.err_Om7, 
                        self.err_K8,self.err_P8,self.err_e8,self.err_om8,self.err_ma8, self.err_incl8, self.err_Om8,
                        self.err_K9,self.err_P9,self.err_e9,self.err_om9,self.err_ma9, self.err_incl9, self.err_Om9,                       
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

def param_gui_K(self):
    
    param_gui_K = [
            self.K1,
            self.K2,
            self.K3,
            self.K4,
            self.K5,
            self.K6, 
            self.K7, 
            self.K8, 
            self.K9
            ]

    return param_gui_K

def param_errors_gui_K(self):

    param_errors_gui_K = [self.err_K1,
                        self.err_K2,
                        self.err_K3,
                        self.err_K4,
                        self.err_K5,
                        self.err_K6,
                        self.err_K7,
                        self.err_K8,
                        self.err_K9                 
                        ]
    return param_errors_gui_K

def use_param_gui_K(self):

    use_param_gui_K  = [self.use_K1,
                      self.use_K2, 
                      self.use_K3,                
                      self.use_K4, 
                      self.use_K5, 
                      self.use_K6,
                      self.use_K7, 
                      self.use_K8, 
                      self.use_K9,                  
                      ]
    return use_param_gui_K

###########################################################################

def param_gui_P(self):
    
    param_gui_P = [
            self.P1, 
            self.P2, 
            self.P3,
            self.P4, 
            self.P5, 
            self.P6, 
            self.P7, 
            self.P8, 
            self.P9
            ]

    return param_gui_P

def param_errors_gui_P(self):

    param_errors_gui_P = [self.err_P1,
                        self.err_P2,
                        self.err_P3,
                        self.err_P4,  
                        self.err_P5,
                        self.err_P6,
                        self.err_P7,
                        self.err_P8,
                        self.err_P9                      
                        ]
    return param_errors_gui_P

def use_param_gui_P(self):

    use_param_gui_P  = [self.use_P1,  
                      self.use_P2,  
                      self.use_P3,                       
                      self.use_P4,     
                      self.use_P5,  
                      self.use_P6,  
                      self.use_P7,      
                      self.use_P8,  
                      self.use_P9                
                      ]
    return use_param_gui_P

###########################################################################

def param_gui_e(self):
    
    param_gui_e = [
            self.e1, 
            self.e2, 
            self.e3,
            self.e4, 
            self.e5, 
            self.e6, 
            self.e7, 
            self.e8, 
            self.e9
            ]

    return param_gui_e

def param_errors_gui_e(self):

    param_errors_gui_e = [self.err_e1,
                        self.err_e2,
                        self.err_e3,
                        self.err_e4,  
                        self.err_e5,
                        self.err_e6,
                        self.err_e7,
                        self.err_e8,
                        self.err_e9                      
                        ]
    return param_errors_gui_e

def use_param_gui_e(self):

    use_param_gui_e  = [self.use_e1,  
                      self.use_e2,  
                      self.use_e3,                       
                      self.use_e4,     
                      self.use_e5,  
                      self.use_e6,  
                      self.use_e7,      
                      self.use_e8,  
                      self.use_e9                
                      ]
    return use_param_gui_e

###########################################################################

def param_gui_om(self):
    
    param_gui_om = [
            self.om1, 
            self.om2, 
            self.om3,
            self.om4, 
            self.om5, 
            self.om6, 
            self.om7, 
            self.om8, 
            self.om9
            ]

    return param_gui_om

def param_errors_gui_om(self):

    param_errors_gui_om = [self.err_om1,
                        self.err_om2,
                        self.err_om3,
                        self.err_om4,  
                        self.err_om5,
                        self.err_om6,
                        self.err_om7,
                        self.err_om8,
                        self.err_om9                      
                        ]
    return param_errors_gui_om

def use_param_gui_om(self):

    use_param_gui_om  = [self.use_om1,  
                      self.use_om2,  
                      self.use_om3,                       
                      self.use_om4,     
                      self.use_om5,  
                      self.use_om6,  
                      self.use_om7,      
                      self.use_om8,  
                      self.use_om9                
                      ]
    return use_param_gui_om

###########################################################################

def param_gui_ma(self):
    
    param_gui_ma = [
            self.ma1, 
            self.ma2, 
            self.ma3,
            self.ma4, 
            self.ma5, 
            self.ma6, 
            self.ma7, 
            self.ma8, 
            self.ma9
            ]

    return param_gui_ma

def param_errors_gui_ma(self):

    param_errors_gui_ma = [self.err_ma1,
                        self.err_ma2,
                        self.err_ma3,
                        self.err_ma4,  
                        self.err_ma5,
                        self.err_ma6,
                        self.err_ma7,
                        self.err_ma8,
                        self.err_ma9                      
                        ]
    return param_errors_gui_ma

def use_param_gui_ma(self):

    use_param_gui_ma  = [self.use_ma1,  
                      self.use_ma2,  
                      self.use_ma3,                       
                      self.use_ma4,     
                      self.use_ma5,  
                      self.use_ma6,  
                      self.use_ma7,      
                      self.use_ma8,  
                      self.use_ma9                
                      ]
    return use_param_gui_ma

###########################################################################

def param_gui_incl(self):
    
    param_gui_incl = [
            self.incl1, 
            self.incl2, 
            self.incl3,
            self.incl4, 
            self.incl5, 
            self.incl6, 
            self.incl7, 
            self.incl8, 
            self.incl9
            ]

    return param_gui_incl

def param_errors_gui_incl(self):

    param_errors_gui_incl = [self.err_incl1,
                        self.err_incl2,
                        self.err_incl3,
                        self.err_incl4,  
                        self.err_incl5,
                        self.err_incl6,
                        self.err_incl7,
                        self.err_incl8,
                        self.err_incl9                      
                        ]
    return param_errors_gui_incl

def use_param_gui_incl(self):

    use_param_gui_incl  = [self.use_incl1,  
                      self.use_incl2,  
                      self.use_incl3,                       
                      self.use_incl4,     
                      self.use_incl5,  
                      self.use_incl6,  
                      self.use_incl7,      
                      self.use_incl8,  
                      self.use_incl9                
                      ]
    return use_param_gui_incl

###########################################################################

def param_gui_Omega(self):
    
    param_gui_Omega = [
            self.Omega1, 
            self.Omega2, 
            self.Omega3,
            self.Omega4, 
            self.Omega5, 
            self.Omega6, 
            self.Omega7, 
            self.Omega8, 
            self.Omega9
            ]

    return param_gui_Omega

def param_errors_gui_Omega(self):

    param_errors_gui_Omega = [self.err_Omega1,
                        self.err_Omega2,
                        self.err_Omega3,
                        self.err_Omega4,  
                        self.err_Omega5,
                        self.err_Omega6,
                        self.err_Omega7,
                        self.err_Omega8,
                        self.err_Omega9                      
                        ]
    return param_errors_gui_Omega

def use_param_gui_Omega(self):

    use_param_gui_Omega  = [self.use_Omega1,  
                      self.use_Omega2,  
                      self.use_Omega3,                       
                      self.use_Omega4,     
                      self.use_Omega5,  
                      self.use_Omega6,  
                      self.use_Omega7,      
                      self.use_Omega8,  
                      self.use_Omega9                
                      ]
    return use_param_gui_Omega


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


def param_gui_ast(self):
    
    param_gui_ast = [
            self.ast_alpha,
            self.ast_delta,
            self.ast_pi,
            self.ast_mu_alpha,       
            self.ast_mu_delta
            ]
    return param_gui_ast

def use_param_gui_ast(self):

    use_param_gui_ast = [
            self.use_ast_alpha,
            self.use_ast_delta,
            self.use_ast_pi,
            self.use_ast_mu_alpha,       
            self.use_ast_mu_delta
            ]
    return use_param_gui_ast

def param_errors_gui_ast(self):

    param_errors_gui_ast = [
            self.err_ast_alpha,
            self.err_ast_delta,
            self.err_ast_pi,
            self.err_ast_mu_alpha,       
            self.err_ast_mu_delta
            ]
    return param_errors_gui_ast



def ast_bounds_gui(self):

    ast_bounds_gui = [
    [self.ast_alpha_min,self.ast_alpha_max],
    [self.ast_delta_min,self.ast_delta_max],
    [self.ast_pi_min,self.ast_pi_max],
    [self.ast_mu_alpha_min,self.ast_mu_alpha_max],
    [self.ast_mu_delta_min,self.ast_mu_delta_max]
    ]  
    return ast_bounds_gui



def ast_norm_pr_gui(self):

    ast_norm_pr_gui = [
    [self.ast_alpha_mean,self.ast_alpha_sigma,self.use_ast_alpha_nr_pr],
    [self.ast_delta_mean,self.ast_delta_sigma,self.use_ast_delta_nr_pr],
    [self.ast_pi_mean,self.ast_pi_sigma,self.use_ast_pi_nr_pr],
    [self.ast_mu_alpha_mean,self.ast_mu_alpha_sigma,self.use_ast_mu_alpha_nr_pr],
    [self.ast_mu_delta_mean,self.ast_mu_delta_sigma,self.use_ast_mu_delta_nr_pr]
    ]
    return ast_norm_pr_gui



def ast_jeff_pr_gui(self):

    ast_jeff_pr_gui = [
    [self.ast_alpha_alpha,self.ast_alpha_beta,self.use_ast_alpha_jeff_pr],
    [self.ast_delta_alpha,self.ast_delta_beta,self.use_ast_delta_jeff_pr],
    [self.ast_pi_alpha,self.ast_pi_beta,self.use_ast_pi_jeff_pr],
    [self.ast_mu_alpha_alpha,self.ast_mu_alpha_beta,self.use_ast_mu_alpha_jeff_pr],
    [self.ast_mu_delta_alpha,self.ast_mu_delta_beta,self.use_ast_mu_delta_jeff_pr]
    ]
    return ast_jeff_pr_gui


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


def err_t0(self):

    err_t0 = [self.err_t0_1,self.err_t0_2,self.err_t0_3,
              self.err_t0_4,self.err_t0_5,self.err_t0_6,
              self.err_t0_7,self.err_t0_8,self.err_t0_9,
            ]
    return err_t0


def err_pl_rad(self):

    err_pl_rad = [self.err_pl_rad_1,self.err_pl_rad_2,self.err_pl_rad_3,
                  self.err_pl_rad_4,self.err_pl_rad_5,self.err_pl_rad_6,
                  self.err_pl_rad_7,self.err_pl_rad_8,self.err_pl_rad_9,
            ]
    return err_pl_rad


def err_a_sol(self):

    err_a_sol = [self.err_a_sol_1,self.err_a_sol_2,self.err_a_sol_3,
                 self.err_a_sol_4,self.err_a_sol_5,self.err_a_sol_6,
                 self.err_a_sol_7,self.err_a_sol_8,self.err_a_sol_9,
            ]
    return err_a_sol



###########################  RV data #########################################


def rvs_data_gui(self):
    rvs_data_gui = [
            self.Data1,self.Data2,self.Data3,self.Data4,self.Data5,
            self.Data6,self.Data7,self.Data8,self.Data9,self.Data10,
            self.Data11,self.Data12,self.Data13,self.Data14,self.Data15,
            self.Data16,self.Data17,self.Data18,self.Data19,self.Data20
            ]
    return rvs_data_gui



def rvs_data_jitter_gui(self):
    rvs_data_jitter_gui = [
            self.jitter_Data1,self.jitter_Data2,self.jitter_Data3,self.jitter_Data4,self.jitter_Data5,
            self.jitter_Data6,self.jitter_Data7,self.jitter_Data8,self.jitter_Data9,self.jitter_Data10,
            self.jitter_Data11,self.jitter_Data12,self.jitter_Data13,self.jitter_Data14,self.jitter_Data15,
            self.jitter_Data16,self.jitter_Data17,self.jitter_Data18,self.jitter_Data19,self.jitter_Data20
            ]
    return rvs_data_jitter_gui


def use_data_offset_gui(self):

    use_data_offset_gui = [self.use_offset_Data1,self.use_offset_Data2,self.use_offset_Data3,self.use_offset_Data4, self.use_offset_Data5,
                           self.use_offset_Data6,self.use_offset_Data7,self.use_offset_Data8,self.use_offset_Data9,self.use_offset_Data10,
                           self.use_offset_Data11,self.use_offset_Data12,self.use_offset_Data13,self.use_offset_Data14, self.use_offset_Data15,
                           self.use_offset_Data16,self.use_offset_Data17,self.use_offset_Data18,self.use_offset_Data19,self.use_offset_Data20
                           ] 
    return use_data_offset_gui


def use_data_jitter_gui(self):

    use_data_jitter_gui = [self.use_jitter_Data1,self.use_jitter_Data2,self.use_jitter_Data3,self.use_jitter_Data4,self.use_jitter_Data5,
                           self.use_jitter_Data6,self.use_jitter_Data7,self.use_jitter_Data8,self.use_jitter_Data9,self.use_jitter_Data10,
                           self.use_jitter_Data11,self.use_jitter_Data12,self.use_jitter_Data13,self.use_jitter_Data14,self.use_jitter_Data15,
                           self.use_jitter_Data16,self.use_jitter_Data17,self.use_jitter_Data18,self.use_jitter_Data19,self.use_jitter_Data20
                           ]
    return use_data_jitter_gui


def data_errors_gui(self):

    data_errors_gui = [
            self.err_Data1,self.err_Data2,self.err_Data3,self.err_Data4,self.err_Data5,
            self.err_Data6,self.err_Data7,self.err_Data8,self.err_Data9,self.err_Data10,
            self.err_Data11,self.err_Data12,self.err_Data13,self.err_Data14,self.err_Data15,
            self.err_Data16,self.err_Data17,self.err_Data18,self.err_Data19,self.err_Data20
            ]
    return data_errors_gui


def data_errors_jitter_gui(self):

    data_errors_jitter_gui = [
            self.err_jitter_Data1,self.err_jitter_Data2,self.err_jitter_Data3,self.err_jitter_Data4,self.err_jitter_Data5,
            self.err_jitter_Data6,self.err_jitter_Data7,self.err_jitter_Data8,self.err_jitter_Data9,self.err_jitter_Data10,
            self.err_jitter_Data11,self.err_jitter_Data12,self.err_jitter_Data13,self.err_jitter_Data14,self.err_jitter_Data15,
            self.err_jitter_Data16,self.err_jitter_Data17,self.err_jitter_Data18,self.err_jitter_Data19,self.err_jitter_Data20
            ]

    return data_errors_jitter_gui


######################### Transit data ###########################

def tra_data_gui(self):

    tra_data_gui = [
            self.trans_Data_1,self.trans_Data_2,self.trans_Data_3,self.trans_Data_4,self.trans_Data_5,
            self.trans_Data_6,self.trans_Data_7,self.trans_Data_8,self.trans_Data_9,self.trans_Data_10,
            self.trans_Data_11,self.trans_Data_12,self.trans_Data_13,self.trans_Data_14,self.trans_Data_15,
            self.trans_Data_16,self.trans_Data_17,self.trans_Data_18,self.trans_Data_19,self.trans_Data_20
            ]
    return tra_data_gui

def use_tra_data_offset_gui(self):

    use_tra_data_offset_gui = [
            self.use_offset_trans_Data_1,self.use_offset_trans_Data_2,self.use_offset_trans_Data_3,self.use_offset_trans_Data_4,
            self.use_offset_trans_Data_5,self.use_offset_trans_Data_6,self.use_offset_trans_Data_7,self.use_offset_trans_Data_8,
            self.use_offset_trans_Data_9,self.use_offset_trans_Data_10,
            self.use_offset_trans_Data_11,self.use_offset_trans_Data_12,self.use_offset_trans_Data_13,self.use_offset_trans_Data_14,
            self.use_offset_trans_Data_15,self.use_offset_trans_Data_16,self.use_offset_trans_Data_17,self.use_offset_trans_Data_18,
            self.use_offset_trans_Data_19,self.use_offset_trans_Data_20
            ]

    return use_tra_data_offset_gui

def tra_data_errors_gui(self):

    tra_data_errors_gui = [
            self.err_trans_Data_1,self.err_trans_Data_2,self.err_trans_Data_3,self.err_trans_Data_4,self.err_trans_Data_5,
            self.err_trans_Data_6,self.err_trans_Data_7,self.err_trans_Data_8,self.err_trans_Data_9,self.err_trans_Data_10,
            self.err_trans_Data_11,self.err_trans_Data_12,self.err_trans_Data_13,self.err_trans_Data_14,self.err_trans_Data_15,
            self.err_trans_Data_16,self.err_trans_Data_17,self.err_trans_Data_18,self.err_trans_Data_19,self.err_trans_Data_20
            ]
    
    return tra_data_errors_gui

           
def tra_data_jitter_gui(self):

    tra_data_jitter_gui = [
            self.jitter_trans_Data_1,self.jitter_trans_Data_2,self.jitter_trans_Data_3,self.jitter_trans_Data_4,self.jitter_trans_Data_5,
            self.jitter_trans_Data_6,self.jitter_trans_Data_7,self.jitter_trans_Data_8,self.jitter_trans_Data_9,self.jitter_trans_Data_10,
            self.jitter_trans_Data_11,self.jitter_trans_Data_12,self.jitter_trans_Data_13,self.jitter_trans_Data_14,self.jitter_trans_Data_15,
            self.jitter_trans_Data_16,self.jitter_trans_Data_17,self.jitter_trans_Data_18,self.jitter_trans_Data_19,self.jitter_trans_Data_20

            ]
    return tra_data_jitter_gui
 
    
def use_tra_data_jitter_gui(self):
    
    use_tra_data_jitter_gui = [
            self.use_jitter_trans_Data_1,self.use_jitter_trans_Data_2,self.use_jitter_trans_Data_3,self.use_jitter_trans_Data_4,
            self.use_jitter_trans_Data_5,self.use_jitter_trans_Data_6,self.use_jitter_trans_Data_7,self.use_jitter_trans_Data_8,
            self.use_jitter_trans_Data_9,self.use_jitter_trans_Data_10,
            self.use_jitter_trans_Data_11,self.use_jitter_trans_Data_12,self.use_jitter_trans_Data_13,self.use_jitter_trans_Data_14,
            self.use_jitter_trans_Data_15,self.use_jitter_trans_Data_16,self.use_jitter_trans_Data_17,self.use_jitter_trans_Data_18,
            self.use_jitter_trans_Data_19,self.use_jitter_trans_Data_20
            ]
    
    return use_tra_data_jitter_gui


def tra_data_errors_jitter_gui(self):

    tra_data_errors_jitter_gui = [
            self.err_jitter_trans_Data_1,self.err_jitter_trans_Data_2,self.err_jitter_trans_Data_3,self.err_jitter_trans_Data_4,
            self.err_jitter_trans_Data_5,self.err_jitter_trans_Data_6,self.err_jitter_trans_Data_7,self.err_jitter_trans_Data_8,
            self.err_jitter_trans_Data_9,self.err_jitter_trans_Data_10,
            self.err_jitter_trans_Data_11,self.err_jitter_trans_Data_12,self.err_jitter_trans_Data_13,self.err_jitter_trans_Data_14,
            self.err_jitter_trans_Data_15,self.err_jitter_trans_Data_16,self.err_jitter_trans_Data_17,self.err_jitter_trans_Data_18,
            self.err_jitter_trans_Data_19,self.err_jitter_trans_Data_20
            ]
    
    return tra_data_errors_jitter_gui



def tra_data_lin_trend_gui(self):

    tra_data_lin_trend_gui = [
            self.tra_lin_trend_1,self.tra_lin_trend_2,self.tra_lin_trend_3,self.tra_lin_trend_4,self.tra_lin_trend_5,
            self.tra_lin_trend_6,self.tra_lin_trend_7,self.tra_lin_trend_8,self.tra_lin_trend_9,self.tra_lin_trend_10,
            self.tra_lin_trend_11,self.tra_lin_trend_12,self.tra_lin_trend_13,self.tra_lin_trend_14,self.tra_lin_trend_15,
            self.tra_lin_trend_16,self.tra_lin_trend_17,self.tra_lin_trend_18,self.tra_lin_trend_19,self.tra_lin_trend_20
            ]
    return tra_data_lin_trend_gui

def use_tra_data_lin_trend_gui(self):

    use_tra_data_lin_trend_gui = [
            self.use_tra_lin_trend_1,self.use_tra_lin_trend_2,self.use_tra_lin_trend_3,self.use_tra_lin_trend_4,
            self.use_tra_lin_trend_5,self.use_tra_lin_trend_6,self.use_tra_lin_trend_7,self.use_tra_lin_trend_8,
            self.use_tra_lin_trend_9,self.use_tra_lin_trend_10,
            self.use_tra_lin_trend_11,self.use_tra_lin_trend_12,self.use_tra_lin_trend_13,self.use_tra_lin_trend_14,
            self.use_tra_lin_trend_15,self.use_tra_lin_trend_16,self.use_tra_lin_trend_17,self.use_tra_lin_trend_18,
            self.use_tra_lin_trend_19,self.use_tra_lin_trend_20
            ]

    return use_tra_data_lin_trend_gui

def err_tra_data_lin_trend_gui(self):

    err_tra_data_lin_trend_gui = [
            self.err_tra_lin_trend_1,self.err_tra_lin_trend_2,self.err_tra_lin_trend_3,self.err_tra_lin_trend_4,self.err_tra_lin_trend_5,
            self.err_tra_lin_trend_6,self.err_tra_lin_trend_7,self.err_tra_lin_trend_8,self.err_tra_lin_trend_9,self.err_tra_lin_trend_10,
            self.err_tra_lin_trend_11,self.err_tra_lin_trend_12,self.err_tra_lin_trend_13,self.err_tra_lin_trend_14,self.err_tra_lin_trend_15,
            self.err_tra_lin_trend_16,self.err_tra_lin_trend_17,self.err_tra_lin_trend_18,self.err_tra_lin_trend_19,self.err_tra_lin_trend_20
            ]
    
    return err_tra_data_lin_trend_gui


def tra_data_quad_trend_gui(self):

    tra_data_quad_trend_gui = [
            self.tra_quad_trend_1,self.tra_quad_trend_2,self.tra_quad_trend_3,self.tra_quad_trend_4,self.tra_quad_trend_5,
            self.tra_quad_trend_6,self.tra_quad_trend_7,self.tra_quad_trend_8,self.tra_quad_trend_9,self.tra_quad_trend_10,
            self.tra_quad_trend_11,self.tra_quad_trend_12,self.tra_quad_trend_13,self.tra_quad_trend_14,self.tra_quad_trend_15,
            self.tra_quad_trend_16,self.tra_quad_trend_17,self.tra_quad_trend_18,self.tra_quad_trend_19,self.tra_quad_trend_20
            ]
    return tra_data_quad_trend_gui

def use_tra_data_quad_trend_gui(self):

    use_tra_data_quad_trend_gui = [
            self.use_tra_quad_trend_1,self.use_tra_quad_trend_2,self.use_tra_quad_trend_3,self.use_tra_quad_trend_4,
            self.use_tra_quad_trend_5,self.use_tra_quad_trend_6,self.use_tra_quad_trend_7,self.use_tra_quad_trend_8,
            self.use_tra_quad_trend_9,self.use_tra_quad_trend_10,
            self.use_tra_quad_trend_11,self.use_tra_quad_trend_12,self.use_tra_quad_trend_13,self.use_tra_quad_trend_14,
            self.use_tra_quad_trend_15,self.use_tra_quad_trend_16,self.use_tra_quad_trend_17,self.use_tra_quad_trend_18,
            self.use_tra_quad_trend_19,self.use_tra_quad_trend_20
            ]

    return use_tra_data_quad_trend_gui

def err_tra_data_quad_trend_gui(self):

    err_tra_data_quad_trend_gui = [
            self.err_tra_quad_trend_1,self.err_tra_quad_trend_2,self.err_tra_quad_trend_3,self.err_tra_quad_trend_4,self.err_tra_quad_trend_5,
            self.err_tra_quad_trend_6,self.err_tra_quad_trend_7,self.err_tra_quad_trend_8,self.err_tra_quad_trend_9,self.err_tra_quad_trend_10,
            self.err_tra_quad_trend_11,self.err_tra_quad_trend_12,self.err_tra_quad_trend_13,self.err_tra_quad_trend_14,self.err_tra_quad_trend_15,
            self.err_tra_quad_trend_16,self.err_tra_quad_trend_17,self.err_tra_quad_trend_18,self.err_tra_quad_trend_19,self.err_tra_quad_trend_20
            ]
    
    return err_tra_data_quad_trend_gui



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
    [self.Data6_min,self.Data6_max], [self.Data7_min,self.Data7_max], [self.Data8_min,self.Data8_max], [self.Data9_min,self.Data9_max], [self.Data10_min,self.Data10_max],
    [self.Data11_min,self.Data11_max], [self.Data12_min,self.Data12_max], [self.Data13_min,self.Data13_max], [self.Data14_min,self.Data14_max], [self.Data15_min,self.Data15_max],   
    [self.Data16_min,self.Data16_max], [self.Data17_min,self.Data17_max], [self.Data18_min,self.Data18_max], [self.Data19_min,self.Data19_max], [self.Data20_min,self.Data20_max]
    ]
    
    return offset_bounds_gui
        


def jitter_bounds_gui(self):  
      
    jitter_bounds_gui = [
    [self.jitter1_min,self.jitter1_max], [self.jitter2_min,self.jitter2_max], [self.jitter3_min,self.jitter3_max], [self.jitter4_min,self.jitter4_max], [self.jitter5_min,self.jitter5_max],   
    [self.jitter6_min,self.jitter6_max], [self.jitter7_min,self.jitter7_max], [self.jitter8_min,self.jitter8_max], [self.jitter9_min,self.jitter9_max], [self.jitter10_min,self.Data10_max],
    [self.jitter11_min,self.jitter11_max], [self.jitter12_min,self.jitter12_max], [self.jitter13_min,self.jitter13_max], [self.jitter14_min,self.jitter14_max], [self.jitter15_min,self.jitter15_max],   
    [self.jitter16_min,self.jitter16_max], [self.jitter17_min,self.jitter17_max], [self.jitter18_min,self.jitter18_max], [self.jitter19_min,self.jitter19_max], [self.jitter20_min,self.Data20_max]
    ]  
    
    return jitter_bounds_gui




def offset_bounds_gui_tra(self):
    
    offset_bounds_gui_tra = [
    [self.tra_Data_min_1,self.tra_Data_max_1], [self.tra_Data_min_2,self.tra_Data_max_2], [self.tra_Data_min_3,self.tra_Data_max_3],
    [self.tra_Data_min_4,self.tra_Data_max_4], [self.tra_Data_min_5,self.tra_Data_max_5], [self.tra_Data_min_6,self.tra_Data_max_6],
    [self.tra_Data_min_7,self.tra_Data_max_7], [self.tra_Data_min_8,self.tra_Data_max_8], [self.tra_Data_min_9,self.tra_Data_max_9],
    [self.tra_Data_min_10,self.tra_Data_max_10],
    [self.tra_Data_min_11,self.tra_Data_max_11], [self.tra_Data_min_12,self.tra_Data_max_12], [self.tra_Data_min_13,self.tra_Data_max_13],
    [self.tra_Data_min_14,self.tra_Data_max_14], [self.tra_Data_min_15,self.tra_Data_max_15], [self.tra_Data_min_16,self.tra_Data_max_16],
    [self.tra_Data_min_17,self.tra_Data_max_17], [self.tra_Data_min_18,self.tra_Data_max_18], [self.tra_Data_min_19,self.tra_Data_max_19],
    [self.tra_Data_min_20,self.tra_Data_max_20]
    ]

    return offset_bounds_gui_tra



def jitter_bounds_gui_tra(self):  
      
    jitter_bounds_gui_tra = [
            [self.tra_jitter_min_1,self.tra_jitter_max_1],[self.tra_jitter_min_2,self.tra_jitter_max_2],
            [self.tra_jitter_min_3,self.tra_jitter_max_3],[self.tra_jitter_min_4,self.tra_jitter_max_4],
            [self.tra_jitter_min_5,self.tra_jitter_max_5],[self.tra_jitter_min_6,self.tra_jitter_max_6],
            [self.tra_jitter_min_7,self.tra_jitter_max_7],[self.tra_jitter_min_8,self.tra_jitter_max_8],
            [self.tra_jitter_min_9,self.tra_jitter_max_9],[self.tra_jitter_min_10,self.tra_jitter_max_10],
            [self.tra_jitter_min_11,self.tra_jitter_max_11],[self.tra_jitter_min_12,self.tra_jitter_max_12],
            [self.tra_jitter_min_13,self.tra_jitter_max_13],[self.tra_jitter_min_14,self.tra_jitter_max_14],
            [self.tra_jitter_min_15,self.tra_jitter_max_15],[self.tra_jitter_min_16,self.tra_jitter_max_16],
            [self.tra_jitter_min_17,self.tra_jitter_max_17],[self.tra_jitter_min_18,self.tra_jitter_max_18],
            [self.tra_jitter_min_19,self.tra_jitter_max_19],[self.tra_jitter_min_20,self.tra_jitter_max_20]
    ]

    return jitter_bounds_gui_tra



################### h = e sin(w) ########################

def h_bounds_gui(self):

    h_bounds_gui = [
    [self.h_min_1,self.h_max_1], [self.h_min_2,self.h_max_2], 
    [self.h_min_3,self.h_max_3], [self.h_min_4,self.h_max_4], 
    [self.h_min_5,self.h_max_5], [self.h_min_6,self.h_max_6], 
    [self.h_min_7,self.h_max_7], [self.h_min_8,self.h_max_8], 
    [self.h_min_9,self.h_max_9] 
    ]  
    return h_bounds_gui


def h_norm_pr_gui(self):

    h_norm_pr_gui = [
    [self.h_mean_1,self.h_sigma_1,self.use_h_norm_pr_1], 
    [self.h_mean_2,self.h_sigma_2,self.use_h_norm_pr_2], 
    [self.h_mean_3,self.h_sigma_3,self.use_h_norm_pr_3], 
    [self.h_mean_4,self.h_sigma_4,self.use_h_norm_pr_4], 
    [self.h_mean_5,self.h_sigma_5,self.use_h_norm_pr_5], 
    [self.h_mean_6,self.h_sigma_6,self.use_h_norm_pr_6], 
    [self.h_mean_7,self.h_sigma_7,self.use_h_norm_pr_7], 
    [self.h_mean_8,self.h_sigma_8,self.use_h_norm_pr_8], 
    [self.h_mean_9,self.h_sigma_9,self.use_h_norm_pr_9]
    ]

    return h_norm_pr_gui

def h_jeff_pr_gui(self):

    h_jeff_pr_gui = [
    [self.h_alpha_1,self.h_beta_1,self.use_h_jeff_pr_1], 
    [self.h_alpha_2,self.h_beta_2,self.use_h_jeff_pr_2], 
    [self.h_alpha_3,self.h_beta_3,self.use_h_jeff_pr_3], 
    [self.h_alpha_4,self.h_beta_4,self.use_h_jeff_pr_4], 
    [self.h_alpha_5,self.h_beta_5,self.use_h_jeff_pr_5], 
    [self.h_alpha_6,self.h_beta_6,self.use_h_jeff_pr_6], 
    [self.h_alpha_7,self.h_beta_7,self.use_h_jeff_pr_7], 
    [self.h_alpha_8,self.h_beta_8,self.use_h_jeff_pr_8], 
    [self.h_alpha_9,self.h_beta_9,self.use_h_jeff_pr_9]
    ]

    return h_jeff_pr_gui



################### k = e cos(w) ########################

def k_bounds_gui(self):

    k_bounds_gui = [
    [self.k_min_1,self.k_max_1], [self.k_min_2,self.k_max_2], 
    [self.k_min_3,self.k_max_3], [self.k_min_4,self.k_max_4], 
    [self.k_min_5,self.k_max_5], [self.k_min_6,self.k_max_6], 
    [self.k_min_7,self.k_max_7], [self.k_min_8,self.k_max_8], 
    [self.k_min_9,self.k_max_9] 
    ]  
    return k_bounds_gui


def k_norm_pr_gui(self):

    k_norm_pr_gui = [
    [self.k_mean_1,self.k_sigma_1,self.use_k_norm_pr_1], 
    [self.k_mean_2,self.k_sigma_2,self.use_k_norm_pr_2], 
    [self.k_mean_3,self.k_sigma_3,self.use_k_norm_pr_3], 
    [self.k_mean_4,self.k_sigma_4,self.use_k_norm_pr_4], 
    [self.k_mean_5,self.k_sigma_5,self.use_k_norm_pr_5], 
    [self.k_mean_6,self.k_sigma_6,self.use_k_norm_pr_6], 
    [self.k_mean_7,self.k_sigma_7,self.use_k_norm_pr_7], 
    [self.k_mean_8,self.k_sigma_8,self.use_k_norm_pr_8], 
    [self.k_mean_9,self.k_sigma_9,self.use_k_norm_pr_9]
    ]

    return k_norm_pr_gui

def k_jeff_pr_gui(self):

    k_jeff_pr_gui = [
    [self.k_alpha_1,self.k_beta_1,self.use_k_jeff_pr_1], 
    [self.k_alpha_2,self.k_beta_2,self.use_k_jeff_pr_2], 
    [self.k_alpha_3,self.k_beta_3,self.use_k_jeff_pr_3], 
    [self.k_alpha_4,self.k_beta_4,self.use_k_jeff_pr_4], 
    [self.k_alpha_5,self.k_beta_5,self.use_k_jeff_pr_5], 
    [self.k_alpha_6,self.k_beta_6,self.use_k_jeff_pr_6], 
    [self.k_alpha_7,self.k_beta_7,self.use_k_jeff_pr_7], 
    [self.k_alpha_8,self.k_beta_8,self.use_k_jeff_pr_8], 
    [self.k_alpha_9,self.k_beta_9,self.use_k_jeff_pr_9]
    ]

    return k_jeff_pr_gui



################### lambda = M0 + w ########################

def lambda_bounds_gui(self):

    lambda_bounds_gui = [
    [self.lambda_min_1,self.lambda_max_1], [self.lambda_min_2,self.lambda_max_2], 
    [self.lambda_min_3,self.lambda_max_3], [self.lambda_min_4,self.lambda_max_4], 
    [self.lambda_min_5,self.lambda_max_5], [self.lambda_min_6,self.lambda_max_6], 
    [self.lambda_min_7,self.lambda_max_7], [self.lambda_min_8,self.lambda_max_8], 
    [self.lambda_min_9,self.lambda_max_9] 
    ]  
    return lambda_bounds_gui


def lambda_norm_pr_gui(self):

    lambda_norm_pr_gui = [
    [self.lambda_mean_1,self.lambda_sigma_1,self.use_lambda_norm_pr_1], 
    [self.lambda_mean_2,self.lambda_sigma_2,self.use_lambda_norm_pr_2], 
    [self.lambda_mean_3,self.lambda_sigma_3,self.use_lambda_norm_pr_3], 
    [self.lambda_mean_4,self.lambda_sigma_4,self.use_lambda_norm_pr_4], 
    [self.lambda_mean_5,self.lambda_sigma_5,self.use_lambda_norm_pr_5], 
    [self.lambda_mean_6,self.lambda_sigma_6,self.use_lambda_norm_pr_6], 
    [self.lambda_mean_7,self.lambda_sigma_7,self.use_lambda_norm_pr_7], 
    [self.lambda_mean_8,self.lambda_sigma_8,self.use_lambda_norm_pr_8], 
    [self.lambda_mean_9,self.lambda_sigma_9,self.use_lambda_norm_pr_9]
    ]

    return lambda_norm_pr_gui

def lambda_jeff_pr_gui(self):

    lambda_jeff_pr_gui = [
    [self.lambda_alpha_1,self.lambda_beta_1,self.use_lambda_jeff_pr_1], 
    [self.lambda_alpha_2,self.lambda_beta_2,self.use_lambda_jeff_pr_2], 
    [self.lambda_alpha_3,self.lambda_beta_3,self.use_lambda_jeff_pr_3], 
    [self.lambda_alpha_4,self.lambda_beta_4,self.use_lambda_jeff_pr_4], 
    [self.lambda_alpha_5,self.lambda_beta_5,self.use_lambda_jeff_pr_5], 
    [self.lambda_alpha_6,self.lambda_beta_6,self.use_lambda_jeff_pr_6], 
    [self.lambda_alpha_7,self.lambda_beta_7,self.use_lambda_jeff_pr_7], 
    [self.lambda_alpha_8,self.lambda_beta_8,self.use_lambda_jeff_pr_8], 
    [self.lambda_alpha_9,self.lambda_beta_9,self.use_lambda_jeff_pr_9]
    ]

    return lambda_jeff_pr_gui



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



def om_dot_bounds_gui2(self):

    om_dot_bounds_gui = [
    [self.omega_dot_min_1,self.omega_dot_max_1,self.use_omega_dot_bound_1], 
    [self.omega_dot_min_2,self.omega_dot_max_2,self.use_omega_dot_bound_2], 
    [self.omega_dot_min_3,self.omega_dot_max_3,self.use_omega_dot_bound_3], 
    [self.omega_dot_min_4,self.omega_dot_max_4,self.use_omega_dot_bound_4], 
    [self.omega_dot_min_5,self.omega_dot_max_5,self.use_omega_dot_bound_5], 
    [self.omega_dot_min_6,self.omega_dot_max_6,self.use_omega_dot_bound_6], 
    [self.omega_dot_min_7,self.omega_dot_max_7,self.use_omega_dot_bound_7], 
    [self.omega_dot_min_8,self.omega_dot_max_8,self.use_omega_dot_bound_8], 
    [self.omega_dot_min_9,self.omega_dot_max_9,self.use_omega_dot_bound_9]
    ]

    return om_dot_bounds_gui



def om_dot_norm_pr_gui(self):

    om_dot_norm_pr_gui = [
    [self.omega_dot_mean_1,self.omega_dot_sigma_1,self.use_omega_dot_norm_pr_1], 
    [self.omega_dot_mean_2,self.omega_dot_sigma_2,self.use_omega_dot_norm_pr_2], 
    [self.omega_dot_mean_3,self.omega_dot_sigma_3,self.use_omega_dot_norm_pr_3], 
    [self.omega_dot_mean_4,self.omega_dot_sigma_4,self.use_omega_dot_norm_pr_4], 
    [self.omega_dot_mean_5,self.omega_dot_sigma_5,self.use_omega_dot_norm_pr_5], 
    [self.omega_dot_mean_6,self.omega_dot_sigma_6,self.use_omega_dot_norm_pr_6], 
    [self.omega_dot_mean_7,self.omega_dot_sigma_7,self.use_omega_dot_norm_pr_7], 
    [self.omega_dot_mean_8,self.omega_dot_sigma_8,self.use_omega_dot_norm_pr_8], 
    [self.omega_dot_mean_9,self.omega_dot_sigma_9,self.use_omega_dot_norm_pr_9]
    ]

    return om_dot_norm_pr_gui



def om_dot_jeff_pr_gui(self):

    om_dot_jeff_pr_gui = [
    [self.omega_dot_alpha_1,self.omega_dot_beta_1,self.use_omega_dot_jeff_pr_1], 
    [self.omega_dot_alpha_2,self.omega_dot_beta_2,self.use_omega_dot_jeff_pr_2], 
    [self.omega_dot_alpha_3,self.omega_dot_beta_3,self.use_omega_dot_jeff_pr_3], 
    [self.omega_dot_alpha_4,self.omega_dot_beta_4,self.use_omega_dot_jeff_pr_4], 
    [self.omega_dot_alpha_5,self.omega_dot_beta_5,self.use_omega_dot_jeff_pr_5], 
    [self.omega_dot_alpha_6,self.omega_dot_beta_6,self.use_omega_dot_jeff_pr_6], 
    [self.omega_dot_alpha_7,self.omega_dot_beta_7,self.use_omega_dot_jeff_pr_7], 
    [self.omega_dot_alpha_8,self.omega_dot_beta_8,self.use_omega_dot_jeff_pr_8], 
    [self.omega_dot_alpha_9,self.omega_dot_beta_9,self.use_omega_dot_jeff_pr_9]
    ]

    return om_dot_jeff_pr_gui


################### Tra. reg.  ########################
 
def data_tra_reg_group(self):

    data_tra_reg_group = [
            [self.tra_reg_bjd_1,self.tra_reg_airmass_1],[self.tra_reg_bjd_2,self.tra_reg_airmass_2],
            [self.tra_reg_bjd_3,self.tra_reg_airmass_3],[self.tra_reg_bjd_4,self.tra_reg_airmass_4],
            [self.tra_reg_bjd_5,self.tra_reg_airmass_5],[self.tra_reg_bjd_6,self.tra_reg_airmass_6],
            [self.tra_reg_bjd_7,self.tra_reg_airmass_7],[self.tra_reg_bjd_8,self.tra_reg_airmass_8],
            [self.tra_reg_bjd_9,self.tra_reg_airmass_9],[self.tra_reg_bjd_10,self.tra_reg_airmass_10],
            [self.tra_reg_bjd_11,self.tra_reg_airmass_11],[self.tra_reg_bjd_12,self.tra_reg_airmass_12],
            [self.tra_reg_bjd_13,self.tra_reg_airmass_13],[self.tra_reg_bjd_14,self.tra_reg_airmass_14],
            [self.tra_reg_bjd_15,self.tra_reg_airmass_15],[self.tra_reg_bjd_16,self.tra_reg_airmass_16],
            [self.tra_reg_bjd_17,self.tra_reg_airmass_17],[self.tra_reg_bjd_18,self.tra_reg_airmass_18],
            [self.tra_reg_bjd_19,self.tra_reg_airmass_19],[self.tra_reg_bjd_20,self.tra_reg_airmass_20]
            ]
    return data_tra_reg_group


################### LD  ########################
 
def data_ld_group(self):

    data_ld_group = [
            self.LD_group_1,self.LD_group_2,self.LD_group_3,self.LD_group_4,self.LD_group_5,
            self.LD_group_6,self.LD_group_7,self.LD_group_8,self.LD_group_9,self.LD_group_10,
            self.LD_group_11,self.LD_group_12,self.LD_group_13,self.LD_group_14,self.LD_group_15,
            self.LD_group_16,self.LD_group_17,self.LD_group_18,self.LD_group_19,self.LD_group_20
            ]
    return data_ld_group


def use_uni_ld_models(self):

    use_uni_ld_models = [
            self.use_uniform_ld_1,self.use_uniform_ld_2,self.use_uniform_ld_3,self.use_uniform_ld_4,self.use_uniform_ld_5,
            self.use_uniform_ld_6,self.use_uniform_ld_7,self.use_uniform_ld_8,self.use_uniform_ld_9,self.use_uniform_ld_10,
            self.use_uniform_ld_11,self.use_uniform_ld_12,self.use_uniform_ld_13,self.use_uniform_ld_14,self.use_uniform_ld_15,
            self.use_uniform_ld_16,self.use_uniform_ld_17,self.use_uniform_ld_18,self.use_uniform_ld_19,self.use_uniform_ld_20

            ]
    return use_uni_ld_models


def use_lin_ld_models(self):

    use_lin_ld_models = [
            self.use_linear_ld_1,self.use_linear_ld_2,self.use_linear_ld_3,self.use_linear_ld_4,self.use_linear_ld_5,
            self.use_linear_ld_6,self.use_linear_ld_7,self.use_linear_ld_8,self.use_linear_ld_9,self.use_linear_ld_10,
            self.use_linear_ld_11,self.use_linear_ld_12,self.use_linear_ld_13,self.use_linear_ld_14,self.use_linear_ld_15,
            self.use_linear_ld_16,self.use_linear_ld_17,self.use_linear_ld_18,self.use_linear_ld_19,self.use_linear_ld_20
            ]
    return use_lin_ld_models


def use_quad_ld_models(self):
    
    use_quad_ld_models =[
            self.use_quadratic_ld_1,self.use_quadratic_ld_2,self.use_quadratic_ld_3,self.use_quadratic_ld_4,self.use_quadratic_ld_5,
            self.use_quadratic_ld_6,self.use_quadratic_ld_7,self.use_quadratic_ld_8,self.use_quadratic_ld_9,self.use_quadratic_ld_10,
            self.use_quadratic_ld_11,self.use_quadratic_ld_12,self.use_quadratic_ld_13,self.use_quadratic_ld_14,self.use_quadratic_ld_15,
            self.use_quadratic_ld_16,self.use_quadratic_ld_17,self.use_quadratic_ld_18,self.use_quadratic_ld_19,self.use_quadratic_ld_20
            ] 
    return use_quad_ld_models


def use_nonlin_ld_models(self):
    
    use_nonlin_ld_models = [
            self.use_nonlinear_ld_1,self.use_nonlinear_ld_2,self.use_nonlinear_ld_3,self.use_nonlinear_ld_4,self.use_nonlinear_ld_5,
            self.use_nonlinear_ld_6,self.use_nonlinear_ld_7,self.use_nonlinear_ld_8,self.use_nonlinear_ld_9,self.use_nonlinear_ld_10,
            self.use_nonlinear_ld_11,self.use_nonlinear_ld_12,self.use_nonlinear_ld_13,self.use_nonlinear_ld_14,self.use_nonlinear_ld_15,
            self.use_nonlinear_ld_16,self.use_nonlinear_ld_17,self.use_nonlinear_ld_18,self.use_nonlinear_ld_19,self.use_nonlinear_ld_20
            ] 
    return use_nonlin_ld_models


def lin_u(self):

    lin_u = [self.u1_linear_1,self.u1_linear_2,self.u1_linear_3,self.u1_linear_4,self.u1_linear_5,
             self.u1_linear_6,self.u1_linear_7,self.u1_linear_8,self.u1_linear_9,self.u1_linear_10,
             self.u1_linear_11,self.u1_linear_12,self.u1_linear_13,self.u1_linear_14,self.u1_linear_15,
             self.u1_linear_16,self.u1_linear_17,self.u1_linear_18,self.u1_linear_19,self.u1_linear_20

             ]
    return lin_u

def use_lin_u(self):

    use_lin_u = [
            self.use_u1_linear_1,self.use_u1_linear_2,self.use_u1_linear_3,self.use_u1_linear_4,self.use_u1_linear_5,
            self.use_u1_linear_6,self.use_u1_linear_7,self.use_u1_linear_8,self.use_u1_linear_9,self.use_u1_linear_10,
            self.use_u1_linear_11,self.use_u1_linear_12,self.use_u1_linear_13,self.use_u1_linear_14,self.use_u1_linear_15,
            self.use_u1_linear_16,self.use_u1_linear_17,self.use_u1_linear_18,self.use_u1_linear_19,self.use_u1_linear_20
            ]
    return use_lin_u

def err_lin_u(self):

    err_lin_u = [self.err_u1_linear_1,self.err_u1_linear_2,self.err_u1_linear_3,self.err_u1_linear_4,self.err_u1_linear_5,
                 self.err_u1_linear_6,self.err_u1_linear_7,self.err_u1_linear_8,self.err_u1_linear_9,self.err_u1_linear_10,
                 self.err_u1_linear_11,self.err_u1_linear_12,self.err_u1_linear_13,self.err_u1_linear_14,self.err_u1_linear_15,
                 self.err_u1_linear_16,self.err_u1_linear_17,self.err_u1_linear_18,self.err_u1_linear_19,self.err_u1_linear_20
             ]
    return err_lin_u


def quad_u1(self):

    quad_u1 = [
            self.u1_quadratic_1,self.u1_quadratic_2,self.u1_quadratic_3,self.u1_quadratic_4,self.u1_quadratic_5,
            self.u1_quadratic_6,self.u1_quadratic_7,self.u1_quadratic_8,self.u1_quadratic_9,self.u1_quadratic_10,
            self.u1_quadratic_11,self.u1_quadratic_12,self.u1_quadratic_13,self.u1_quadratic_14,self.u1_quadratic_15,
            self.u1_quadratic_16,self.u1_quadratic_17,self.u1_quadratic_18,self.u1_quadratic_19,self.u1_quadratic_20
             ]
    return quad_u1

def use_quad_u1(self):

    use_quad_u1 =  [
            self.use_u1_quadratic_1,self.use_u1_quadratic_2,self.use_u1_quadratic_3,self.use_u1_quadratic_4,self.use_u1_quadratic_5,
            self.use_u1_quadratic_6,self.use_u1_quadratic_7,self.use_u1_quadratic_8,self.use_u1_quadratic_9,self.use_u1_quadratic_10,
            self.use_u1_quadratic_11,self.use_u1_quadratic_12,self.use_u1_quadratic_13,self.use_u1_quadratic_14,self.use_u1_quadratic_15,
            self.use_u1_quadratic_16,self.use_u1_quadratic_17,self.use_u1_quadratic_18,self.use_u1_quadratic_19,self.use_u1_quadratic_20
             ]
    return use_quad_u1

def err_quad_u1(self):

    err_quad_u1 =  [
            self.err_u1_quadratic_1,self.err_u1_quadratic_2,self.err_u1_quadratic_3,self.err_u1_quadratic_4,self.err_u1_quadratic_5,
            self.err_u1_quadratic_6,self.err_u1_quadratic_7,self.err_u1_quadratic_8,self.err_u1_quadratic_9,self.err_u1_quadratic_10,
            self.err_u1_quadratic_11,self.err_u1_quadratic_12,self.err_u1_quadratic_13,self.err_u1_quadratic_14,self.err_u1_quadratic_15,
            self.err_u1_quadratic_16,self.err_u1_quadratic_17,self.err_u1_quadratic_18,self.err_u1_quadratic_19,self.err_u1_quadratic_20

             ]
    return err_quad_u1


def quad_u2(self):

    quad_u2 = [
            self.u2_quadratic_1,self.u2_quadratic_2,self.u2_quadratic_3,self.u2_quadratic_4,self.u2_quadratic_5,
            self.u2_quadratic_6,self.u2_quadratic_7,self.u2_quadratic_8,self.u2_quadratic_9,self.u2_quadratic_10,
            self.u2_quadratic_11,self.u2_quadratic_12,self.u2_quadratic_13,self.u2_quadratic_14,self.u2_quadratic_15,
            self.u2_quadratic_16,self.u2_quadratic_17,self.u2_quadratic_18,self.u2_quadratic_19,self.u2_quadratic_20
             ]
    return quad_u2

def use_quad_u2(self):

    use_quad_u2 =  [
            self.use_u2_quadratic_1,self.use_u2_quadratic_2,self.use_u2_quadratic_3,self.use_u2_quadratic_4,self.use_u2_quadratic_5,
            self.use_u2_quadratic_6,self.use_u2_quadratic_7,self.use_u2_quadratic_8,self.use_u2_quadratic_9,self.use_u2_quadratic_10,
            self.use_u2_quadratic_11,self.use_u2_quadratic_12,self.use_u2_quadratic_13,self.use_u2_quadratic_14,self.use_u2_quadratic_15,
            self.use_u2_quadratic_16,self.use_u2_quadratic_17,self.use_u2_quadratic_18,self.use_u2_quadratic_19,self.use_u2_quadratic_20

             ]
    return use_quad_u2

def err_quad_u2(self):

    err_quad_u2 =  [
            self.err_u2_quadratic_1,self.err_u2_quadratic_2,self.err_u2_quadratic_3,self.err_u2_quadratic_4,self.err_u2_quadratic_5,
            self.err_u2_quadratic_6,self.err_u2_quadratic_7,self.err_u2_quadratic_8,self.err_u2_quadratic_9,self.err_u2_quadratic_10,
            self.err_u2_quadratic_11,self.err_u2_quadratic_12,self.err_u2_quadratic_13,self.err_u2_quadratic_14,self.err_u2_quadratic_15,
            self.err_u2_quadratic_16,self.err_u2_quadratic_17,self.err_u2_quadratic_18,self.err_u2_quadratic_19,self.err_u2_quadratic_20
             ]
    return err_quad_u2
    


def nonlin_u1(self):

    nonlin_u1 = [
            self.u1_nonlin_1,self.u1_nonlin_2,self.u1_nonlin_3,self.u1_nonlin_4,self.u1_nonlin_5,
            self.u1_nonlin_6,self.u1_nonlin_7,self.u1_nonlin_8,self.u1_nonlin_9,self.u1_nonlin_10,
            self.u1_nonlin_11,self.u1_nonlin_12,self.u1_nonlin_13,self.u1_nonlin_14,self.u1_nonlin_15,
            self.u1_nonlin_16,self.u1_nonlin_17,self.u1_nonlin_18,self.u1_nonlin_19,self.u1_nonlin_20
             ]
    return nonlin_u1

def use_nonlin_u1(self):

    use_nonlin_u1 =  [
            self.use_u1_nonlin_1,self.use_u1_nonlin_2,self.use_u1_nonlin_3,self.use_u1_nonlin_4,self.use_u1_nonlin_5,
            self.use_u1_nonlin_6,self.use_u1_nonlin_7,self.use_u1_nonlin_8,self.use_u1_nonlin_9,self.use_u1_nonlin_10,
            self.use_u1_nonlin_11,self.use_u1_nonlin_12,self.use_u1_nonlin_13,self.use_u1_nonlin_14,self.use_u1_nonlin_15,
            self.use_u1_nonlin_16,self.use_u1_nonlin_17,self.use_u1_nonlin_18,self.use_u1_nonlin_19,self.use_u1_nonlin_20
             ]
    return use_nonlin_u1

def err_nonlin_u1(self):

    err_nonlin_u1 =  [
            self.err_u1_nonlin_1,self.err_u1_nonlin_2,self.err_u1_nonlin_3,self.err_u1_nonlin_4,self.err_u1_nonlin_5,
            self.err_u1_nonlin_6,self.err_u1_nonlin_7,self.err_u1_nonlin_8,self.err_u1_nonlin_9,self.err_u1_nonlin_10,
            self.err_u1_nonlin_11,self.err_u1_nonlin_12,self.err_u1_nonlin_13,self.err_u1_nonlin_14,self.err_u1_nonlin_15,
            self.err_u1_nonlin_16,self.err_u1_nonlin_17,self.err_u1_nonlin_18,self.err_u1_nonlin_19,self.err_u1_nonlin_20
             ]
    return err_nonlin_u1

def nonlin_u2(self):

    nonlin_u2 = [
            self.u2_nonlin_1,self.u2_nonlin_2,self.u2_nonlin_3,self.u2_nonlin_4,self.u2_nonlin_5,
            self.u2_nonlin_6,self.u2_nonlin_7,self.u2_nonlin_8,self.u2_nonlin_9,self.u2_nonlin_10,
            self.u2_nonlin_11,self.u2_nonlin_12,self.u2_nonlin_13,self.u2_nonlin_14,self.u2_nonlin_15,
            self.u2_nonlin_16,self.u2_nonlin_17,self.u2_nonlin_18,self.u2_nonlin_19,self.u2_nonlin_20
             ]
    return nonlin_u2

def use_nonlin_u2(self):

    use_nonlin_u2 =  [
            self.use_u2_nonlin_1,self.use_u2_nonlin_2,self.use_u2_nonlin_3,self.use_u2_nonlin_4,self.use_u2_nonlin_5,
            self.use_u2_nonlin_6,self.use_u2_nonlin_7,self.use_u2_nonlin_8,self.use_u2_nonlin_9,self.use_u2_nonlin_10,
            self.use_u2_nonlin_11,self.use_u2_nonlin_12,self.use_u2_nonlin_13,self.use_u2_nonlin_14,self.use_u2_nonlin_15,
            self.use_u2_nonlin_16,self.use_u2_nonlin_17,self.use_u2_nonlin_18,self.use_u2_nonlin_19,self.use_u2_nonlin_20
             ]
    return use_nonlin_u2

def err_nonlin_u2(self):

    err_nonlin_u2 =  [
            self.err_u2_nonlin_1,self.err_u2_nonlin_2,self.err_u2_nonlin_3,self.err_u2_nonlin_4,self.err_u2_nonlin_5,
            self.err_u2_nonlin_6,self.err_u2_nonlin_7,self.err_u2_nonlin_8,self.err_u2_nonlin_9,self.err_u2_nonlin_10,
            self.err_u2_nonlin_11,self.err_u2_nonlin_12,self.err_u2_nonlin_13,self.err_u2_nonlin_14,self.err_u2_nonlin_15,
            self.err_u2_nonlin_16,self.err_u2_nonlin_17,self.err_u2_nonlin_18,self.err_u2_nonlin_19,self.err_u2_nonlin_20
             ]
    return err_nonlin_u2

def nonlin_u3(self):

    nonlin_u3 = [
            self.u3_nonlin_1,self.u3_nonlin_2,self.u3_nonlin_3,self.u3_nonlin_4,self.u3_nonlin_5,
            self.u3_nonlin_6,self.u3_nonlin_7,self.u3_nonlin_8,self.u3_nonlin_9,self.u3_nonlin_10,
            self.u3_nonlin_11,self.u3_nonlin_12,self.u3_nonlin_13,self.u3_nonlin_14,self.u3_nonlin_15,
            self.u3_nonlin_16,self.u3_nonlin_17,self.u3_nonlin_18,self.u3_nonlin_19,self.u3_nonlin_20
             ]
    return nonlin_u3

def use_nonlin_u3(self):

    use_nonlin_u3 =  [
            self.use_u3_nonlin_1,self.use_u3_nonlin_2,self.use_u3_nonlin_3,self.use_u3_nonlin_4,self.use_u3_nonlin_5,
            self.use_u3_nonlin_6,self.use_u3_nonlin_7,self.use_u3_nonlin_8,self.use_u3_nonlin_9,self.use_u3_nonlin_10,
            self.use_u3_nonlin_11,self.use_u3_nonlin_12,self.use_u3_nonlin_13,self.use_u3_nonlin_14,self.use_u3_nonlin_15,
            self.use_u3_nonlin_16,self.use_u3_nonlin_17,self.use_u3_nonlin_18,self.use_u3_nonlin_19,self.use_u3_nonlin_20
             ]
    return use_nonlin_u3

def err_nonlin_u3(self):

    err_nonlin_u3 =  [
            self.err_u3_nonlin_1,self.err_u3_nonlin_2,self.err_u3_nonlin_3,self.err_u3_nonlin_4,self.err_u3_nonlin_5,
            self.err_u3_nonlin_6,self.err_u3_nonlin_7,self.err_u3_nonlin_8,self.err_u3_nonlin_9,self.err_u3_nonlin_10,
            self.err_u3_nonlin_11,self.err_u3_nonlin_12,self.err_u3_nonlin_13,self.err_u3_nonlin_14,self.err_u3_nonlin_15,
            self.err_u3_nonlin_16,self.err_u3_nonlin_17,self.err_u3_nonlin_18,self.err_u3_nonlin_19,self.err_u3_nonlin_20
             ]
    return err_nonlin_u3

def nonlin_u4(self):

    nonlin_u4 = [
            self.u4_nonlin_1,self.u4_nonlin_2,self.u4_nonlin_3,self.u4_nonlin_4,self.u4_nonlin_5,
            self.u4_nonlin_6,self.u4_nonlin_7,self.u4_nonlin_8,self.u4_nonlin_9,self.u4_nonlin_10,
            self.u4_nonlin_11,self.u4_nonlin_12,self.u4_nonlin_13,self.u4_nonlin_14,self.u4_nonlin_15,
            self.u4_nonlin_16,self.u4_nonlin_17,self.u4_nonlin_18,self.u4_nonlin_19,self.u4_nonlin_20
             ]
    return nonlin_u4

def use_nonlin_u4(self):

    use_nonlin_u4 =  [
            self.use_u4_nonlin_1,self.use_u4_nonlin_2,self.use_u4_nonlin_3,self.use_u4_nonlin_4,self.use_u4_nonlin_5,
            self.use_u4_nonlin_6,self.use_u4_nonlin_7,self.use_u4_nonlin_8,self.use_u4_nonlin_9,self.use_u4_nonlin_10,
            self.use_u4_nonlin_11,self.use_u4_nonlin_12,self.use_u4_nonlin_13,self.use_u4_nonlin_14,self.use_u4_nonlin_15,
            self.use_u4_nonlin_16,self.use_u4_nonlin_17,self.use_u4_nonlin_18,self.use_u4_nonlin_19,self.use_u4_nonlin_20
             ]
    return use_nonlin_u4

def err_nonlin_u4(self):

    err_nonlin_u4 =  [
            self.err_u4_nonlin_1,self.err_u4_nonlin_2,self.err_u4_nonlin_3,self.err_u4_nonlin_4,self.err_u4_nonlin_5,
            self.err_u4_nonlin_6,self.err_u4_nonlin_7,self.err_u4_nonlin_8,self.err_u4_nonlin_9,self.err_u4_nonlin_10,
            self.err_u4_nonlin_11,self.err_u4_nonlin_12,self.err_u4_nonlin_13,self.err_u4_nonlin_14,self.err_u4_nonlin_15,
            self.err_u4_nonlin_16,self.err_u4_nonlin_17,self.err_u4_nonlin_18,self.err_u4_nonlin_19,self.err_u4_nonlin_20
             ]
    return err_nonlin_u4

def ld_u1_bounds_gui(self):
    
        ld_u1_bounds_gui = [
        [self.u1_min_1,self.u1_max_1],[self.u1_min_2,self.u1_max_2],[self.u1_min_3,self.u1_max_3],
        [self.u1_min_4,self.u1_max_4],[self.u1_min_5,self.u1_max_5],[self.u1_min_6,self.u1_max_6],
        [self.u1_min_7,self.u1_max_7],[self.u1_min_8,self.u1_max_8],[self.u1_min_9,self.u1_max_9],
        [self.u1_min_10,self.u1_max_10],
        [self.u1_min_11,self.u1_max_11],[self.u1_min_12,self.u1_max_12],[self.u1_min_13,self.u1_max_13],
        [self.u1_min_14,self.u1_max_14],[self.u1_min_15,self.u1_max_15],[self.u1_min_16,self.u1_max_16],
        [self.u1_min_17,self.u1_max_17],[self.u1_min_18,self.u1_max_18],[self.u1_min_19,self.u1_max_19],
        [self.u1_min_20,self.u1_max_20],
        ]         
        return ld_u1_bounds_gui

def ld_u2_bounds_gui(self):
    
        ld_u2_bounds_gui = [
        [self.u2_min_1,self.u2_max_1],[self.u2_min_2,self.u2_max_2],[self.u2_min_3,self.u2_max_3],
        [self.u2_min_4,self.u2_max_4],[self.u2_min_5,self.u2_max_5],[self.u2_min_6,self.u2_max_6],
        [self.u2_min_7,self.u2_max_7],[self.u2_min_8,self.u2_max_8],[self.u2_min_9,self.u2_max_9],
        [self.u2_min_10,self.u2_max_10],
        [self.u2_min_11,self.u2_max_11],[self.u2_min_12,self.u2_max_12],[self.u2_min_13,self.u2_max_13],
        [self.u2_min_14,self.u2_max_14],[self.u2_min_15,self.u2_max_15],[self.u2_min_16,self.u2_max_16],
        [self.u2_min_17,self.u2_max_17],[self.u2_min_18,self.u2_max_18],[self.u2_min_19,self.u2_max_19],
        [self.u2_min_20,self.u2_max_20]
        ]         
        return ld_u2_bounds_gui

def ld_u3_bounds_gui(self):
    
        ld_u3_bounds_gui = [
        [self.u3_min_1,self.u3_max_1],[self.u3_min_2,self.u3_max_2],[self.u3_min_3,self.u3_max_3],
        [self.u3_min_4,self.u3_max_4],[self.u3_min_5,self.u3_max_5],[self.u3_min_6,self.u3_max_6],
        [self.u3_min_7,self.u3_max_7],[self.u3_min_8,self.u3_max_8],[self.u3_min_9,self.u3_max_9],
        [self.u3_min_10,self.u3_max_10],
        [self.u3_min_11,self.u3_max_11],[self.u3_min_12,self.u3_max_12],[self.u3_min_13,self.u3_max_13],
        [self.u3_min_14,self.u3_max_14],[self.u3_min_15,self.u3_max_15],[self.u3_min_16,self.u3_max_16],
        [self.u3_min_17,self.u3_max_17],[self.u3_min_18,self.u3_max_18],[self.u3_min_19,self.u3_max_19],
        [self.u3_min_20,self.u3_max_20]
        ]         
        return ld_u3_bounds_gui


def ld_u4_bounds_gui(self):
    
        ld_u4_bounds_gui = [
        [self.u4_min_1,self.u4_max_1],[self.u4_min_2,self.u4_max_2],[self.u4_min_3,self.u4_max_3],
        [self.u4_min_4,self.u4_max_4],[self.u4_min_5,self.u4_max_5],[self.u4_min_6,self.u4_max_6],
        [self.u4_min_7,self.u4_max_7],[self.u4_min_8,self.u4_max_8],[self.u4_min_9,self.u4_max_9],
        [self.u4_min_10,self.u4_max_10],
        [self.u4_min_11,self.u4_max_11],[self.u4_min_12,self.u4_max_12],[self.u4_min_13,self.u4_max_13],
        [self.u4_min_14,self.u4_max_14],[self.u4_min_15,self.u4_max_15],[self.u4_min_16,self.u4_max_16],
        [self.u4_min_17,self.u4_max_17],[self.u4_min_18,self.u4_max_18],[self.u4_min_19,self.u4_max_19],
        [self.u4_min_20,self.u4_max_20]
        ]         
        return ld_u4_bounds_gui
  

############################

def ld_u1_norm_pr_gui(self):
    
        ld_u1_norm_pr_gui = [
        [self.u1_mean_1,self.u1_sigma_1,self.use_u1_mean_nr_1],
        [self.u1_mean_2,self.u1_sigma_2,self.use_u1_mean_nr_2],
        [self.u1_mean_3,self.u1_sigma_3,self.use_u1_mean_nr_3],
        [self.u1_mean_4,self.u1_sigma_4,self.use_u1_mean_nr_4],
        [self.u1_mean_5,self.u1_sigma_5,self.use_u1_mean_nr_5],
        [self.u1_mean_6,self.u1_sigma_6,self.use_u1_mean_nr_6],
        [self.u1_mean_7,self.u1_sigma_7,self.use_u1_mean_nr_7],
        [self.u1_mean_8,self.u1_sigma_8,self.use_u1_mean_nr_8],
        [self.u1_mean_9,self.u1_sigma_9,self.use_u1_mean_nr_9],
        [self.u1_mean_10,self.u1_sigma_10,self.use_u1_mean_nr_10],
        [self.u1_mean_11,self.u1_sigma_11,self.use_u1_mean_nr_11],
        [self.u1_mean_12,self.u1_sigma_12,self.use_u1_mean_nr_12],
        [self.u1_mean_13,self.u1_sigma_13,self.use_u1_mean_nr_13],
        [self.u1_mean_14,self.u1_sigma_14,self.use_u1_mean_nr_14],
        [self.u1_mean_15,self.u1_sigma_15,self.use_u1_mean_nr_15],
        [self.u1_mean_16,self.u1_sigma_16,self.use_u1_mean_nr_16],
        [self.u1_mean_17,self.u1_sigma_17,self.use_u1_mean_nr_17],
        [self.u1_mean_18,self.u1_sigma_18,self.use_u1_mean_nr_18],
        [self.u1_mean_19,self.u1_sigma_19,self.use_u1_mean_nr_19],
        [self.u1_mean_20,self.u1_sigma_20,self.use_u1_mean_nr_20],
        ]         
        return ld_u1_norm_pr_gui

def ld_u2_norm_pr_gui(self):
    
        ld_u2_norm_pr_gui = [
        [self.u2_mean_1,self.u2_sigma_1,self.use_u2_mean_nr_1],
        [self.u2_mean_2,self.u2_sigma_2,self.use_u2_mean_nr_2],
        [self.u2_mean_3,self.u2_sigma_3,self.use_u2_mean_nr_3],
        [self.u2_mean_4,self.u2_sigma_4,self.use_u2_mean_nr_4],
        [self.u2_mean_5,self.u2_sigma_5,self.use_u2_mean_nr_5],
        [self.u2_mean_6,self.u2_sigma_6,self.use_u2_mean_nr_6],
        [self.u2_mean_7,self.u2_sigma_7,self.use_u2_mean_nr_7],
        [self.u2_mean_8,self.u2_sigma_8,self.use_u2_mean_nr_8],
        [self.u2_mean_9,self.u2_sigma_9,self.use_u2_mean_nr_9],
        [self.u2_mean_10,self.u2_sigma_10,self.use_u2_mean_nr_10],
        [self.u2_mean_11,self.u2_sigma_11,self.use_u2_mean_nr_11],
        [self.u2_mean_12,self.u2_sigma_12,self.use_u2_mean_nr_12],
        [self.u2_mean_13,self.u2_sigma_13,self.use_u2_mean_nr_13],
        [self.u2_mean_14,self.u2_sigma_14,self.use_u2_mean_nr_14],
        [self.u2_mean_15,self.u2_sigma_15,self.use_u2_mean_nr_15],
        [self.u2_mean_16,self.u2_sigma_16,self.use_u2_mean_nr_16],
        [self.u2_mean_17,self.u2_sigma_17,self.use_u2_mean_nr_17],
        [self.u2_mean_18,self.u2_sigma_18,self.use_u2_mean_nr_18],
        [self.u2_mean_19,self.u2_sigma_19,self.use_u2_mean_nr_19],
        [self.u2_mean_20,self.u2_sigma_20,self.use_u2_mean_nr_20],
        ]         
        return ld_u2_norm_pr_gui

def ld_u3_norm_pr_gui(self):
    
        ld_u3_norm_pr_gui = [
        [self.u3_mean_1,self.u3_sigma_1,self.use_u3_mean_nr_1],
        [self.u3_mean_2,self.u3_sigma_2,self.use_u3_mean_nr_2],
        [self.u3_mean_3,self.u3_sigma_3,self.use_u3_mean_nr_3],
        [self.u3_mean_4,self.u3_sigma_4,self.use_u3_mean_nr_4],
        [self.u3_mean_5,self.u3_sigma_5,self.use_u3_mean_nr_5],
        [self.u3_mean_6,self.u3_sigma_6,self.use_u3_mean_nr_6],
        [self.u3_mean_7,self.u3_sigma_7,self.use_u3_mean_nr_7],
        [self.u3_mean_8,self.u3_sigma_8,self.use_u3_mean_nr_8],
        [self.u3_mean_9,self.u3_sigma_9,self.use_u3_mean_nr_9],
        [self.u3_mean_10,self.u3_sigma_10,self.use_u3_mean_nr_10],
        [self.u3_mean_11,self.u3_sigma_11,self.use_u3_mean_nr_11],
        [self.u3_mean_12,self.u3_sigma_12,self.use_u3_mean_nr_12],
        [self.u3_mean_13,self.u3_sigma_13,self.use_u3_mean_nr_13],
        [self.u3_mean_14,self.u3_sigma_14,self.use_u3_mean_nr_14],
        [self.u3_mean_15,self.u3_sigma_15,self.use_u3_mean_nr_15],
        [self.u3_mean_16,self.u3_sigma_16,self.use_u3_mean_nr_16],
        [self.u3_mean_17,self.u3_sigma_17,self.use_u3_mean_nr_17],
        [self.u3_mean_18,self.u3_sigma_18,self.use_u3_mean_nr_18],
        [self.u3_mean_19,self.u3_sigma_19,self.use_u3_mean_nr_19],
        [self.u3_mean_20,self.u3_sigma_20,self.use_u3_mean_nr_20],
        ]         
        return ld_u3_norm_pr_gui


def ld_u4_norm_pr_gui(self):
    
        ld_u4_norm_pr_gui = [
        [self.u4_mean_1,self.u4_sigma_1,self.use_u4_mean_nr_1],
        [self.u4_mean_2,self.u4_sigma_2,self.use_u4_mean_nr_2],
        [self.u4_mean_3,self.u4_sigma_3,self.use_u4_mean_nr_3],
        [self.u4_mean_4,self.u4_sigma_4,self.use_u4_mean_nr_4],
        [self.u4_mean_5,self.u4_sigma_5,self.use_u4_mean_nr_5],
        [self.u4_mean_6,self.u4_sigma_6,self.use_u4_mean_nr_6],
        [self.u4_mean_7,self.u4_sigma_7,self.use_u4_mean_nr_7],
        [self.u4_mean_8,self.u4_sigma_8,self.use_u4_mean_nr_8],
        [self.u4_mean_9,self.u4_sigma_9,self.use_u4_mean_nr_9],
        [self.u4_mean_10,self.u4_sigma_10,self.use_u4_mean_nr_10],
        [self.u4_mean_11,self.u4_sigma_11,self.use_u4_mean_nr_11],
        [self.u4_mean_12,self.u4_sigma_12,self.use_u4_mean_nr_12],
        [self.u4_mean_13,self.u4_sigma_13,self.use_u4_mean_nr_13],
        [self.u4_mean_14,self.u4_sigma_14,self.use_u4_mean_nr_14],
        [self.u4_mean_15,self.u4_sigma_15,self.use_u4_mean_nr_15],
        [self.u4_mean_16,self.u4_sigma_16,self.use_u4_mean_nr_16],
        [self.u4_mean_17,self.u4_sigma_17,self.use_u4_mean_nr_17],
        [self.u4_mean_18,self.u4_sigma_18,self.use_u4_mean_nr_18],
        [self.u4_mean_19,self.u4_sigma_19,self.use_u4_mean_nr_19],
        [self.u4_mean_20,self.u4_sigma_20,self.use_u4_mean_nr_20],
        ]         
        return ld_u4_norm_pr_gui


############################

def ld_u1_jeff_pr_gui(self):
    
        ld_u1_jeff_pr_gui = [
        [self.u1_alpha_1,self.u1_beta_1,self.use_u1_jeff_1],
        [self.u1_alpha_2,self.u1_beta_2,self.use_u1_jeff_2],
        [self.u1_alpha_3,self.u1_beta_3,self.use_u1_jeff_3],
        [self.u1_alpha_4,self.u1_beta_4,self.use_u1_jeff_4],
        [self.u1_alpha_5,self.u1_beta_5,self.use_u1_jeff_5],
        [self.u1_alpha_6,self.u1_beta_6,self.use_u1_jeff_6],
        [self.u1_alpha_7,self.u1_beta_7,self.use_u1_jeff_7],
        [self.u1_alpha_8,self.u1_beta_8,self.use_u1_jeff_8],
        [self.u1_alpha_9,self.u1_beta_9,self.use_u1_jeff_9],
        [self.u1_alpha_10,self.u1_beta_10,self.use_u1_jeff_10],
        [self.u1_alpha_11,self.u1_beta_11,self.use_u1_jeff_11],
        [self.u1_alpha_12,self.u1_beta_12,self.use_u1_jeff_12],
        [self.u1_alpha_13,self.u1_beta_13,self.use_u1_jeff_13],
        [self.u1_alpha_14,self.u1_beta_14,self.use_u1_jeff_14],
        [self.u1_alpha_15,self.u1_beta_15,self.use_u1_jeff_15],
        [self.u1_alpha_16,self.u1_beta_16,self.use_u1_jeff_16],
        [self.u1_alpha_17,self.u1_beta_17,self.use_u1_jeff_17],
        [self.u1_alpha_18,self.u1_beta_18,self.use_u1_jeff_18],
        [self.u1_alpha_19,self.u1_beta_19,self.use_u1_jeff_19],
        [self.u1_alpha_20,self.u1_beta_20,self.use_u1_jeff_20],
        ]         
        return ld_u1_jeff_pr_gui

def ld_u2_jeff_pr_gui(self):
    
        ld_u2_jeff_pr_gui = [
        [self.u2_alpha_1,self.u2_beta_1,self.use_u2_jeff_1],
        [self.u2_alpha_2,self.u2_beta_2,self.use_u2_jeff_2],
        [self.u2_alpha_3,self.u2_beta_3,self.use_u2_jeff_3],
        [self.u2_alpha_4,self.u2_beta_4,self.use_u2_jeff_4],
        [self.u2_alpha_5,self.u2_beta_5,self.use_u2_jeff_5],
        [self.u2_alpha_6,self.u2_beta_6,self.use_u2_jeff_6],
        [self.u2_alpha_7,self.u2_beta_7,self.use_u2_jeff_7],
        [self.u2_alpha_8,self.u2_beta_8,self.use_u2_jeff_8],
        [self.u2_alpha_9,self.u2_beta_9,self.use_u2_jeff_9],
        [self.u2_alpha_10,self.u2_beta_10,self.use_u2_jeff_10],
        [self.u2_alpha_11,self.u2_beta_11,self.use_u2_jeff_11],
        [self.u2_alpha_12,self.u2_beta_12,self.use_u2_jeff_12],
        [self.u2_alpha_13,self.u2_beta_13,self.use_u2_jeff_13],
        [self.u2_alpha_14,self.u2_beta_14,self.use_u2_jeff_14],
        [self.u2_alpha_15,self.u2_beta_15,self.use_u2_jeff_15],
        [self.u2_alpha_16,self.u2_beta_16,self.use_u2_jeff_16],
        [self.u2_alpha_17,self.u2_beta_17,self.use_u2_jeff_17],
        [self.u2_alpha_18,self.u2_beta_18,self.use_u2_jeff_18],
        [self.u2_alpha_19,self.u2_beta_19,self.use_u2_jeff_19],
        [self.u2_alpha_20,self.u2_beta_20,self.use_u2_jeff_20],
        ]         
        return ld_u2_jeff_pr_gui

def ld_u3_jeff_pr_gui(self):
    
        ld_u3_jeff_pr_gui = [
        [self.u3_alpha_1,self.u3_beta_1,self.use_u3_jeff_1],
        [self.u3_alpha_2,self.u3_beta_2,self.use_u3_jeff_2],
        [self.u3_alpha_3,self.u3_beta_3,self.use_u3_jeff_3],
        [self.u3_alpha_4,self.u3_beta_4,self.use_u3_jeff_4],
        [self.u3_alpha_5,self.u3_beta_5,self.use_u3_jeff_5],
        [self.u3_alpha_6,self.u3_beta_6,self.use_u3_jeff_6],
        [self.u3_alpha_7,self.u3_beta_7,self.use_u3_jeff_7],
        [self.u3_alpha_8,self.u3_beta_8,self.use_u3_jeff_8],
        [self.u3_alpha_9,self.u3_beta_9,self.use_u3_jeff_9],
        [self.u3_alpha_10,self.u3_beta_10,self.use_u3_jeff_10],
        [self.u3_alpha_11,self.u3_beta_11,self.use_u3_jeff_11],
        [self.u3_alpha_12,self.u3_beta_12,self.use_u3_jeff_12],
        [self.u3_alpha_13,self.u3_beta_13,self.use_u3_jeff_13],
        [self.u3_alpha_14,self.u3_beta_14,self.use_u3_jeff_14],
        [self.u3_alpha_15,self.u3_beta_15,self.use_u3_jeff_15],
        [self.u3_alpha_16,self.u3_beta_16,self.use_u3_jeff_16],
        [self.u3_alpha_17,self.u3_beta_17,self.use_u3_jeff_17],
        [self.u3_alpha_18,self.u3_beta_18,self.use_u3_jeff_18],
        [self.u3_alpha_19,self.u3_beta_19,self.use_u3_jeff_19],
        [self.u3_alpha_20,self.u3_beta_20,self.use_u3_jeff_20],
        ]         
        return ld_u3_jeff_pr_gui


def ld_u4_jeff_pr_gui(self):
    
        ld_u4_jeff_pr_gui = [
        [self.u4_alpha_1,self.u4_beta_1,self.use_u4_jeff_1],
        [self.u4_alpha_2,self.u4_beta_2,self.use_u4_jeff_2],
        [self.u4_alpha_3,self.u4_beta_3,self.use_u4_jeff_3],
        [self.u4_alpha_4,self.u4_beta_4,self.use_u4_jeff_4],
        [self.u4_alpha_5,self.u4_beta_5,self.use_u4_jeff_5],
        [self.u4_alpha_6,self.u4_beta_6,self.use_u4_jeff_6],
        [self.u4_alpha_7,self.u4_beta_7,self.use_u4_jeff_7],
        [self.u4_alpha_8,self.u4_beta_8,self.use_u4_jeff_8],
        [self.u4_alpha_9,self.u4_beta_9,self.use_u4_jeff_9],
        [self.u4_alpha_10,self.u4_beta_10,self.use_u4_jeff_10],
        [self.u4_alpha_11,self.u4_beta_11,self.use_u4_jeff_11],
        [self.u4_alpha_12,self.u4_beta_12,self.use_u4_jeff_12],
        [self.u4_alpha_13,self.u4_beta_13,self.use_u4_jeff_13],
        [self.u4_alpha_14,self.u4_beta_14,self.use_u4_jeff_14],
        [self.u4_alpha_15,self.u4_beta_15,self.use_u4_jeff_15],
        [self.u4_alpha_16,self.u4_beta_16,self.use_u4_jeff_16],
        [self.u4_alpha_17,self.u4_beta_17,self.use_u4_jeff_17],
        [self.u4_alpha_18,self.u4_beta_18,self.use_u4_jeff_18],
        [self.u4_alpha_19,self.u4_beta_19,self.use_u4_jeff_19],
        [self.u4_alpha_20,self.u4_beta_20,self.use_u4_jeff_20],
        ]         
        return ld_u4_jeff_pr_gui
    
################# Bounds (Flat Prior) ################
    
    
def tra_lintr_bounds_gui(self):


        tra_lintr_bounds_gui = [
        [self.tra_lin_trend_bound_min_1,self.tra_lin_trend_bound_max_1], 
        [self.tra_lin_trend_bound_min_2,self.tra_lin_trend_bound_max_2], 
        [self.tra_lin_trend_bound_min_3,self.tra_lin_trend_bound_max_3], 
        [self.tra_lin_trend_bound_min_4,self.tra_lin_trend_bound_max_4], 
        [self.tra_lin_trend_bound_min_5,self.tra_lin_trend_bound_max_5], 
        [self.tra_lin_trend_bound_min_6,self.tra_lin_trend_bound_max_6], 
        [self.tra_lin_trend_bound_min_7,self.tra_lin_trend_bound_max_7], 
        [self.tra_lin_trend_bound_min_8,self.tra_lin_trend_bound_max_8], 
        [self.tra_lin_trend_bound_min_9,self.tra_lin_trend_bound_max_9], 
        [self.tra_lin_trend_bound_min_10,self.tra_lin_trend_bound_max_10],
        [self.tra_lin_trend_bound_min_11,self.tra_lin_trend_bound_max_11], 
        [self.tra_lin_trend_bound_min_12,self.tra_lin_trend_bound_max_12], 
        [self.tra_lin_trend_bound_min_13,self.tra_lin_trend_bound_max_13], 
        [self.tra_lin_trend_bound_min_14,self.tra_lin_trend_bound_max_14], 
        [self.tra_lin_trend_bound_min_15,self.tra_lin_trend_bound_max_15], 
        [self.tra_lin_trend_bound_min_16,self.tra_lin_trend_bound_max_16], 
        [self.tra_lin_trend_bound_min_17,self.tra_lin_trend_bound_max_17], 
        [self.tra_lin_trend_bound_min_18,self.tra_lin_trend_bound_max_18], 
        [self.tra_lin_trend_bound_min_19,self.tra_lin_trend_bound_max_19], 
        [self.tra_lin_trend_bound_min_20,self.tra_lin_trend_bound_max_20]
        ]
        
        return tra_lintr_bounds_gui
    
def tra_quadtr_bounds_gui(self):


        tra_quadtr_bounds_gui = [
        [self.tra_quad_trend_bound_min_1,self.tra_quad_trend_bound_max_1], 
        [self.tra_quad_trend_bound_min_2,self.tra_quad_trend_bound_max_2], 
        [self.tra_quad_trend_bound_min_3,self.tra_quad_trend_bound_max_3], 
        [self.tra_quad_trend_bound_min_4,self.tra_quad_trend_bound_max_4], 
        [self.tra_quad_trend_bound_min_5,self.tra_quad_trend_bound_max_5], 
        [self.tra_quad_trend_bound_min_6,self.tra_quad_trend_bound_max_6], 
        [self.tra_quad_trend_bound_min_7,self.tra_quad_trend_bound_max_7], 
        [self.tra_quad_trend_bound_min_8,self.tra_quad_trend_bound_max_8], 
        [self.tra_quad_trend_bound_min_9,self.tra_quad_trend_bound_max_9], 
        [self.tra_quad_trend_bound_min_10,self.tra_quad_trend_bound_max_10],
        [self.tra_quad_trend_bound_min_11,self.tra_quad_trend_bound_max_11], 
        [self.tra_quad_trend_bound_min_12,self.tra_quad_trend_bound_max_12], 
        [self.tra_quad_trend_bound_min_13,self.tra_quad_trend_bound_max_13], 
        [self.tra_quad_trend_bound_min_14,self.tra_quad_trend_bound_max_14], 
        [self.tra_quad_trend_bound_min_15,self.tra_quad_trend_bound_max_15], 
        [self.tra_quad_trend_bound_min_16,self.tra_quad_trend_bound_max_16], 
        [self.tra_quad_trend_bound_min_17,self.tra_quad_trend_bound_max_17], 
        [self.tra_quad_trend_bound_min_18,self.tra_quad_trend_bound_max_18], 
        [self.tra_quad_trend_bound_min_19,self.tra_quad_trend_bound_max_19], 
        [self.tra_quad_trend_bound_min_20,self.tra_quad_trend_bound_max_20]
        ]
        
        return tra_quadtr_bounds_gui





################# Normal Prior ################
    


def param_nr_priors_gui(self):
    
        param_nr_priors_gui = [
        [self.K_mean_1,self.K_sigma_1,self.use_K_norm_pr_1],[self.P_mean_1,self.P_sigma_1,self.use_P_norm_pr_1], [self.e_mean_1,self.e_sigma_1,self.use_e_norm_pr_1],[self.om_mean_1,self.om_sigma_1,self.use_om_norm_pr_1], [self.ma_mean_1,self.ma_sigma_1,self.use_ma_norm_pr_1],[self.incl_mean_1,self.incl_sigma_1,self.use_incl_norm_pr_1], [self.Omega_mean_1,self.Omega_sigma_1, self.use_Omega_norm_pr_1],[self.t0_mean_1,self.t0_sigma_1, self.use_t0_norm_pr_1],[self.pl_rad_mean_1,self.pl_rad_sigma_1,self.use_pl_rad_norm_pr_1],[self.a_sol_mean_1,self.a_sol_sigma_1,self.use_a_sol_norm_pr_1],
        [self.K_mean_2,self.K_sigma_2,self.use_K_norm_pr_2],[self.P_mean_2,self.P_sigma_2,self.use_P_norm_pr_2], [self.e_mean_2,self.e_sigma_2,self.use_e_norm_pr_2],[self.om_mean_2,self.om_sigma_2,self.use_om_norm_pr_2], [self.ma_mean_2,self.ma_sigma_2,self.use_ma_norm_pr_2],[self.incl_mean_2,self.incl_sigma_2,self.use_incl_norm_pr_2], [self.Omega_mean_2,self.Omega_sigma_2, self.use_Omega_norm_pr_2],[self.t0_mean_2,self.t0_sigma_2, self.use_t0_norm_pr_2],[self.pl_rad_mean_2,self.pl_rad_sigma_2,self.use_pl_rad_norm_pr_2],[self.a_sol_mean_2,self.a_sol_sigma_2,self.use_a_sol_norm_pr_2],
        [self.K_mean_3,self.K_sigma_3,self.use_K_norm_pr_3],[self.P_mean_3,self.P_sigma_3,self.use_P_norm_pr_3], [self.e_mean_3,self.e_sigma_3,self.use_e_norm_pr_3],[self.om_mean_3,self.om_sigma_3,self.use_om_norm_pr_3], [self.ma_mean_3,self.ma_sigma_3,self.use_ma_norm_pr_3],[self.incl_mean_3,self.incl_sigma_3,self.use_incl_norm_pr_3], [self.Omega_mean_3,self.Omega_sigma_3, self.use_Omega_norm_pr_3],[self.t0_mean_3,self.t0_sigma_3, self.use_t0_norm_pr_3],[self.pl_rad_mean_3,self.pl_rad_sigma_3,self.use_pl_rad_norm_pr_3],[self.a_sol_mean_3,self.a_sol_sigma_3,self.use_a_sol_norm_pr_3],
        [self.K_mean_4,self.K_sigma_4,self.use_K_norm_pr_4],[self.P_mean_4,self.P_sigma_4,self.use_P_norm_pr_4], [self.e_mean_4,self.e_sigma_4,self.use_e_norm_pr_4],[self.om_mean_4,self.om_sigma_4,self.use_om_norm_pr_4], [self.ma_mean_4,self.ma_sigma_4,self.use_ma_norm_pr_4],[self.incl_mean_4,self.incl_sigma_4,self.use_incl_norm_pr_4], [self.Omega_mean_4,self.Omega_sigma_4, self.use_Omega_norm_pr_4],[self.t0_mean_4,self.t0_sigma_4, self.use_t0_norm_pr_4],[self.pl_rad_mean_4,self.pl_rad_sigma_4,self.use_pl_rad_norm_pr_4],[self.a_sol_mean_4,self.a_sol_sigma_4,self.use_a_sol_norm_pr_4],
        [self.K_mean_5,self.K_sigma_5,self.use_K_norm_pr_5],[self.P_mean_5,self.P_sigma_5,self.use_P_norm_pr_5], [self.e_mean_5,self.e_sigma_5,self.use_e_norm_pr_5],[self.om_mean_5,self.om_sigma_5,self.use_om_norm_pr_5], [self.ma_mean_5,self.ma_sigma_5,self.use_ma_norm_pr_5],[self.incl_mean_5,self.incl_sigma_5,self.use_incl_norm_pr_5], [self.Omega_mean_5,self.Omega_sigma_5, self.use_Omega_norm_pr_5],[self.t0_mean_5,self.t0_sigma_5, self.use_t0_norm_pr_5],[self.pl_rad_mean_5,self.pl_rad_sigma_5,self.use_pl_rad_norm_pr_5],[self.a_sol_mean_5,self.a_sol_sigma_5,self.use_a_sol_norm_pr_5],
        [self.K_mean_6,self.K_sigma_6,self.use_K_norm_pr_6],[self.P_mean_6,self.P_sigma_6,self.use_P_norm_pr_6], [self.e_mean_6,self.e_sigma_6,self.use_e_norm_pr_6],[self.om_mean_6,self.om_sigma_6,self.use_om_norm_pr_6], [self.ma_mean_6,self.ma_sigma_6,self.use_ma_norm_pr_6],[self.incl_mean_6,self.incl_sigma_6,self.use_incl_norm_pr_6], [self.Omega_mean_6,self.Omega_sigma_6, self.use_Omega_norm_pr_6],[self.t0_mean_6,self.t0_sigma_6, self.use_t0_norm_pr_6],[self.pl_rad_mean_6,self.pl_rad_sigma_6,self.use_pl_rad_norm_pr_6],[self.a_sol_mean_6,self.a_sol_sigma_6,self.use_a_sol_norm_pr_6],
        [self.K_mean_7,self.K_sigma_7,self.use_K_norm_pr_7],[self.P_mean_7,self.P_sigma_7,self.use_P_norm_pr_7], [self.e_mean_7,self.e_sigma_7,self.use_e_norm_pr_7],[self.om_mean_7,self.om_sigma_7,self.use_om_norm_pr_7], [self.ma_mean_7,self.ma_sigma_7,self.use_ma_norm_pr_7],[self.incl_mean_7,self.incl_sigma_7,self.use_incl_norm_pr_7], [self.Omega_mean_7,self.Omega_sigma_7, self.use_Omega_norm_pr_7],[self.t0_mean_7,self.t0_sigma_7, self.use_t0_norm_pr_7],[self.pl_rad_mean_7,self.pl_rad_sigma_7,self.use_pl_rad_norm_pr_7],[self.a_sol_mean_7,self.a_sol_sigma_7,self.use_a_sol_norm_pr_7],
        [self.K_mean_8,self.K_sigma_8,self.use_K_norm_pr_8],[self.P_mean_8,self.P_sigma_8,self.use_P_norm_pr_8], [self.e_mean_8,self.e_sigma_8,self.use_e_norm_pr_8],[self.om_mean_8,self.om_sigma_8,self.use_om_norm_pr_8], [self.ma_mean_8,self.ma_sigma_8,self.use_ma_norm_pr_8],[self.incl_mean_8,self.incl_sigma_8,self.use_incl_norm_pr_8], [self.Omega_mean_8,self.Omega_sigma_8, self.use_Omega_norm_pr_8],[self.t0_mean_8,self.t0_sigma_8, self.use_t0_norm_pr_8],[self.pl_rad_mean_8,self.pl_rad_sigma_8,self.use_pl_rad_norm_pr_8],[self.a_sol_mean_8,self.a_sol_sigma_8,self.use_a_sol_norm_pr_8],
        [self.K_mean_9,self.K_sigma_9,self.use_K_norm_pr_9],[self.P_mean_9,self.P_sigma_9,self.use_P_norm_pr_9], [self.e_mean_9,self.e_sigma_9,self.use_e_norm_pr_9],[self.om_mean_9,self.om_sigma_9,self.use_om_norm_pr_9], [self.ma_mean_9,self.ma_sigma_9,self.use_ma_norm_pr_9],[self.incl_mean_9,self.incl_sigma_9,self.use_incl_norm_pr_9], [self.Omega_mean_9,self.Omega_sigma_9, self.use_Omega_norm_pr_9],[self.t0_mean_9,self.t0_sigma_9, self.use_t0_norm_pr_9],[self.pl_rad_mean_9,self.pl_rad_sigma_9,self.use_pl_rad_norm_pr_9],[self.a_sol_mean_9,self.a_sol_sigma_9,self.use_a_sol_norm_pr_9],
        ]
        return param_nr_priors_gui



def offset_nr_priors_gui(self):


        offset_nr_priors_gui = [
        [self.RV_Data_mean_1,self.RV_Data_sigma_1,self.use_rvoff_nr_1], 
        [self.RV_Data_mean_2,self.RV_Data_sigma_2,self.use_rvoff_nr_2], 
        [self.RV_Data_mean_3,self.RV_Data_sigma_3,self.use_rvoff_nr_3], 
        [self.RV_Data_mean_4,self.RV_Data_sigma_4,self.use_rvoff_nr_4], 
        [self.RV_Data_mean_5,self.RV_Data_sigma_5,self.use_rvoff_nr_5], 
        [self.RV_Data_mean_6,self.RV_Data_sigma_6,self.use_rvoff_nr_6], 
        [self.RV_Data_mean_7,self.RV_Data_sigma_7,self.use_rvoff_nr_7], 
        [self.RV_Data_mean_8,self.RV_Data_sigma_8,self.use_rvoff_nr_8], 
        [self.RV_Data_mean_9,self.RV_Data_sigma_9,self.use_rvoff_nr_9], 
        [self.RV_Data_mean_10,self.RV_Data_sigma_10,self.use_rvoff_nr_10],
        [self.RV_Data_mean_11,self.RV_Data_sigma_11,self.use_rvoff_nr_11],
        [self.RV_Data_mean_12,self.RV_Data_sigma_12,self.use_rvoff_nr_12],
        [self.RV_Data_mean_13,self.RV_Data_sigma_13,self.use_rvoff_nr_13],
        [self.RV_Data_mean_14,self.RV_Data_sigma_14,self.use_rvoff_nr_14],
        [self.RV_Data_mean_15,self.RV_Data_sigma_15,self.use_rvoff_nr_15],
        [self.RV_Data_mean_16,self.RV_Data_sigma_16,self.use_rvoff_nr_16],
        [self.RV_Data_mean_17,self.RV_Data_sigma_17,self.use_rvoff_nr_17],
        [self.RV_Data_mean_18,self.RV_Data_sigma_18,self.use_rvoff_nr_18],
        [self.RV_Data_mean_19,self.RV_Data_sigma_19,self.use_rvoff_nr_19],
        [self.RV_Data_mean_20,self.RV_Data_sigma_20,self.use_rvoff_nr_20]
        ]
        
        return offset_nr_priors_gui
        
    
    
def jitter_nr_priors_gui(self):

        jitter_nr_priors_gui = [
        [self.RV_jitter_mean_1,self.RV_jitter_sigma_1,self.use_rvjitt_nr_1], 
        [self.RV_jitter_mean_2,self.RV_jitter_sigma_2,self.use_rvjitt_nr_2], 
        [self.RV_jitter_mean_3,self.RV_jitter_sigma_3,self.use_rvjitt_nr_3], 
        [self.RV_jitter_mean_4,self.RV_jitter_sigma_4,self.use_rvjitt_nr_4], 
        [self.RV_jitter_mean_5,self.RV_jitter_sigma_5,self.use_rvjitt_nr_5], 
        [self.RV_jitter_mean_6,self.RV_jitter_sigma_6,self.use_rvjitt_nr_6],
        [self.RV_jitter_mean_7,self.RV_jitter_sigma_7,self.use_rvjitt_nr_7], 
        [self.RV_jitter_mean_8,self.RV_jitter_sigma_8,self.use_rvjitt_nr_8], 
        [self.RV_jitter_mean_9,self.RV_jitter_sigma_9,self.use_rvjitt_nr_9], 
        [self.RV_jitter_mean_10,self.RV_jitter_sigma_10,self.use_rvjitt_nr_10],
        [self.RV_jitter_mean_11,self.RV_jitter_sigma_11,self.use_rvjitt_nr_11],
        [self.RV_jitter_mean_12,self.RV_jitter_sigma_12,self.use_rvjitt_nr_12],
        [self.RV_jitter_mean_13,self.RV_jitter_sigma_13,self.use_rvjitt_nr_13],
        [self.RV_jitter_mean_14,self.RV_jitter_sigma_14,self.use_rvjitt_nr_14],
        [self.RV_jitter_mean_15,self.RV_jitter_sigma_15,self.use_rvjitt_nr_15],
        [self.RV_jitter_mean_16,self.RV_jitter_sigma_16,self.use_rvjitt_nr_16],
        [self.RV_jitter_mean_17,self.RV_jitter_sigma_17,self.use_rvjitt_nr_17],
        [self.RV_jitter_mean_18,self.RV_jitter_sigma_18,self.use_rvjitt_nr_18],
        [self.RV_jitter_mean_19,self.RV_jitter_sigma_19,self.use_rvjitt_nr_19],
        [self.RV_jitter_mean_20,self.RV_jitter_sigma_20,self.use_rvjitt_nr_20]
        ]
        
        return jitter_nr_priors_gui
    
    
def offset_nr_priors_gui_tra(self):


        offset_nr_priors_gui_tra = [
        [self.tra_Data_mean_1,self.tra_Data_sigma_1,self.use_traoff_nr_1], 
        [self.tra_Data_mean_2,self.tra_Data_sigma_2,self.use_traoff_nr_2], 
        [self.tra_Data_mean_3,self.tra_Data_sigma_3,self.use_traoff_nr_3], 
        [self.tra_Data_mean_4,self.tra_Data_sigma_4,self.use_traoff_nr_4], 
        [self.tra_Data_mean_5,self.tra_Data_sigma_5,self.use_traoff_nr_5], 
        [self.tra_Data_mean_6,self.tra_Data_sigma_6,self.use_traoff_nr_6], 
        [self.tra_Data_mean_7,self.tra_Data_sigma_7,self.use_traoff_nr_7], 
        [self.tra_Data_mean_8,self.tra_Data_sigma_8,self.use_traoff_nr_8], 
        [self.tra_Data_mean_9,self.tra_Data_sigma_9,self.use_traoff_nr_9], 
        [self.tra_Data_mean_10,self.tra_Data_sigma_10,self.use_traoff_nr_10],
        [self.tra_Data_mean_11,self.tra_Data_sigma_11,self.use_traoff_nr_11], 
        [self.tra_Data_mean_12,self.tra_Data_sigma_12,self.use_traoff_nr_12], 
        [self.tra_Data_mean_13,self.tra_Data_sigma_13,self.use_traoff_nr_13], 
        [self.tra_Data_mean_14,self.tra_Data_sigma_14,self.use_traoff_nr_14], 
        [self.tra_Data_mean_15,self.tra_Data_sigma_15,self.use_traoff_nr_15], 
        [self.tra_Data_mean_16,self.tra_Data_sigma_16,self.use_traoff_nr_16], 
        [self.tra_Data_mean_17,self.tra_Data_sigma_17,self.use_traoff_nr_17], 
        [self.tra_Data_mean_18,self.tra_Data_sigma_18,self.use_traoff_nr_18], 
        [self.tra_Data_mean_19,self.tra_Data_sigma_19,self.use_traoff_nr_19], 
        [self.tra_Data_mean_20,self.tra_Data_sigma_20,self.use_traoff_nr_20]
        ]
        
        return offset_nr_priors_gui_tra
        
    
    
def jitter_nr_priors_gui_tra(self):

        jitter_nr_priors_gui_tra = [
        [self.tra_jitter_mean_1,self.tra_jitter_sigma_1,self.use_trajitt_nr_1], 
        [self.tra_jitter_mean_2,self.tra_jitter_sigma_2,self.use_trajitt_nr_2], 
        [self.tra_jitter_mean_3,self.tra_jitter_sigma_3,self.use_trajitt_nr_3], 
        [self.tra_jitter_mean_4,self.tra_jitter_sigma_4,self.use_trajitt_nr_4], 
        [self.tra_jitter_mean_5,self.tra_jitter_sigma_5,self.use_trajitt_nr_5], 
        [self.tra_jitter_mean_6,self.tra_jitter_sigma_6,self.use_trajitt_nr_6],
        [self.tra_jitter_mean_7,self.tra_jitter_sigma_7,self.use_trajitt_nr_7], 
        [self.tra_jitter_mean_8,self.tra_jitter_sigma_8,self.use_trajitt_nr_8], 
        [self.tra_jitter_mean_9,self.tra_jitter_sigma_9,self.use_trajitt_nr_9], 
        [self.tra_jitter_mean_10,self.tra_jitter_sigma_10,self.use_trajitt_nr_10],
        [self.tra_jitter_mean_11,self.tra_jitter_sigma_11,self.use_trajitt_nr_11], 
        [self.tra_jitter_mean_12,self.tra_jitter_sigma_12,self.use_trajitt_nr_12], 
        [self.tra_jitter_mean_13,self.tra_jitter_sigma_13,self.use_trajitt_nr_13], 
        [self.tra_jitter_mean_14,self.tra_jitter_sigma_14,self.use_trajitt_nr_14], 
        [self.tra_jitter_mean_15,self.tra_jitter_sigma_15,self.use_trajitt_nr_15], 
        [self.tra_jitter_mean_16,self.tra_jitter_sigma_16,self.use_trajitt_nr_16],
        [self.tra_jitter_mean_17,self.tra_jitter_sigma_17,self.use_trajitt_nr_17], 
        [self.tra_jitter_mean_18,self.tra_jitter_sigma_18,self.use_trajitt_nr_18], 
        [self.tra_jitter_mean_19,self.tra_jitter_sigma_19,self.use_trajitt_nr_19], 
        [self.tra_jitter_mean_20,self.tra_jitter_sigma_20,self.use_trajitt_nr_20]

        ]
        
        return jitter_nr_priors_gui_tra


def tra_lin_trend_nr_priors_gui(self):

        tra_lin_trend_nr_priors_gui = [
        [self.tra_lin_trend_mean_1,self.tra_lin_trend_sigma_1,self.use_tra_lin_mean_nr_1], 
        [self.tra_lin_trend_mean_2,self.tra_lin_trend_sigma_2,self.use_tra_lin_mean_nr_2], 
        [self.tra_lin_trend_mean_3,self.tra_lin_trend_sigma_3,self.use_tra_lin_mean_nr_3], 
        [self.tra_lin_trend_mean_4,self.tra_lin_trend_sigma_4,self.use_tra_lin_mean_nr_4], 
        [self.tra_lin_trend_mean_5,self.tra_lin_trend_sigma_5,self.use_tra_lin_mean_nr_5], 
        [self.tra_lin_trend_mean_6,self.tra_lin_trend_sigma_6,self.use_tra_lin_mean_nr_6],
        [self.tra_lin_trend_mean_7,self.tra_lin_trend_sigma_7,self.use_tra_lin_mean_nr_7], 
        [self.tra_lin_trend_mean_8,self.tra_lin_trend_sigma_8,self.use_tra_lin_mean_nr_8], 
        [self.tra_lin_trend_mean_9,self.tra_lin_trend_sigma_9,self.use_tra_lin_mean_nr_9], 
        [self.tra_lin_trend_mean_10,self.tra_lin_trend_sigma_10,self.use_tra_lin_mean_nr_10],
        [self.tra_lin_trend_mean_11,self.tra_lin_trend_sigma_11,self.use_tra_lin_mean_nr_11], 
        [self.tra_lin_trend_mean_12,self.tra_lin_trend_sigma_12,self.use_tra_lin_mean_nr_12], 
        [self.tra_lin_trend_mean_13,self.tra_lin_trend_sigma_13,self.use_tra_lin_mean_nr_13], 
        [self.tra_lin_trend_mean_14,self.tra_lin_trend_sigma_14,self.use_tra_lin_mean_nr_14], 
        [self.tra_lin_trend_mean_15,self.tra_lin_trend_sigma_15,self.use_tra_lin_mean_nr_15], 
        [self.tra_lin_trend_mean_16,self.tra_lin_trend_sigma_16,self.use_tra_lin_mean_nr_16],
        [self.tra_lin_trend_mean_17,self.tra_lin_trend_sigma_17,self.use_tra_lin_mean_nr_17], 
        [self.tra_lin_trend_mean_18,self.tra_lin_trend_sigma_18,self.use_tra_lin_mean_nr_18], 
        [self.tra_lin_trend_mean_19,self.tra_lin_trend_sigma_19,self.use_tra_lin_mean_nr_19], 
        [self.tra_lin_trend_mean_20,self.tra_lin_trend_sigma_20,self.use_tra_lin_mean_nr_20]
        ]
        
        return tra_lin_trend_nr_priors_gui
    
    
def tra_quad_trend_nr_priors_gui(self):

        tra_quad_trend_nr_priors_gui = [
        [self.tra_quad_trend_mean_1,self.tra_quad_trend_sigma_1,self.use_tra_quad_mean_nr_1], 
        [self.tra_quad_trend_mean_2,self.tra_quad_trend_sigma_2,self.use_tra_quad_mean_nr_2], 
        [self.tra_quad_trend_mean_3,self.tra_quad_trend_sigma_3,self.use_tra_quad_mean_nr_3], 
        [self.tra_quad_trend_mean_4,self.tra_quad_trend_sigma_4,self.use_tra_quad_mean_nr_4], 
        [self.tra_quad_trend_mean_5,self.tra_quad_trend_sigma_5,self.use_tra_quad_mean_nr_5], 
        [self.tra_quad_trend_mean_6,self.tra_quad_trend_sigma_6,self.use_tra_quad_mean_nr_6],
        [self.tra_quad_trend_mean_7,self.tra_quad_trend_sigma_7,self.use_tra_quad_mean_nr_7], 
        [self.tra_quad_trend_mean_8,self.tra_quad_trend_sigma_8,self.use_tra_quad_mean_nr_8], 
        [self.tra_quad_trend_mean_9,self.tra_quad_trend_sigma_9,self.use_tra_quad_mean_nr_9], 
        [self.tra_quad_trend_mean_10,self.tra_quad_trend_sigma_10,self.use_tra_quad_mean_nr_10],
        [self.tra_quad_trend_mean_11,self.tra_quad_trend_sigma_11,self.use_tra_quad_mean_nr_11], 
        [self.tra_quad_trend_mean_12,self.tra_quad_trend_sigma_12,self.use_tra_quad_mean_nr_12], 
        [self.tra_quad_trend_mean_13,self.tra_quad_trend_sigma_13,self.use_tra_quad_mean_nr_13], 
        [self.tra_quad_trend_mean_14,self.tra_quad_trend_sigma_14,self.use_tra_quad_mean_nr_14], 
        [self.tra_quad_trend_mean_15,self.tra_quad_trend_sigma_15,self.use_tra_quad_mean_nr_15], 
        [self.tra_quad_trend_mean_16,self.tra_quad_trend_sigma_16,self.use_tra_quad_mean_nr_16],
        [self.tra_quad_trend_mean_17,self.tra_quad_trend_sigma_17,self.use_tra_quad_mean_nr_17], 
        [self.tra_quad_trend_mean_18,self.tra_quad_trend_sigma_18,self.use_tra_quad_mean_nr_18], 
        [self.tra_quad_trend_mean_19,self.tra_quad_trend_sigma_19,self.use_tra_quad_mean_nr_19], 
        [self.tra_quad_trend_mean_20,self.tra_quad_trend_sigma_20,self.use_tra_quad_mean_nr_20]

        ]
        
        return tra_quad_trend_nr_priors_gui
    

################# Jeff Prior ################
    

def param_jeff_priors_gui(self):

        param_jeff_priors_gui = [
        [self.K_jeff_alpha_1,self.K_jeff_beta_1,self.use_K_jeff_pr_1],[self.P_jeff_alpha_1,self.P_jeff_beta_1,self.use_P_jeff_pr_1], [self.e_jeff_alpha_1,self.e_jeff_beta_1,self.use_e_jeff_pr_1],[self.om_jeff_alpha_1,self.om_jeff_beta_1,self.use_om_jeff_pr_1], [self.ma_jeff_alpha_1,self.ma_jeff_beta_1,self.use_ma_jeff_pr_1],[self.incl_jeff_alpha_1,self.incl_jeff_beta_1,self.use_incl_jeff_pr_1], [self.Omega_jeff_alpha_1,self.Omega_jeff_beta_1, self.use_Omega_jeff_pr_1],[self.t0_jeff_alpha_1,self.t0_jeff_beta_1, self.use_t0_jeff_pr_1],[self.pl_rad_jeff_alpha_1,self.pl_rad_jeff_beta_1,self.use_pl_rad_jeff_pr_1],[self.a_sol_jeff_alpha_1,self.a_sol_jeff_beta_1,self.use_a_sol_jeff_pr_1],
        [self.K_jeff_alpha_2,self.K_jeff_beta_2,self.use_K_jeff_pr_2],[self.P_jeff_alpha_2,self.P_jeff_beta_2,self.use_P_jeff_pr_2], [self.e_jeff_alpha_2,self.e_jeff_beta_2,self.use_e_jeff_pr_2],[self.om_jeff_alpha_2,self.om_jeff_beta_2,self.use_om_jeff_pr_2], [self.ma_jeff_alpha_2,self.ma_jeff_beta_2,self.use_ma_jeff_pr_2],[self.incl_jeff_alpha_2,self.incl_jeff_beta_2,self.use_incl_jeff_pr_2], [self.Omega_jeff_alpha_2,self.Omega_jeff_beta_2, self.use_Omega_jeff_pr_2],[self.t0_jeff_alpha_2,self.t0_jeff_beta_2, self.use_t0_jeff_pr_2],[self.pl_rad_jeff_alpha_2,self.pl_rad_jeff_beta_2,self.use_pl_rad_jeff_pr_2],[self.a_sol_jeff_alpha_2,self.a_sol_jeff_beta_2,self.use_a_sol_jeff_pr_2],
        [self.K_jeff_alpha_3,self.K_jeff_beta_3,self.use_K_jeff_pr_3],[self.P_jeff_alpha_3,self.P_jeff_beta_3,self.use_P_jeff_pr_3], [self.e_jeff_alpha_3,self.e_jeff_beta_3,self.use_e_jeff_pr_3],[self.om_jeff_alpha_3,self.om_jeff_beta_3,self.use_om_jeff_pr_3], [self.ma_jeff_alpha_3,self.ma_jeff_beta_3,self.use_ma_jeff_pr_3],[self.incl_jeff_alpha_3,self.incl_jeff_beta_3,self.use_incl_jeff_pr_3], [self.Omega_jeff_alpha_3,self.Omega_jeff_beta_3, self.use_Omega_jeff_pr_3],[self.t0_jeff_alpha_3,self.t0_jeff_beta_3, self.use_t0_jeff_pr_3],[self.pl_rad_jeff_alpha_3,self.pl_rad_jeff_beta_3,self.use_pl_rad_jeff_pr_3],[self.a_sol_jeff_alpha_3,self.a_sol_jeff_beta_3,self.use_a_sol_jeff_pr_3],
        [self.K_jeff_alpha_4,self.K_jeff_beta_4,self.use_K_jeff_pr_4],[self.P_jeff_alpha_4,self.P_jeff_beta_4,self.use_P_jeff_pr_4], [self.e_jeff_alpha_4,self.e_jeff_beta_4,self.use_e_jeff_pr_4],[self.om_jeff_alpha_4,self.om_jeff_beta_4,self.use_om_jeff_pr_4], [self.ma_jeff_alpha_4,self.ma_jeff_beta_4,self.use_ma_jeff_pr_4],[self.incl_jeff_alpha_4,self.incl_jeff_beta_4,self.use_incl_jeff_pr_4], [self.Omega_jeff_alpha_4,self.Omega_jeff_beta_4, self.use_Omega_jeff_pr_4],[self.t0_jeff_alpha_4,self.t0_jeff_beta_4, self.use_t0_jeff_pr_4],[self.pl_rad_jeff_alpha_4,self.pl_rad_jeff_beta_4,self.use_pl_rad_jeff_pr_4],[self.a_sol_jeff_alpha_4,self.a_sol_jeff_beta_4,self.use_a_sol_jeff_pr_4],
        [self.K_jeff_alpha_5,self.K_jeff_beta_5,self.use_K_jeff_pr_5],[self.P_jeff_alpha_5,self.P_jeff_beta_5,self.use_P_jeff_pr_5], [self.e_jeff_alpha_5,self.e_jeff_beta_5,self.use_e_jeff_pr_5],[self.om_jeff_alpha_5,self.om_jeff_beta_5,self.use_om_jeff_pr_5], [self.ma_jeff_alpha_5,self.ma_jeff_beta_5,self.use_ma_jeff_pr_5],[self.incl_jeff_alpha_5,self.incl_jeff_beta_5,self.use_incl_jeff_pr_5], [self.Omega_jeff_alpha_5,self.Omega_jeff_beta_5, self.use_Omega_jeff_pr_5],[self.t0_jeff_alpha_5,self.t0_jeff_beta_5, self.use_t0_jeff_pr_5],[self.pl_rad_jeff_alpha_5,self.pl_rad_jeff_beta_5,self.use_pl_rad_jeff_pr_5],[self.a_sol_jeff_alpha_5,self.a_sol_jeff_beta_5,self.use_a_sol_jeff_pr_5],
        [self.K_jeff_alpha_6,self.K_jeff_beta_6,self.use_K_jeff_pr_6],[self.P_jeff_alpha_6,self.P_jeff_beta_6,self.use_P_jeff_pr_6], [self.e_jeff_alpha_6,self.e_jeff_beta_6,self.use_e_jeff_pr_6],[self.om_jeff_alpha_6,self.om_jeff_beta_6,self.use_om_jeff_pr_6], [self.ma_jeff_alpha_6,self.ma_jeff_beta_6,self.use_ma_jeff_pr_6],[self.incl_jeff_alpha_6,self.incl_jeff_beta_6,self.use_incl_jeff_pr_6], [self.Omega_jeff_alpha_6,self.Omega_jeff_beta_6, self.use_Omega_jeff_pr_6],[self.t0_jeff_alpha_6,self.t0_jeff_beta_6, self.use_t0_jeff_pr_6],[self.pl_rad_jeff_alpha_6,self.pl_rad_jeff_beta_6,self.use_pl_rad_jeff_pr_6],[self.a_sol_jeff_alpha_6,self.a_sol_jeff_beta_6,self.use_a_sol_jeff_pr_6],
        [self.K_jeff_alpha_7,self.K_jeff_beta_7,self.use_K_jeff_pr_7],[self.P_jeff_alpha_7,self.P_jeff_beta_7,self.use_P_jeff_pr_7], [self.e_jeff_alpha_7,self.e_jeff_beta_7,self.use_e_jeff_pr_7],[self.om_jeff_alpha_7,self.om_jeff_beta_7,self.use_om_jeff_pr_7], [self.ma_jeff_alpha_7,self.ma_jeff_beta_7,self.use_ma_jeff_pr_7],[self.incl_jeff_alpha_7,self.incl_jeff_beta_7,self.use_incl_jeff_pr_7], [self.Omega_jeff_alpha_7,self.Omega_jeff_beta_7, self.use_Omega_jeff_pr_7],[self.t0_jeff_alpha_7,self.t0_jeff_beta_7, self.use_t0_jeff_pr_7],[self.pl_rad_jeff_alpha_7,self.pl_rad_jeff_beta_7,self.use_pl_rad_jeff_pr_7],[self.a_sol_jeff_alpha_7,self.a_sol_jeff_beta_7,self.use_a_sol_jeff_pr_7],
        [self.K_jeff_alpha_8,self.K_jeff_beta_8,self.use_K_jeff_pr_8],[self.P_jeff_alpha_8,self.P_jeff_beta_8,self.use_P_jeff_pr_8], [self.e_jeff_alpha_8,self.e_jeff_beta_8,self.use_e_jeff_pr_8],[self.om_jeff_alpha_8,self.om_jeff_beta_8,self.use_om_jeff_pr_8], [self.ma_jeff_alpha_8,self.ma_jeff_beta_8,self.use_ma_jeff_pr_8],[self.incl_jeff_alpha_8,self.incl_jeff_beta_8,self.use_incl_jeff_pr_8], [self.Omega_jeff_alpha_8,self.Omega_jeff_beta_8, self.use_Omega_jeff_pr_8],[self.t0_jeff_alpha_8,self.t0_jeff_beta_8, self.use_t0_jeff_pr_8],[self.pl_rad_jeff_alpha_8,self.pl_rad_jeff_beta_8,self.use_pl_rad_jeff_pr_8],[self.a_sol_jeff_alpha_8,self.a_sol_jeff_beta_8,self.use_a_sol_jeff_pr_8],
        [self.K_jeff_alpha_9,self.K_jeff_beta_9,self.use_K_jeff_pr_9],[self.P_jeff_alpha_9,self.P_jeff_beta_9,self.use_P_jeff_pr_9], [self.e_jeff_alpha_9,self.e_jeff_beta_9,self.use_e_jeff_pr_9],[self.om_jeff_alpha_9,self.om_jeff_beta_9,self.use_om_jeff_pr_9], [self.ma_jeff_alpha_9,self.ma_jeff_beta_9,self.use_ma_jeff_pr_9],[self.incl_jeff_alpha_9,self.incl_jeff_beta_9,self.use_incl_jeff_pr_9], [self.Omega_jeff_alpha_9,self.Omega_jeff_beta_9, self.use_Omega_jeff_pr_9],[self.t0_jeff_alpha_9,self.t0_jeff_beta_9, self.use_t0_jeff_pr_9],[self.pl_rad_jeff_alpha_9,self.pl_rad_jeff_beta_9,self.use_pl_rad_jeff_pr_9],[self.a_sol_jeff_alpha_9,self.a_sol_jeff_beta_9,self.use_a_sol_jeff_pr_9],
        ]

        return param_jeff_priors_gui
    
    
    
    
def offset_jeff_priors_gui(self):


        offset_jeff_priors_gui = [
        [self.RV_Data_jeff_alpha_1,self.RV_Data_jeff_beta_1,self.use_rvoff_jeff_1], 
        [self.RV_Data_jeff_alpha_2,self.RV_Data_jeff_beta_2,self.use_rvoff_jeff_2], 
        [self.RV_Data_jeff_alpha_3,self.RV_Data_jeff_beta_3,self.use_rvoff_jeff_3], 
        [self.RV_Data_jeff_alpha_4,self.RV_Data_jeff_beta_4,self.use_rvoff_jeff_4], 
        [self.RV_Data_jeff_alpha_5,self.RV_Data_jeff_beta_5,self.use_rvoff_jeff_5], 
        [self.RV_Data_jeff_alpha_6,self.RV_Data_jeff_beta_6,self.use_rvoff_jeff_6], 
        [self.RV_Data_jeff_alpha_7,self.RV_Data_jeff_beta_7,self.use_rvoff_jeff_7], 
        [self.RV_Data_jeff_alpha_8,self.RV_Data_jeff_beta_8,self.use_rvoff_jeff_8], 
        [self.RV_Data_jeff_alpha_9,self.RV_Data_jeff_beta_9,self.use_rvoff_jeff_9], 
        [self.RV_Data_jeff_alpha_10,self.RV_Data_jeff_beta_10,self.use_rvoff_jeff_10],
        [self.RV_Data_jeff_alpha_11,self.RV_Data_jeff_beta_11,self.use_rvoff_jeff_11],
        [self.RV_Data_jeff_alpha_12,self.RV_Data_jeff_beta_12,self.use_rvoff_jeff_12],
        [self.RV_Data_jeff_alpha_13,self.RV_Data_jeff_beta_13,self.use_rvoff_jeff_13],
        [self.RV_Data_jeff_alpha_14,self.RV_Data_jeff_beta_14,self.use_rvoff_jeff_14],
        [self.RV_Data_jeff_alpha_15,self.RV_Data_jeff_beta_15,self.use_rvoff_jeff_15],
        [self.RV_Data_jeff_alpha_16,self.RV_Data_jeff_beta_16,self.use_rvoff_jeff_16],
        [self.RV_Data_jeff_alpha_17,self.RV_Data_jeff_beta_17,self.use_rvoff_jeff_17],
        [self.RV_Data_jeff_alpha_18,self.RV_Data_jeff_beta_18,self.use_rvoff_jeff_18],
        [self.RV_Data_jeff_alpha_19,self.RV_Data_jeff_beta_19,self.use_rvoff_jeff_19],
        [self.RV_Data_jeff_alpha_20,self.RV_Data_jeff_beta_20,self.use_rvoff_jeff_20]
        ]
        
        return offset_jeff_priors_gui
        
    
    
def jitter_jeff_priors_gui(self):

        jitter_jeff_priors_gui = [
        [self.RV_jitter_jeff_alpha_1,self.RV_jitter_jeff_beta_1,self.use_rvjitt_jeff_1], 
        [self.RV_jitter_jeff_alpha_2,self.RV_jitter_jeff_beta_2,self.use_rvjitt_jeff_2], 
        [self.RV_jitter_jeff_alpha_3,self.RV_jitter_jeff_beta_3,self.use_rvjitt_jeff_3], 
        [self.RV_jitter_jeff_alpha_4,self.RV_jitter_jeff_beta_4,self.use_rvjitt_jeff_4], 
        [self.RV_jitter_jeff_alpha_5,self.RV_jitter_jeff_beta_5,self.use_rvjitt_jeff_5], 
        [self.RV_jitter_jeff_alpha_6,self.RV_jitter_jeff_beta_6,self.use_rvjitt_jeff_6],
        [self.RV_jitter_jeff_alpha_7,self.RV_jitter_jeff_beta_7,self.use_rvjitt_jeff_7], 
        [self.RV_jitter_jeff_alpha_8,self.RV_jitter_jeff_beta_8,self.use_rvjitt_jeff_8], 
        [self.RV_jitter_jeff_alpha_9,self.RV_jitter_jeff_beta_9,self.use_rvjitt_jeff_9], 
        [self.RV_jitter_jeff_alpha_10,self.RV_jitter_jeff_beta_10,self.use_rvjitt_jeff_10],
        [self.RV_jitter_jeff_alpha_11,self.RV_jitter_jeff_beta_11,self.use_rvjitt_jeff_11],
        [self.RV_jitter_jeff_alpha_12,self.RV_jitter_jeff_beta_12,self.use_rvjitt_jeff_12],
        [self.RV_jitter_jeff_alpha_13,self.RV_jitter_jeff_beta_13,self.use_rvjitt_jeff_13],
        [self.RV_jitter_jeff_alpha_14,self.RV_jitter_jeff_beta_14,self.use_rvjitt_jeff_14],
        [self.RV_jitter_jeff_alpha_15,self.RV_jitter_jeff_beta_15,self.use_rvjitt_jeff_15],
        [self.RV_jitter_jeff_alpha_16,self.RV_jitter_jeff_beta_16,self.use_rvjitt_jeff_16],
        [self.RV_jitter_jeff_alpha_17,self.RV_jitter_jeff_beta_17,self.use_rvjitt_jeff_17],
        [self.RV_jitter_jeff_alpha_19,self.RV_jitter_jeff_beta_18,self.use_rvjitt_jeff_18],
        [self.RV_jitter_jeff_alpha_19,self.RV_jitter_jeff_beta_19,self.use_rvjitt_jeff_19],
        [self.RV_jitter_jeff_alpha_20,self.RV_jitter_jeff_beta_20,self.use_rvjitt_jeff_20]
        ]
        
        return jitter_jeff_priors_gui




def offset_jeff_priors_gui_tra(self):


        offset_jeff_priors_gui_tra = [
        [self.tra_Data_alpha_1,self.tra_Data_beta_1,self.use_traoff_jeff_1], 
        [self.tra_Data_alpha_2,self.tra_Data_beta_2,self.use_traoff_jeff_2], 
        [self.tra_Data_alpha_3,self.tra_Data_beta_3,self.use_traoff_jeff_3], 
        [self.tra_Data_alpha_4,self.tra_Data_beta_4,self.use_traoff_jeff_4], 
        [self.tra_Data_alpha_5,self.tra_Data_beta_5,self.use_traoff_jeff_5], 
        [self.tra_Data_alpha_6,self.tra_Data_beta_6,self.use_traoff_jeff_6], 
        [self.tra_Data_alpha_7,self.tra_Data_beta_7,self.use_traoff_jeff_7], 
        [self.tra_Data_alpha_8,self.tra_Data_beta_8,self.use_traoff_jeff_8], 
        [self.tra_Data_alpha_9,self.tra_Data_beta_9,self.use_traoff_jeff_9], 
        [self.tra_Data_alpha_10,self.tra_Data_beta_10,self.use_traoff_jeff_10],
        [self.tra_Data_alpha_11,self.tra_Data_beta_11,self.use_traoff_jeff_11], 
        [self.tra_Data_alpha_12,self.tra_Data_beta_12,self.use_traoff_jeff_12], 
        [self.tra_Data_alpha_13,self.tra_Data_beta_13,self.use_traoff_jeff_13], 
        [self.tra_Data_alpha_14,self.tra_Data_beta_14,self.use_traoff_jeff_14], 
        [self.tra_Data_alpha_15,self.tra_Data_beta_15,self.use_traoff_jeff_15], 
        [self.tra_Data_alpha_16,self.tra_Data_beta_16,self.use_traoff_jeff_16], 
        [self.tra_Data_alpha_17,self.tra_Data_beta_17,self.use_traoff_jeff_17], 
        [self.tra_Data_alpha_18,self.tra_Data_beta_18,self.use_traoff_jeff_18], 
        [self.tra_Data_alpha_19,self.tra_Data_beta_19,self.use_traoff_jeff_19], 
        [self.tra_Data_alpha_20,self.tra_Data_beta_20,self.use_traoff_jeff_20]

        ]
        
        return offset_jeff_priors_gui_tra
        
    
    
def jitter_jeff_priors_gui_tra(self):

        jitter_jeff_priors_gui_tra = [
        [self.tra_jitter_alpha_1,self.tra_jitter_beta_1,self.use_trajitt_jeff_1], 
        [self.tra_jitter_alpha_2,self.tra_jitter_beta_2,self.use_trajitt_jeff_2], 
        [self.tra_jitter_alpha_3,self.tra_jitter_beta_3,self.use_trajitt_jeff_3], 
        [self.tra_jitter_alpha_4,self.tra_jitter_beta_4,self.use_trajitt_jeff_4], 
        [self.tra_jitter_alpha_5,self.tra_jitter_beta_5,self.use_trajitt_jeff_5], 
        [self.tra_jitter_alpha_6,self.tra_jitter_beta_6,self.use_trajitt_jeff_6],
        [self.tra_jitter_alpha_7,self.tra_jitter_beta_7,self.use_trajitt_jeff_7], 
        [self.tra_jitter_alpha_8,self.tra_jitter_beta_8,self.use_trajitt_jeff_8], 
        [self.tra_jitter_alpha_9,self.tra_jitter_beta_9,self.use_trajitt_jeff_9], 
        [self.tra_jitter_alpha_10,self.tra_jitter_beta_10,self.use_trajitt_jeff_10],
        [self.tra_jitter_alpha_11,self.tra_jitter_beta_11,self.use_trajitt_jeff_11], 
        [self.tra_jitter_alpha_12,self.tra_jitter_beta_12,self.use_trajitt_jeff_12], 
        [self.tra_jitter_alpha_13,self.tra_jitter_beta_13,self.use_trajitt_jeff_13], 
        [self.tra_jitter_alpha_14,self.tra_jitter_beta_14,self.use_trajitt_jeff_14], 
        [self.tra_jitter_alpha_15,self.tra_jitter_beta_15,self.use_trajitt_jeff_15], 
        [self.tra_jitter_alpha_16,self.tra_jitter_beta_16,self.use_trajitt_jeff_16],
        [self.tra_jitter_alpha_17,self.tra_jitter_beta_17,self.use_trajitt_jeff_17], 
        [self.tra_jitter_alpha_18,self.tra_jitter_beta_18,self.use_trajitt_jeff_18], 
        [self.tra_jitter_alpha_19,self.tra_jitter_beta_19,self.use_trajitt_jeff_19], 
        [self.tra_jitter_alpha_20,self.tra_jitter_beta_20,self.use_trajitt_jeff_20]

        ]
        
        return jitter_jeff_priors_gui_tra

 
def tra_lin_trend_jeff_priors_gui(self):

        tra_lin_trend_jeff_priors_gui = [
        [self.tra_lin_trend_alpha_1,self.tra_lin_trend_beta_1,self.use_tra_lin_jeff_1], 
        [self.tra_lin_trend_alpha_2,self.tra_lin_trend_beta_2,self.use_tra_lin_jeff_2], 
        [self.tra_lin_trend_alpha_3,self.tra_lin_trend_beta_3,self.use_tra_lin_jeff_3], 
        [self.tra_lin_trend_alpha_4,self.tra_lin_trend_beta_4,self.use_tra_lin_jeff_4], 
        [self.tra_lin_trend_alpha_5,self.tra_lin_trend_beta_5,self.use_tra_lin_jeff_5], 
        [self.tra_lin_trend_alpha_6,self.tra_lin_trend_beta_6,self.use_tra_lin_jeff_6],
        [self.tra_lin_trend_alpha_7,self.tra_lin_trend_beta_7,self.use_tra_lin_jeff_7], 
        [self.tra_lin_trend_alpha_8,self.tra_lin_trend_beta_8,self.use_tra_lin_jeff_8], 
        [self.tra_lin_trend_alpha_9,self.tra_lin_trend_beta_9,self.use_tra_lin_jeff_9], 
        [self.tra_lin_trend_alpha_10,self.tra_lin_trend_beta_10,self.use_tra_lin_jeff_10],
        [self.tra_lin_trend_alpha_11,self.tra_lin_trend_beta_11,self.use_tra_lin_jeff_11], 
        [self.tra_lin_trend_alpha_12,self.tra_lin_trend_beta_12,self.use_tra_lin_jeff_12], 
        [self.tra_lin_trend_alpha_13,self.tra_lin_trend_beta_13,self.use_tra_lin_jeff_13], 
        [self.tra_lin_trend_alpha_14,self.tra_lin_trend_beta_14,self.use_tra_lin_jeff_14], 
        [self.tra_lin_trend_alpha_15,self.tra_lin_trend_beta_15,self.use_tra_lin_jeff_15], 
        [self.tra_lin_trend_alpha_16,self.tra_lin_trend_beta_16,self.use_tra_lin_jeff_16],
        [self.tra_lin_trend_alpha_17,self.tra_lin_trend_beta_17,self.use_tra_lin_jeff_17], 
        [self.tra_lin_trend_alpha_18,self.tra_lin_trend_beta_18,self.use_tra_lin_jeff_18], 
        [self.tra_lin_trend_alpha_19,self.tra_lin_trend_beta_19,self.use_tra_lin_jeff_19], 
        [self.tra_lin_trend_alpha_20,self.tra_lin_trend_beta_20,self.use_tra_lin_jeff_20]
        ]
        
        return tra_lin_trend_jeff_priors_gui
    
    
def tra_quad_trend_jeff_priors_gui(self):

        tra_quad_trend_jeff_priors_gui = [
        [self.tra_quad_trend_alpha_1,self.tra_quad_trend_beta_1,self.use_tra_quad_jeff_1], 
        [self.tra_quad_trend_alpha_2,self.tra_quad_trend_beta_2,self.use_tra_quad_jeff_2], 
        [self.tra_quad_trend_alpha_3,self.tra_quad_trend_beta_3,self.use_tra_quad_jeff_3], 
        [self.tra_quad_trend_alpha_4,self.tra_quad_trend_beta_4,self.use_tra_quad_jeff_4], 
        [self.tra_quad_trend_alpha_5,self.tra_quad_trend_beta_5,self.use_tra_quad_jeff_5], 
        [self.tra_quad_trend_alpha_6,self.tra_quad_trend_beta_6,self.use_tra_quad_jeff_6],
        [self.tra_quad_trend_alpha_7,self.tra_quad_trend_beta_7,self.use_tra_quad_jeff_7], 
        [self.tra_quad_trend_alpha_8,self.tra_quad_trend_beta_8,self.use_tra_quad_jeff_8], 
        [self.tra_quad_trend_alpha_9,self.tra_quad_trend_beta_9,self.use_tra_quad_jeff_9], 
        [self.tra_quad_trend_alpha_10,self.tra_quad_trend_beta_10,self.use_tra_quad_jeff_10],
        [self.tra_quad_trend_alpha_11,self.tra_quad_trend_beta_11,self.use_tra_quad_jeff_11], 
        [self.tra_quad_trend_alpha_12,self.tra_quad_trend_beta_12,self.use_tra_quad_jeff_12], 
        [self.tra_quad_trend_alpha_13,self.tra_quad_trend_beta_13,self.use_tra_quad_jeff_13], 
        [self.tra_quad_trend_alpha_14,self.tra_quad_trend_beta_14,self.use_tra_quad_jeff_14], 
        [self.tra_quad_trend_alpha_15,self.tra_quad_trend_beta_15,self.use_tra_quad_jeff_15], 
        [self.tra_quad_trend_alpha_16,self.tra_quad_trend_beta_16,self.use_tra_quad_jeff_16],
        [self.tra_quad_trend_alpha_17,self.tra_quad_trend_beta_17,self.use_tra_quad_jeff_17], 
        [self.tra_quad_trend_alpha_18,self.tra_quad_trend_beta_18,self.use_tra_quad_jeff_18], 
        [self.tra_quad_trend_alpha_19,self.tra_quad_trend_beta_19,self.use_tra_quad_jeff_19], 
        [self.tra_quad_trend_alpha_20,self.tra_quad_trend_beta_20,self.use_tra_quad_jeff_20]
        ]
        
        return tra_quad_trend_jeff_priors_gui
        
    
 
 
################### RV GP ########################

def gp_rot_params(self):

    gp_rot_params = [
            self.GP_rot_kernel_Amp,
            self.GP_rot_kernel_time_sc,
            self.GP_rot_kernel_Per,
            self.GP_rot_kernel_fact
            ]

    return gp_rot_params


def use_gp_rot_params(self):
    use_gp_rot_params = [
            self.use_GP_rot_kernel_Amp,
            self.use_GP_rot_kernel_time_sc,
            self.use_GP_rot_kernel_Per,
            self.use_GP_rot_kernel_fact
            ]
    return use_gp_rot_params

def gp_rot_errors_gui(self):

    gp_rot_errors_gui = [
            self.err_rot_kernel_Amp,
            self.err_rot_kernel_time_sc,
            self.err_rot_kernel_Per,
            self.err_rot_kernel_fact
            ]
        
    return gp_rot_errors_gui



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


def gp_double_sho_params(self):

    gp_double_sho_params = [
            self.GP_double_sho_kernel_sigma,
            self.GP_double_sho_kernel_P,
            self.GP_double_sho_kernel_Q0,
            self.GP_double_sho_kernel_dQ,
            self.GP_double_sho_kernel_f,
            ]

    return gp_double_sho_params


def use_gp_double_sho_params(self):

    use_gp_double_sho_params = [
            self.use_GP_double_sho_kernel_sigma,
            self.use_GP_double_sho_kernel_P,
            self.use_GP_double_sho_kernel_Q0,
            self.use_GP_double_sho_kernel_dQ,
            self.use_GP_double_sho_kernel_f,
            ]
    return use_gp_double_sho_params


def gp_double_sho_errors_gui(self):

    gp_double_sho_errors_gui = [
            self.err_GP_double_sho_kernel_sigma,
            self.err_GP_double_sho_kernel_P,
            self.err_GP_double_sho_kernel_Q0,
            self.err_GP_double_sho_kernel_dQ,
            self.err_GP_double_sho_kernel_f,
            ]
    
    return gp_double_sho_errors_gui


##################### Matern 3/2 ##############################################

def gp_mat_params(self):

    gp_mat_params = [
            self.GP_mat_kernel_sigma,
            self.GP_mat_kernel_rho,
            self.GP_mat_kernel_eps
            ]

    return gp_mat_params




def use_gp_mat_params(self):

    use_gp_mat_params = [
            self.use_GP_mat_kernel_sigma,
            self.use_GP_mat_kernel_rho,
            self.use_GP_mat_kernel_eps
            ]

    return use_gp_mat_params


def gp_mat_errors_gui(self):

    gp_mat_errors_gui = [
            self.err_mat_kernel_sigma,
            self.err_mat_kernel_rho,
            self.err_mat_kernel_eps
            ]
    
    return gp_mat_errors_gui



##################### RV Real term DRW ##############################################

def gp_drw_params(self):

    gp_drw_params = [
            self.GP_drw_kernel_a,
            self.GP_drw_kernel_c
            ]

    return gp_drw_params




def use_gp_drw_params(self):

    use_gp_drw_params = [
            self.use_GP_drw_kernel_a,
            self.use_GP_drw_kernel_c
            ]

    return use_gp_drw_params


def gp_drw_errors_gui(self):

    gp_drw_errors_gui = [
            self.err_drw_kernel_a,
            self.err_drw_kernel_c
            ]
    
    return gp_drw_errors_gui


################### Tra GP ########################


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


def tra_gp_rot_errors_gui(self):

    tra_gp_rot_errors_gui = [
            self.err_tra_rot_kernel_Amp,
            self.err_tra_rot_kernel_time_sc,
            self.err_tra_rot_kernel_Per,
            self.err_tra_rot_kernel_fact
            ]
        
    return tra_gp_rot_errors_gui


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


def gp_tra_sho_errors_gui(self):

    tra_gp_sho_errors_gui = [
            self.err_tra_sho_kernel_S,
            self.err_tra_sho_kernel_Q,
            self.err_tra_sho_kernel_omega
            ]
    
    return tra_gp_sho_errors_gui


def tra_gp_double_sho_params(self):

    tra_gp_double_sho_params = [
            self.tra_GP_double_sho_kernel_sigma,
            self.tra_GP_double_sho_kernel_P,
            self.tra_GP_double_sho_kernel_Q0,
            self.tra_GP_double_sho_kernel_dQ,
            self.tra_GP_double_sho_kernel_f,
            ]

    return tra_gp_double_sho_params


def use_tra_gp_double_sho_params(self):

    use_tra_gp_double_sho_params = [
            self.use_tra_GP_double_sho_kernel_sigma,
            self.use_tra_GP_double_sho_kernel_P,
            self.use_tra_GP_double_sho_kernel_Q0,
            self.use_tra_GP_double_sho_kernel_dQ,
            self.use_tra_GP_double_sho_kernel_f,
            ]
    return use_tra_gp_double_sho_params


def tra_gp_double_sho_errors_gui(self):

    tra_gp_double_sho_errors_gui = [
            self.err_tra_GP_double_sho_kernel_sigma,
            self.err_tra_GP_double_sho_kernel_P,
            self.err_tra_GP_double_sho_kernel_Q0,
            self.err_tra_GP_double_sho_kernel_dQ,
            self.err_tra_GP_double_sho_kernel_f,
            ]
    
    return tra_gp_double_sho_errors_gui


##################### Tra Matern 3/2 ##############################################


def tra_gp_mat_params(self):

    tra_gp_mat_params = [
            self.tra_GP_mat_kernel_sigma,
            self.tra_GP_mat_kernel_rho,
            self.tra_GP_mat_kernel_eps
            ]

    return tra_gp_mat_params

def use_tra_gp_mat_params(self):

    use_tra_gp_mat_params = [
            self.use_tra_GP_mat_kernel_sigma,
            self.use_tra_GP_mat_kernel_rho,
            self.use_tra_GP_mat_kernel_eps
            ]

    return use_tra_gp_mat_params


def tra_gp_mat_errors_gui(self):

    tra_gp_mat_errors_gui = [
            self.err_tra_mat_kernel_sigma,
            self.err_tra_mat_kernel_rho,
            self.err_tra_mat_kernel_eps
            ]
    
    return tra_gp_mat_errors_gui



##################### Tra RealTerm DRW ##############################################


def tra_gp_drw_params(self):

    tra_gp_drw_params = [
            self.tra_GP_drw_kernel_a,
            self.tra_GP_drw_kernel_c
            ]

    return tra_gp_drw_params

def use_tra_gp_drw_params(self):

    use_tra_gp_drw_params = [
            self.use_tra_GP_drw_kernel_a,
            self.use_tra_GP_drw_kernel_c
            ]

    return use_tra_gp_drw_params


def tra_gp_drw_errors_gui(self):

    tra_gp_drw_errors_gui = [
            self.tra_err_drw_kernel_a,
            self.tra_err_drw_kernel_c
            ]
    
    return tra_gp_drw_errors_gui



############### GP priors ######################



def GP_sho_bounds_gui(self):
    
    GP_sho_bounds_gui = [
        [self.GP_sho_kernel_S_min,self.GP_sho_kernel_S_max,self.use_GP_sho_kernel_S_bound],  
        [self.GP_sho_kernel_Q_min,self.GP_sho_kernel_Q_max,self.use_GP_sho_kernel_Q_bound],  
        [self.GP_sho_kernel_omega_min,self.GP_sho_kernel_omega_max,self.use_GP_sho_kernel_omega_bound],  
        ]
        
    return GP_sho_bounds_gui


def GP_double_sho_bounds_gui(self):
    
    GP_double_sho_bounds_gui = [
        [self.GP_double_sho_kernel_sigma_min,self.GP_double_sho_kernel_sigma_max,self.use_GP_double_sho_kernel_sigma_bound],  
        [self.GP_double_sho_kernel_P_min,self.GP_double_sho_kernel_P_max,self.use_GP_double_sho_kernel_P_bound],  
        [self.GP_double_sho_kernel_Q0_min,self.GP_double_sho_kernel_Q0_max,self.use_GP_double_sho_kernel_Q0_bound],  
        [self.GP_double_sho_kernel_dQ_min,self.GP_double_sho_kernel_dQ_max,self.use_GP_double_sho_kernel_dQ_bound],  
        [self.GP_double_sho_kernel_f_min,self.GP_double_sho_kernel_f_max,self.use_GP_double_sho_kernel_f_bound] 
        ]
        
    return GP_double_sho_bounds_gui


def GP_rot_bounds_gui(self):

    GP_rot_bounds_gui = [
        [self.GP_rot_kernel_Amp_min,self.GP_rot_kernel_Amp_max,self.use_GP_rot_kernel_Amp_bound],  
        [self.GP_rot_kernel_time_sc_min,self.GP_rot_kernel_time_sc_max,self.use_GP_rot_kernel_time_sc_bound],  
        [self.GP_rot_kernel_Per_min,self.GP_rot_kernel_Per_max,self.use_GP_rot_kernel_Per_sigma_bound],  
        [self.GP_rot_kernel_fact_min,self.GP_rot_kernel_fact_max,self.use_GP_rot_kernel_fact_bound],  
        ]
        
    return GP_rot_bounds_gui

def GP_mat_bounds_gui(self):

    GP_mat_bounds_gui = [
        [self.GP_mat_kernel_sigma_min,self.GP_mat_kernel_sigma_max,self.use_GP_mat_kernel_sigma_bound],  
        [self.GP_mat_kernel_rho_min,self.GP_mat_kernel_rho_max,self.use_GP_mat_kernel_rho_bound],  
        [self.GP_mat_kernel_eps_min,self.GP_mat_kernel_eps_max,self.use_GP_mat_kernel_sigma_bound],  
        ]
    return GP_mat_bounds_gui


def GP_drw_bounds_gui(self):

    GP_drw_bounds_gui = [
        [self.GP_drw_kernel_a_min,self.GP_drw_kernel_a_max,self.use_GP_drw_kernel_a_bound],  
        [self.GP_drw_kernel_c_min,self.GP_drw_kernel_c_max,self.use_GP_drw_kernel_c_bound]  
        ]
    return GP_drw_bounds_gui



def GP_rot_nr_priors_gui(self):

    GP_rot_nr_priors_gui = [
        [self.GP_rot_kernel_Amp_mean,self.GP_rot_kernel_Amp_sigma,self.use_GP_rot_kernel_Amp_nr_pr],  
        [self.GP_rot_kernel_time_sc_mean,self.GP_rot_kernel_time_sc_sigma,self.use_GP_rot_kernel_time_sc_nr_pr],  
        [self.GP_rot_kernel_Per_mean,self.GP_rot_kernel_Per_sigma,self.use_GP_rot_kernel_Per_sigma_nr_pr],  
        [self.GP_rot_kernel_fact_mean,self.GP_rot_kernel_fact_sigma,self.use_GP_rot_kernel_fact_nr_pr],  
        ]

    return GP_rot_nr_priors_gui

def GP_sho_nr_priors_gui(self):

    GP_sho_nr_priors_gui = [
        [self.GP_sho_kernel_S_mean,self.GP_sho_kernel_S_sigma, self.use_GP_sho_kernel_S_nr_pr],  
        [self.GP_sho_kernel_Q_mean,self.GP_sho_kernel_Q_sigma, self.use_GP_sho_kernel_Q_nr_pr],  
        [self.GP_sho_kernel_omega_mean,self.GP_sho_kernel_omega_sigma, self.use_GP_sho_kernel_omega_nr_pr],  
        ]
    
    return GP_sho_nr_priors_gui
  

def GP_double_sho_nr_priors_gui(self):

 
    GP_double_sho_nr_priors_gui = [
        [self.GP_double_sho_kernel_sigma_mean,self.GP_double_sho_kernel_sigma_sigma,self.use_GP_double_sho_kernel_sigma_nr_pr],  
        [self.GP_double_sho_kernel_P_mean,self.GP_double_sho_kernel_P_sigma,self.use_GP_double_sho_kernel_P_nr_pr],  
        [self.GP_double_sho_kernel_Q0_mean,self.GP_double_sho_kernel_Q0_sigma,self.use_GP_double_sho_kernel_Q0_nr_pr],  
        [self.GP_double_sho_kernel_dQ_mean,self.GP_double_sho_kernel_dQ_sigma,self.use_GP_double_sho_kernel_dQ_nr_pr],  
        [self.GP_double_sho_kernel_f_mean,self.GP_double_sho_kernel_f_sigma,self.use_GP_double_sho_kernel_f_nr_pr] 
        ]    
    return GP_double_sho_nr_priors_gui


def GP_mat_nr_priors_gui(self):

    GP_mat_nr_priors_gui = [
        [self.GP_mat_kernel_sigma_mean,self.GP_mat_kernel_sigma_sigma, self.use_GP_mat_kernel_sigma_nr_pr],  
        [self.GP_mat_kernel_rho_mean,self.GP_mat_kernel_rho_sigma, self.use_GP_mat_kernel_rho_nr_pr],  
        [self.GP_mat_kernel_eps_mean,self.GP_mat_kernel_eps_sigma, self.use_GP_mat_kernel_eps_nr_pr],  
        ]
        
    return GP_mat_nr_priors_gui


def GP_drw_nr_priors_gui(self):

    GP_drw_nr_priors_gui = [
        [self.GP_drw_kernel_a_mean,self.GP_drw_kernel_a_sigma,self.use_GP_drw_kernel_a_nr_pr],  
        [self.GP_drw_kernel_c_mean,self.GP_drw_kernel_c_sigma,self.use_GP_drw_kernel_c_nr_pr],  
        ]
    return GP_drw_nr_priors_gui






def tra_GP_sho_bounds_gui(self):
    
    tra_GP_sho_bounds_gui = [
        [self.tra_GP_sho_kernel_S_min,self.tra_GP_sho_kernel_S_max,self.use_tra_GP_sho_kernel_S_bound],  
        [self.tra_GP_sho_kernel_Q_min,self.tra_GP_sho_kernel_Q_max,self.use_tra_GP_sho_kernel_Q_bound],  
        [self.tra_GP_sho_kernel_omega_min,self.tra_GP_sho_kernel_omega_max,self.use_tra_GP_sho_kernel_omega_bound],  
        ]
        
    return tra_GP_sho_bounds_gui

def tra_GP_double_sho_bounds_gui(self):
    
    tra_GP_double_sho_bounds_gui = [
        [self.tra_GP_double_sho_kernel_sigma_min,self.tra_GP_double_sho_kernel_sigma_max,self.use_tra_GP_double_sho_kernel_sigma_bound],  
        [self.tra_GP_double_sho_kernel_P_min,self.tra_GP_double_sho_kernel_P_max,self.use_tra_GP_double_sho_kernel_P_bound],  
        [self.tra_GP_double_sho_kernel_Q0_min,self.tra_GP_double_sho_kernel_Q0_max,self.use_tra_GP_double_sho_kernel_Q0_bound],  
        [self.tra_GP_double_sho_kernel_dQ_min,self.tra_GP_double_sho_kernel_dQ_max,self.use_tra_GP_double_sho_kernel_dQ_bound],  
        [self.tra_GP_double_sho_kernel_f_min,self.tra_GP_double_sho_kernel_f_max,self.use_tra_GP_double_sho_kernel_f_bound] 
        ]
        
    return tra_GP_double_sho_bounds_gui

def tra_GP_rot_bounds_gui(self):

    tra_GP_rot_bounds_gui = [
        [self.tra_GP_rot_kernel_Amp_min,self.tra_GP_rot_kernel_Amp_max,self.use_tra_GP_rot_kernel_Amp_bound],  
        [self.tra_GP_rot_kernel_time_sc_min,self.tra_GP_rot_kernel_time_sc_max,self.use_tra_GP_rot_kernel_time_sc_bound],  
        [self.tra_GP_rot_kernel_Per_min,self.tra_GP_rot_kernel_Per_max,self.use_tra_GP_rot_kernel_Per_sigma_bound],  
        [self.tra_GP_rot_kernel_fact_min,self.tra_GP_rot_kernel_fact_max,self.use_tra_GP_rot_kernel_fact_bound],  
        ]
        
    return tra_GP_rot_bounds_gui

def tra_GP_mat_bounds_gui(self):

    tra_GP_mat_bounds_gui = [
        [self.tra_GP_mat_kernel_sigma_min,self.tra_GP_mat_kernel_sigma_max,self.use_tra_GP_mat_kernel_sigma_bound],  
        [self.tra_GP_mat_kernel_rho_min,self.tra_GP_mat_kernel_rho_max,self.use_tra_GP_mat_kernel_rho_bound],  
        [self.tra_GP_mat_kernel_eps_min,self.tra_GP_mat_kernel_eps_max,self.use_tra_GP_mat_kernel_sigma_bound],  
        ]
    return tra_GP_mat_bounds_gui
 

def tra_GP_drw_bounds_gui(self):

    tra_GP_drw_bounds_gui = [
        [self.tra_GP_drw_kernel_a_min,self.tra_GP_drw_kernel_a_max,self.use_tra_GP_drw_kernel_a_bound],  
        [self.tra_GP_drw_kernel_c_min,self.tra_GP_drw_kernel_c_max,self.use_tra_GP_drw_kernel_c_bound]  
        ]
    return tra_GP_drw_bounds_gui
 
    

def tra_GP_rot_nr_priors_gui(self):

    tra_GP_rot_nr_priors_gui = [
        [self.tra_GP_rot_kernel_Amp_mean,self.tra_GP_rot_kernel_Amp_sigma,self.use_tra_GP_rot_kernel_Amp_nr_pr],  
        [self.tra_GP_rot_kernel_time_sc_mean,self.tra_GP_rot_kernel_time_sc_sigma,self.use_tra_GP_rot_kernel_time_sc_nr_pr],  
        [self.tra_GP_rot_kernel_Per_mean,self.tra_GP_rot_kernel_Per_sigma,self.use_tra_GP_rot_kernel_Per_sigma_nr_pr],  
        [self.tra_GP_rot_kernel_fact_mean,self.tra_GP_rot_kernel_fact_sigma,self.use_tra_GP_rot_kernel_fact_nr_pr],  
        ]
 
    return tra_GP_rot_nr_priors_gui

def tra_GP_sho_nr_priors_gui(self):
    
    tra_GP_sho_nr_priors_gui = [
        [self.tra_GP_sho_kernel_S_mean,self.tra_GP_sho_kernel_S_sigma, self.use_tra_GP_sho_kernel_S_nr_pr],  
        [self.tra_GP_sho_kernel_Q_mean,self.tra_GP_sho_kernel_Q_sigma, self.use_tra_GP_sho_kernel_Q_nr_pr],  
        [self.tra_GP_sho_kernel_omega_mean,self.tra_GP_sho_kernel_omega_sigma, self.use_tra_GP_sho_kernel_omega_nr_pr],  
        ]
 
    return tra_GP_sho_nr_priors_gui

def tra_GP_double_sho_nr_priors_gui(self):

 
    tra_GP_double_sho_nr_priors_gui = [
        [self.tra_GP_double_sho_kernel_sigma_mean,self.tra_GP_double_sho_kernel_sigma_sigma,self.use_tra_GP_double_sho_kernel_sigma_nr_pr],  
        [self.tra_GP_double_sho_kernel_P_mean,self.tra_GP_double_sho_kernel_P_sigma,self.use_tra_GP_double_sho_kernel_P_nr_pr],  
        [self.tra_GP_double_sho_kernel_Q0_mean,self.tra_GP_double_sho_kernel_Q0_sigma,self.use_tra_GP_double_sho_kernel_Q0_nr_pr],  
        [self.tra_GP_double_sho_kernel_dQ_mean,self.tra_GP_double_sho_kernel_dQ_sigma,self.use_tra_GP_double_sho_kernel_dQ_nr_pr],  
        [self.tra_GP_double_sho_kernel_f_mean,self.tra_GP_double_sho_kernel_f_sigma,self.use_tra_GP_double_sho_kernel_f_nr_pr] 
        ]    
    return tra_GP_double_sho_nr_priors_gui



def tra_GP_mat_nr_priors_gui(self):
            
    tra_GP_mat_nr_priors_gui = [
        [self.tra_GP_mat_kernel_sigma_mean,self.tra_GP_mat_kernel_sigma_sigma, self.use_tra_GP_mat_kernel_sigma_nr_pr],  
        [self.tra_GP_mat_kernel_rho_mean,self.tra_GP_mat_kernel_rho_sigma, self.use_tra_GP_mat_kernel_rho_nr_pr],  
        [self.tra_GP_mat_kernel_eps_mean,self.tra_GP_mat_kernel_eps_sigma, self.use_tra_GP_mat_kernel_eps_nr_pr],  
        ]
    return tra_GP_mat_nr_priors_gui




def tra_GP_drw_nr_priors_gui(self):

    tra_GP_drw_nr_priors_gui = [
        [self.tra_GP_drw_kernel_a_mean,self.tra_GP_drw_kernel_a_sigma,self.use_tra_GP_drw_kernel_a_nr_pr], 
        [self.tra_GP_drw_kernel_c_mean,self.tra_GP_drw_kernel_c_sigma,self.use_tra_GP_drw_kernel_c_nr_pr],  
        ]
    return tra_GP_drw_nr_priors_gui







def GP_rot_jeff_priors_gui(self):

    GP_rot_jeff_priors_gui = [
        [self.GP_rot_kernel_Amp_jeff_alpha,self.GP_rot_kernel_Amp_jeff_beta,self.use_GP_rot_kernel_Amp_jeff_pr],  
        [self.GP_rot_kernel_time_sc_jeff_alpha,self.GP_rot_kernel_time_sc_jeff_beta,self.use_GP_rot_kernel_time_sc_jeff_pr],  
        [self.GP_rot_kernel_Per_jeff_alpha,self.GP_rot_kernel_Per_jeff_beta,self.use_GP_rot_kernel_Per_jeff_pr],  
        [self.GP_rot_kernel_fact_jeff_alpha,self.GP_rot_kernel_fact_jeff_beta,self.use_GP_rot_kernel_fact_jeff_pr],  
        ]
 
    return GP_rot_jeff_priors_gui


              
def GP_sho_jeff_priors_gui(self):

    GP_sho_jeff_priors_gui = [
        [self.GP_sho_kernel_S_jeff_alpha,self.GP_sho_kernel_S_jeff_beta, self.use_GP_sho_kernel_S_jeff_pr],  
        [self.GP_sho_kernel_Q_jeff_alpha,self.GP_sho_kernel_Q_jeff_beta, self.use_GP_sho_kernel_Q_jeff_pr],  
        [self.GP_sho_kernel_omega_jeff_alpha,self.GP_sho_kernel_omega_jeff_beta, self.use_GP_sho_kernel_omega_jeff_pr],  
        ]
 
    return GP_sho_jeff_priors_gui


def GP_double_sho_jeff_priors_gui(self):
    
    GP_double_sho_jeff_priors_gui = [
        [self.GP_double_sho_kernel_sigma_alpha,self.GP_double_sho_kernel_sigma_beta,self.use_GP_double_sho_kernel_sigma_jeff_pr],  
        [self.GP_double_sho_kernel_P_alpha,self.GP_double_sho_kernel_P_beta,self.use_GP_double_sho_kernel_P_jeff_pr],  
        [self.GP_double_sho_kernel_Q0_alpha,self.GP_double_sho_kernel_Q0_beta,self.use_GP_double_sho_kernel_Q0_jeff_pr],  
        [self.GP_double_sho_kernel_dQ_alpha,self.GP_double_sho_kernel_dQ_beta,self.use_GP_double_sho_kernel_dQ_jeff_pr],  
        [self.GP_double_sho_kernel_f_alpha,self.GP_double_sho_kernel_f_beta,self.use_GP_double_sho_kernel_f_jeff_pr] 
        ]
        
    return GP_double_sho_jeff_priors_gui

 
def GP_mat_jeff_priors_gui(self):
            
    GP_mat_jeff_priors_gui = [
        [self.GP_mat_kernel_sigma_alpha,self.GP_mat_kernel_sigma_beta, self.use_GP_mat_kernel_sigma_jeff_pr],  
        [self.GP_mat_kernel_rho_alpha,self.GP_mat_kernel_rho_beta, self.use_GP_mat_kernel_rho_jeff_pr],  
        [self.GP_mat_kernel_eps_alpha,self.GP_mat_kernel_eps_beta, self.use_GP_mat_kernel_eps_jeff_pr],  
        ]
        
    return GP_mat_jeff_priors_gui



def GP_drw_jeff_priors_gui(self):
            
    GP_drw_jeff_priors_gui = [
        [self.GP_drw_kernel_a_alpha,self.GP_drw_kernel_a_beta, self.use_GP_drw_kernel_a_jeff_pr],  
        [self.GP_drw_kernel_c_alpha,self.GP_drw_kernel_c_beta, self.use_GP_drw_kernel_c_jeff_pr]   
        ]
        
    return GP_drw_jeff_priors_gui






def tra_GP_rot_jeff_priors_gui(self):

    tra_GP_rot_jeff_priors_gui = [
        [self.tra_GP_rot_kernel_Amp_jeff_alpha,self.tra_GP_rot_kernel_Amp_jeff_beta,self.use_tra_GP_rot_kernel_Amp_jeff_pr],  
        [self.tra_GP_rot_kernel_time_sc_jeff_alpha,self.tra_GP_rot_kernel_time_sc_jeff_beta,self.use_tra_GP_rot_kernel_time_sc_jeff_pr],  
        [self.tra_GP_rot_kernel_Per_jeff_alpha,self.tra_GP_rot_kernel_Per_jeff_beta,self.use_tra_GP_rot_kernel_Per_jeff_pr],  
        [self.tra_GP_rot_kernel_fact_jeff_alpha,self.tra_GP_rot_kernel_fact_jeff_beta,self.use_tra_GP_rot_kernel_fact_jeff_pr],  
        ]
 
    return tra_GP_rot_jeff_priors_gui
              
def tra_GP_sho_jeff_priors_gui(self):

    tra_GP_sho_jeff_priors_gui = [
        [self.tra_GP_sho_kernel_S_jeff_alpha,self.tra_GP_sho_kernel_S_jeff_beta, self.use_tra_GP_sho_kernel_S_jeff_pr],  
        [self.tra_GP_sho_kernel_Q_jeff_alpha,self.tra_GP_sho_kernel_Q_jeff_beta, self.use_tra_GP_sho_kernel_Q_jeff_pr],  
        [self.tra_GP_sho_kernel_omega_jeff_alpha,self.tra_GP_sho_kernel_omega_jeff_beta, self.use_tra_GP_sho_kernel_omega_jeff_pr],  
        ]
 
    return tra_GP_sho_jeff_priors_gui
 

def tra_GP_double_sho_jeff_priors_gui(self):
    
    tra_GP_double_sho_jeff_priors_gui = [
        [self.tra_GP_double_sho_kernel_sigma_alpha,self.tra_GP_double_sho_kernel_sigma_beta,self.use_tra_GP_double_sho_kernel_sigma_jeff_pr],  
        [self.tra_GP_double_sho_kernel_P_alpha,self.tra_GP_double_sho_kernel_P_beta,self.use_tra_GP_double_sho_kernel_P_jeff_pr],  
        [self.tra_GP_double_sho_kernel_Q0_alpha,self.tra_GP_double_sho_kernel_Q0_beta,self.use_tra_GP_double_sho_kernel_Q0_jeff_pr],  
        [self.tra_GP_double_sho_kernel_dQ_alpha,self.tra_GP_double_sho_kernel_dQ_beta,self.use_tra_GP_double_sho_kernel_dQ_jeff_pr],  
        [self.tra_GP_double_sho_kernel_f_alpha,self.tra_GP_double_sho_kernel_f_beta,self.use_tra_GP_double_sho_kernel_f_jeff_pr] 
        ]
        
    return tra_GP_double_sho_jeff_priors_gui



def tra_GP_mat_jeff_priors_gui(self):
            
    tra_GP_mat_jeff_priors_gui = [
        [self.tra_GP_mat_kernel_sigma_alpha,self.tra_GP_mat_kernel_sigma_beta, self.use_tra_GP_mat_kernel_sigma_jeff_pr],  
        [self.tra_GP_mat_kernel_rho_alpha,self.tra_GP_mat_kernel_rho_beta, self.use_tra_GP_mat_kernel_rho_jeff_pr],  
        [self.tra_GP_mat_kernel_eps_alpha,self.tra_GP_mat_kernel_eps_beta, self.use_tra_GP_mat_kernel_eps_jeff_pr],  
        ]
        
    return tra_GP_mat_jeff_priors_gui


def tra_GP_drw_jeff_priors_gui(self):
            
    tra_GP_drw_jeff_priors_gui = [
        [self.tra_GP_drw_kernel_a_alpha,self.tra_GP_drw_kernel_a_beta, self.use_tra_GP_drw_kernel_a_jeff_pr],  
        [self.tra_GP_drw_kernel_c_alpha,self.tra_GP_drw_kernel_c_beta, self.use_tra_GP_drw_kernel_c_jeff_pr]   
        ]
        
    return tra_GP_drw_jeff_priors_gui






def use_tra_data_GP(self):
            
    use_tra_data_GP = [
        self.use_tra_gp_1,self.use_tra_gp_2, self.use_tra_gp_3,  
        self.use_tra_gp_4,self.use_tra_gp_5, self.use_tra_gp_6,  
        self.use_tra_gp_8,self.use_tra_gp_7, self.use_tra_gp_9,
        self.use_tra_gp_10,
        self.use_tra_gp_11,self.use_tra_gp_12, self.use_tra_gp_13,  
        self.use_tra_gp_14,self.use_tra_gp_15, self.use_tra_gp_16,  
        self.use_tra_gp_18,self.use_tra_gp_17, self.use_tra_gp_19,
        self.use_tra_gp_20,
        ]
        
    return use_tra_data_GP



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
 

def act_sigma_clip(self): 
    act_sigma_clip = [
            [self.act_sigma_clip_1,self.use_act_sigma_clip_1],[self.act_sigma_clip_2,self.use_act_sigma_clip_2],
            [self.act_sigma_clip_3,self.use_act_sigma_clip_3],[self.act_sigma_clip_4,self.use_act_sigma_clip_4],
            [self.act_sigma_clip_5,self.use_act_sigma_clip_5],[self.act_sigma_clip_6,self.use_act_sigma_clip_6],
            [self.act_sigma_clip_7,self.use_act_sigma_clip_7],[self.act_sigma_clip_8,self.use_act_sigma_clip_8],
            [self.act_sigma_clip_9,self.use_act_sigma_clip_9],[self.act_sigma_clip_10,self.use_act_sigma_clip_10]
            ]
    return act_sigma_clip

def act_remove_mean(self): 
    act_remove_mean = [
             self.act_remove_mean_1,self.act_remove_mean_2,self.act_remove_mean_3,
             self.act_remove_mean_4,self.act_remove_mean_5,self.act_remove_mean_6,
             self.act_remove_mean_7,self.act_remove_mean_8,self.act_remove_mean_9,
             self.act_remove_mean_10
            ]
    return act_remove_mean


#def tra_sigma_clip(self): 
#    tra_sigma_clip = [
#            [self.tra_sigma_clip_1,self.use_tra_sigma_clip_1],[self.tra_sigma_clip_2,self.use_tra_sigma_clip_2],
#            [self.tra_sigma_clip_3,self.use_tra_sigma_clip_3],[self.tra_sigma_clip_4,self.use_tra_sigma_clip_4],
#            [self.tra_sigma_clip_5,self.use_tra_sigma_clip_5],[self.tra_sigma_clip_6,self.use_tra_sigma_clip_6],
#            [self.tra_sigma_clip_7,self.use_tra_sigma_clip_7],[self.tra_sigma_clip_8,self.use_tra_sigma_clip_8],
#            [self.tra_sigma_clip_9,self.use_tra_sigma_clip_9],[self.tra_sigma_clip_10,self.use_tra_sigma_clip_10]
#            ]
 
#    return tra_sigma_clip

def tra_norm(self): 
    tra_norm = [
             self.tra_norm_1,self.tra_norm_2,self.tra_norm_3,self.tra_norm_4,self.tra_norm_5,
             self.tra_norm_6,self.tra_norm_7,self.tra_norm_8,self.tra_norm_9,self.tra_norm_10,
             self.tra_norm_11,self.tra_norm_12,self.tra_norm_13,self.tra_norm_14,self.tra_norm_15,
             self.tra_norm_16,self.tra_norm_17,self.tra_norm_18,self.tra_norm_19,self.tra_norm_20
            ]
 
    return tra_norm



def tra_dilution(self): 
    tra_dilution = [
             [self.tra_dilution_1,self.use_tra_dilution_1], [self.tra_dilution_2,self.use_tra_dilution_2],
             [self.tra_dilution_3,self.use_tra_dilution_3], [self.tra_dilution_4,self.use_tra_dilution_4],
             [self.tra_dilution_5,self.use_tra_dilution_5], [self.tra_dilution_6,self.use_tra_dilution_6], 
             [self.tra_dilution_7,self.use_tra_dilution_7], [self.tra_dilution_8,self.use_tra_dilution_8], 
             [self.tra_dilution_9,self.use_tra_dilution_9], [self.tra_dilution_10,self.use_tra_dilution_10],
             [self.tra_dilution_11,self.use_tra_dilution_11], [self.tra_dilution_12,self.use_tra_dilution_12],
             [self.tra_dilution_13,self.use_tra_dilution_13], [self.tra_dilution_14,self.use_tra_dilution_14],
             [self.tra_dilution_15,self.use_tra_dilution_15], [self.tra_dilution_16,self.use_tra_dilution_16], 
             [self.tra_dilution_17,self.use_tra_dilution_17], [self.tra_dilution_18,self.use_tra_dilution_18], 
             [self.tra_dilution_19,self.use_tra_dilution_19], [self.tra_dilution_20,self.use_tra_dilution_20]


            ]
 
    return tra_dilution


def rvs_opt(self): 
    rvs_opt = [
             self.rvs_opt_1,self.rvs_opt_2,self.rvs_opt_3,self.rvs_opt_4,self.rvs_opt_5,
             self.rvs_opt_6,self.rvs_opt_7,self.rvs_opt_8,self.rvs_opt_9,self.rvs_opt_10,
             self.rvs_opt_11,self.rvs_opt_12,self.rvs_opt_13,self.rvs_opt_14,self.rvs_opt_15,
             self.rvs_opt_16,self.rvs_opt_17,self.rvs_opt_18,self.rvs_opt_19,self.rvs_opt_20
            ]
 
    return rvs_opt



def act_opt(self): 
    act_opt = [
             self.act_opt_1,self.act_opt_2,self.act_opt_3,self.act_opt_4,self.act_opt_5,
             self.act_opt_6,self.act_opt_7,self.act_opt_8,self.act_opt_9,self.act_opt_10
            ]
 
    return act_opt


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



def ast_data_to_planet(self): 
    
    ast_data_to_planet = [
            self.ast_data_planet_1,self.ast_data_planet_2,self.ast_data_planet_3,self.ast_data_planet_4,self.ast_data_planet_5,
            self.ast_data_planet_6,self.ast_data_planet_7,self.ast_data_planet_8,self.ast_data_planet_9,self.ast_data_planet_10,
            ]
   
    return ast_data_to_planet


def use_ast_data_to_planet(self): 
    
    use_ast_data_to_planet = [
            self.use_ast_data_1,self.use_ast_data_2,self.use_ast_data_3,self.use_ast_data_4,self.use_ast_data_5,
            self.use_ast_data_6,self.use_ast_data_7,self.use_ast_data_8,self.use_ast_data_9,self.use_ast_data_10,
            ]
   
    return use_ast_data_to_planet
    
    
def ast_data_to_planet_2(self): 
    
    ast_data_to_planet_2 = [
            self.ast_data_planet_hipp_1,self.ast_data_planet_hipp_2
            ]
   
    return ast_data_to_planet_2


def use_ast_data_to_planet_2(self): 
    
    use_ast_data_to_planet_2 = [
            self.use_ast_data_hipp_1,self.use_ast_data_hipp_2
            ]
   
    return use_ast_data_to_planet_2
    
        
    
    
