#!/usr/bin/python



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
        