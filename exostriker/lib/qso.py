#!/usr/bin/python
import gui_groups


def _keep_tabs(tabs, pages):
    """Select QSO pages by identity without shifting shared GUI tab indices."""
    allowed = {page for page, _ in pages}
    for index in range(tabs.count()):
        tabs.setTabVisible(index, tabs.widget(index) in allowed)
    for page, label in pages:
        tabs.setTabText(tabs.indexOf(page), label)
    if tabs.currentWidget() not in allowed:
        tabs.setCurrentWidget(pages[0][0])


def _configure_tabs(self):
    # Keep the original widgets and their indices: other GUI handlers navigate
    # these same pages, and new Exo-Striker tabs must not be relabeled as QSO GP.
    _keep_tabs(self.param_tabs, [
        (self.tab_RV_pl_param, 'Sine param.'), (self.tab_GP_param, 'GP param.'),
        (self.Limits_and_priors, 'Limits and Priors'),
        (self.More_param_tabs, 'Model opt.'), (self.plot_options_tabs, 'Plot opt.')])
    _keep_tabs(self.plot_tabs, [
        (self.tab_timeseries_RV, 'Data'), (self.tab_timeseries_act, 'Test data'),
        (self.Parameter_correlations, 'Sample corr.')])
    _keep_tabs(self.tab_ts_RV, [
        (self.RV_time_series_tab, 'Data and model'),
        (self.RV_timeseries_periodogram, 'GLS'), (self.tab_172, 'GLS o-c'),
        (self.tab_157, 'MLP'), (self.Window_func_periodograms, 'Window (DFT)')])
    _keep_tabs(self.tabWidget_data_section, [
        (self.tab_143, 'Data'), (self.tab_24, 'Test data'),
        (self.tab_52, 'Limits and Priors')])
    _keep_tabs(self.tabWidget_14, [(self.tab_60, 'GP')])
    _keep_tabs(self.tabWidget_34, [
        (self.tab_123, 'Model'), (self.tab_145, 'SciPyOp opt.')])
    _keep_tabs(self.model_param_tab, [
        (self.tab_122, 'Models'), (self.tab_44, 'GLS/MLP opt.'),
        (self.mcmc_param_tab, 'MCMC opt.'), (self.tab_121, 'Nest. Samp. opt.'),
        (self.tab_33, 'Auto fit opt.')])
    _keep_tabs(self.tabWidget, [
        (self.tab_6, 'Sine signal 1'), (self.tab_15, 'Sine signal 2'),
        (self.tab, 'Sine signal 3'), (self.tab_20, 'GP')])
    _keep_tabs(self.tabWidget_8, [(self.tab_2, 'Sine signals 1-3')])
    _keep_tabs(self.plot_opt_tab, [
        (self.tab_plot_axes, 'Axes'), (self.tab_11, 'Data / GP'),
        (self.tab_21, 'Test data'), (self.tab_25, 'GLS/MLP'),
        (self.tab_3, 'Data Insp.'), (self.cornerplot_plot_tab, 'Cornerplot')])
    _keep_tabs(self.tabWidget_10, [
        (self.tab_53, 'Offsets and Jitters'), (self.tab_80, 'Trends')])
    _keep_tabs(self.tabWidget_4, [
        (self.tab_144, 'Offsets and jitters'), (self.tab_156, 'Trends'),
        (self.tab_189, 'Data options')])
    _keep_tabs(self.tabWidget_7, [
        (self.tab_185, 'Phase plots'), (self.tab_186, 'Mirrored plots')])


def qso_mode(self):
    self.setWindowTitle("The QSO-Striker: Modeling the variability of quasar light-curves.")

    _configure_tabs(self)

    self.Planet3_11.setText("Offset min.")
    self.Planet3_13.setText("Offset max.")
    self.label_36.setText("Jitter min.")
    self.label_37.setText("Jitter max.")


    self.data_insp_load_data.setText("Load to Data")
    for name, title in (('RvsWindow', 'Data options'), ('ActivityWindow', 'Test data options')):
        window = getattr(self, name, None)
        if window is not None:
            window.setWindowTitle(title)


    self.use_mix_fitting.setVisible(False)
    self.label.setVisible(False)
    self.label_2.setVisible(False)
    self.label_134.setVisible(False)
    self.label_133.setVisible(False)
    self.mix_pl_1.setVisible(False)
    self.mix_pl_2.setVisible(False)
    self.mix_pl_3.setVisible(False)
    self.mix_pl_4.setVisible(False)
    self.mix_pl_5.setVisible(False)
    self.mix_pl_6.setVisible(False)
    self.mix_pl_7.setVisible(False)
    self.mix_pl_8.setVisible(False)
    self.mix_pl_9.setVisible(False)

    self.label_dyn_model_to_kill.setVisible(False)
    self.dyn_model_to_kill.setVisible(False)
    self.label_dyn_model_accuracy.setVisible(False)
    self.dyn_model_accuracy.setVisible(False)
    self.label_time_step_of_model.setVisible(False)
    self.time_step_model.setVisible(False)
    self.use_optim_stab_constraints.setVisible(False)
    self.optim_AMD.setVisible(False)
    self.optim_Nbody.setVisible(False)
    self.radioButton_omega_dot_free.setVisible(False)
    self.radioButton_omega_dot_GR.setVisible(False)
    self.deattach_omega_dot.setVisible(False)
    self.force_copl_incl.setVisible(False)


    self.actionopen_RVmod_init_file.setVisible(False)
    self.actionOpen_RVbank_file.setVisible(False)
    self.actionOpen_MAROONX_file.setVisible(False)
    self.actionGet_Transit_model.setVisible(False)
    self.actionGet_Transit_data.setVisible(False)
    self.actionGet_Orb_evol.setVisible(False)

    self.actionRV_Auto_fit.setText("Auto fit")

    self.do_RV_GP.setText("Use GP") 

    self.radioButton_ewm.setHidden(True)

    self.radioButton_RV.setHidden(True)
    self.radioButton_transit.setHidden(True)
    self.radioButton_transit_RV.setHidden(True)
    self.radioButton_ttv.setHidden(True)
    self.radioButton_ttv_RV.setHidden(True)
    self.radioButton_hkl.setHidden(True)
    self.radioButton_ast.setHidden(True)
    self.radioButton_ast_RV.setHidden(True)

    self.radioButton_Keplerian.setHidden(True)
    self.radioButton_Dynamical.setHidden(True)



    ############  Bounds ###################

    self.label_ecc_minmax_1.setVisible(False)
    self.use_e_bound_1.setVisible(False)
    self.e_min_1.setVisible(False)
    self.e_max_1.setVisible(False)

    self.label_omega_minmax_1.setVisible(False)
    self.use_om_bound_1.setVisible(False)
    self.om_min_1.setVisible(False)
    self.om_max_1.setVisible(False)

    self.label_incl_minmax_1.setVisible(False)
    self.use_inc_bound_1.setVisible(False)
    self.incl_min_1.setVisible(False)
    self.incl_max_1.setVisible(False)

    self.label_Omega_minmax_1.setVisible(False)
    self.use_Omega_bound_1.setVisible(False)
    self.Omega_min_1.setVisible(False)
    self.Omega_max_1.setVisible(False)

    self.label_t0_minmax_1.setVisible(False)
    self.use_t0_bound_1.setVisible(False)
    self.t0_min_1.setVisible(False)
    self.t0_max_1.setVisible(False) 

    self.label_pl_rad_minmax_1.setVisible(False)
    self.use_pl_rad_bound_1.setVisible(False)
    self.pl_rad_min_1.setVisible(False)
    self.pl_rad_max_1.setVisible(False) 

    self.label_a_sol_minmax_1.setVisible(False)
    self.use_a_sol_bound_1.setVisible(False)
    self.a_sol_min_1.setVisible(False)
    self.a_sol_max_1.setVisible(False) 

    self.label_K_minmax_28.setVisible(False)
    self.use_omega_dot_bound_1.setVisible(False)
    self.omega_dot_min_1.setVisible(False)
    self.omega_dot_max_1.setVisible(False) 


    self.label_K_minmax_29.setVisible(False)
    self.use_h_bound_1.setVisible(False)
    self.h_min_1.setVisible(False)
    self.h_max_1.setVisible(False) 

    self.label_K_minmax_30.setVisible(False)
    self.use_k_bound_1.setVisible(False)
    self.k_min_1.setVisible(False)
    self.k_max_1.setVisible(False) 

    self.label_K_minmax_31.setVisible(False) 
    self.use_lambda_bound_1.setVisible(False)
    self.lambda_min_1.setVisible(False)
    self.lambda_max_1.setVisible(False) 

    self.Planet1_32.setVisible(False) 
    self.Planet2_30.setVisible(False) 

    self.label_ecc_minmax_2.setVisible(False)
    self.use_e_bound_2.setVisible(False)
    self.e_min_2.setVisible(False)
    self.e_max_2.setVisible(False)

    self.label_omega_minmax_2.setVisible(False)
    self.use_om_bound_2.setVisible(False)
    self.om_min_2.setVisible(False)
    self.om_max_2.setVisible(False)

    self.label_incl_minmax_2.setVisible(False)
    self.use_inc_bound_2.setVisible(False)
    self.incl_min_2.setVisible(False)
    self.incl_max_2.setVisible(False)

    self.label_Omega_minmax_2.setVisible(False)
    self.use_Omega_bound_2.setVisible(False)
    self.Omega_min_2.setVisible(False)
    self.Omega_max_2.setVisible(False)

    self.label_t0_minmax_2.setVisible(False)
    self.use_t0_bound_2.setVisible(False)
    self.t0_min_2.setVisible(False)
    self.t0_max_2.setVisible(False) 

    self.label_pl_rad_minmax_2.setVisible(False)
    self.use_pl_rad_bound_2.setVisible(False)
    self.pl_rad_min_2.setVisible(False)
    self.pl_rad_max_2.setVisible(False) 

    self.label_a_sol_minmax_2.setVisible(False)
    self.use_a_sol_bound_2.setVisible(False)
    self.a_sol_min_2.setVisible(False)
    self.a_sol_max_2.setVisible(False) 

    self.label_K_minmax_32.setVisible(False)
    self.use_omega_dot_bound_2.setVisible(False)
    self.omega_dot_min_2.setVisible(False)
    self.omega_dot_max_2.setVisible(False) 


    self.label_K_minmax_33.setVisible(False)
    self.use_h_bound_2.setVisible(False)
    self.h_min_2.setVisible(False)
    self.h_max_2.setVisible(False) 

    self.label_K_minmax_34.setVisible(False)
    self.use_k_bound_2.setVisible(False)
    self.k_min_2.setVisible(False)
    self.k_max_2.setVisible(False) 

    self.label_K_minmax_35.setVisible(False) 
    self.use_lambda_bound_2.setVisible(False)
    self.lambda_min_2.setVisible(False)
    self.lambda_max_2.setVisible(False) 

    self.Planet1_33.setVisible(False) 
    self.Planet2_31.setVisible(False) 


    self.label_ecc_minmax_3.setVisible(False)
    self.use_e_bound_3.setVisible(False)
    self.e_min_3.setVisible(False)
    self.e_max_3.setVisible(False)

    self.label_omega_minmax_3.setVisible(False)
    self.use_om_bound_3.setVisible(False)
    self.om_min_3.setVisible(False)
    self.om_max_3.setVisible(False)

    self.label_incl_minmax_3.setVisible(False)
    self.use_inc_bound_3.setVisible(False)
    self.incl_min_3.setVisible(False)
    self.incl_max_3.setVisible(False)

    self.label_Omega_minmax_3.setVisible(False)
    self.use_Omega_bound_3.setVisible(False)
    self.Omega_min_3.setVisible(False)
    self.Omega_max_3.setVisible(False)

    self.label_t0_minmax_3.setVisible(False)
    self.use_t0_bound_3.setVisible(False)
    self.t0_min_3.setVisible(False)
    self.t0_max_3.setVisible(False) 

    self.label_pl_rad_minmax_3.setVisible(False)
    self.use_pl_rad_bound_3.setVisible(False)
    self.pl_rad_min_3.setVisible(False)
    self.pl_rad_max_3.setVisible(False) 

    self.label_a_sol_minmax_3.setVisible(False)
    self.use_a_sol_bound_3.setVisible(False)
    self.a_sol_min_3.setVisible(False)
    self.a_sol_max_3.setVisible(False) 

    self.label_K_minmax_36.setVisible(False)
    self.use_omega_dot_bound_3.setVisible(False)
    self.omega_dot_min_3.setVisible(False)
    self.omega_dot_max_3.setVisible(False) 


    self.label_K_minmax_37.setVisible(False)
    self.use_h_bound_3.setVisible(False)
    self.h_min_3.setVisible(False)
    self.h_max_3.setVisible(False) 

    self.label_K_minmax_38.setVisible(False)
    self.use_k_bound_3.setVisible(False)
    self.k_min_3.setVisible(False)
    self.k_max_3.setVisible(False) 

    self.label_K_minmax_39.setVisible(False) 
    self.use_lambda_bound_3.setVisible(False)
    self.lambda_min_3.setVisible(False)
    self.lambda_max_3.setVisible(False) 

    self.Planet1_34.setVisible(False) 
    self.Planet2_32.setVisible(False) 

    ############  Bounds End ###################


   ############  Gauss ###################

    self.label_ecc_minmax_10.setVisible(False)
    self.use_e_norm_pr_1.setVisible(False)
    self.e_mean_1.setVisible(False)
    self.e_sigma_1.setVisible(False)

    self.label_omega_minmax_10.setVisible(False)
    self.use_om_norm_pr_1.setVisible(False)
    self.om_mean_1.setVisible(False)
    self.om_sigma_1.setVisible(False)

    self.label_incl_minmax_10.setVisible(False)
    self.use_incl_norm_pr_1.setVisible(False)
    self.incl_mean_1.setVisible(False)
    self.incl_sigma_1.setVisible(False)

    self.label_Omega_minmax_10.setVisible(False)
    self.use_Omega_norm_pr_1.setVisible(False)
    self.Omega_mean_1.setVisible(False)
    self.Omega_sigma_1.setVisible(False)

    self.label_t0_minmax_10.setVisible(False)
    self.use_t0_norm_pr_1.setVisible(False)
    self.t0_mean_1.setVisible(False)
    self.t0_sigma_1.setVisible(False) 

    self.label_pl_rad_minmax_10.setVisible(False)
    self.use_pl_rad_norm_pr_1.setVisible(False)
    self.pl_rad_mean_1.setVisible(False)
    self.pl_rad_sigma_1.setVisible(False) 

    self.label_a_sol_minmax_10.setVisible(False)
    self.use_a_sol_norm_pr_1.setVisible(False)
    self.a_sol_mean_1.setVisible(False)
    self.a_sol_sigma_1.setVisible(False) 

    self.label_K_minmax_66.setVisible(False)
    self.use_omega_dot_norm_pr_1.setVisible(False)
    self.omega_dot_mean_1.setVisible(False)
    self.omega_dot_sigma_1.setVisible(False) 


    self.label_K_minmax_65.setVisible(False)
    self.use_h_norm_pr_1.setVisible(False)
    self.h_mean_1.setVisible(False)
    self.h_sigma_1.setVisible(False) 

    self.label_K_minmax_64.setVisible(False)
    self.use_k_norm_pr_1.setVisible(False)
    self.k_mean_1.setVisible(False)
    self.k_sigma_1.setVisible(False) 


    self.label_K_minmax_67.setVisible(False) 
    self.use_lambda_norm_pr_1.setVisible(False)
    self.lambda_mean_1.setVisible(False)
    self.lambda_sigma_1.setVisible(False) 



    self.label_ecc_minmax_11.setVisible(False)
    self.use_e_norm_pr_2.setVisible(False)
    self.e_mean_2.setVisible(False)
    self.e_sigma_2.setVisible(False)

    self.label_omega_minmax_11.setVisible(False)
    self.use_om_norm_pr_2.setVisible(False)
    self.om_mean_2.setVisible(False)
    self.om_sigma_2.setVisible(False)

    self.label_incl_minmax_11.setVisible(False)
    self.use_incl_norm_pr_2.setVisible(False)
    self.incl_mean_2.setVisible(False)
    self.incl_sigma_2.setVisible(False)

    self.label_Omega_minmax_11.setVisible(False)
    self.use_Omega_norm_pr_2.setVisible(False)
    self.Omega_mean_2.setVisible(False)
    self.Omega_sigma_2.setVisible(False)

    self.label_t0_minmax_11.setVisible(False)
    self.use_t0_norm_pr_2.setVisible(False)
    self.t0_mean_2.setVisible(False)
    self.t0_sigma_2.setVisible(False) 

    self.label_pl_rad_minmax_11.setVisible(False)
    self.use_pl_rad_norm_pr_2.setVisible(False)
    self.pl_rad_mean_2.setVisible(False)
    self.pl_rad_sigma_2.setVisible(False) 

    self.label_a_sol_minmax_11.setVisible(False)
    self.use_a_sol_norm_pr_2.setVisible(False)
    self.a_sol_mean_2.setVisible(False)
    self.a_sol_sigma_2.setVisible(False) 

    self.label_K_minmax_70.setVisible(False)
    self.use_omega_dot_norm_pr_2.setVisible(False)
    self.omega_dot_mean_2.setVisible(False)
    self.omega_dot_sigma_2.setVisible(False) 


    self.label_K_minmax_69.setVisible(False)
    self.use_h_norm_pr_2.setVisible(False)
    self.h_mean_2.setVisible(False)
    self.h_sigma_2.setVisible(False) 

    self.label_K_minmax_68.setVisible(False)
    self.use_k_norm_pr_2.setVisible(False)
    self.k_mean_2.setVisible(False)
    self.k_sigma_2.setVisible(False) 


    self.label_K_minmax_71.setVisible(False) 
    self.use_lambda_norm_pr_2.setVisible(False)
    self.lambda_mean_2.setVisible(False)
    self.lambda_sigma_2.setVisible(False) 


    self.label_ecc_minmax_12.setVisible(False)
    self.use_e_norm_pr_3.setVisible(False)
    self.e_mean_3.setVisible(False)
    self.e_sigma_3.setVisible(False)

    self.label_omega_minmax_12.setVisible(False)
    self.use_om_norm_pr_3.setVisible(False)
    self.om_mean_3.setVisible(False)
    self.om_sigma_3.setVisible(False)

    self.label_incl_minmax_12.setVisible(False)
    self.use_incl_norm_pr_3.setVisible(False)
    self.incl_mean_3.setVisible(False)
    self.incl_sigma_3.setVisible(False)

    self.label_Omega_minmax_12.setVisible(False)
    self.use_Omega_norm_pr_3.setVisible(False)
    self.Omega_mean_3.setVisible(False)
    self.Omega_sigma_3.setVisible(False)

    self.label_t0_minmax_12.setVisible(False)
    self.use_t0_norm_pr_3.setVisible(False)
    self.t0_mean_3.setVisible(False)
    self.t0_sigma_3.setVisible(False) 

    self.label_pl_rad_minmax_12.setVisible(False)
    self.use_pl_rad_norm_pr_3.setVisible(False)
    self.pl_rad_mean_3.setVisible(False)
    self.pl_rad_sigma_3.setVisible(False) 

    self.label_a_sol_minmax_12.setVisible(False)
    self.use_a_sol_norm_pr_3.setVisible(False)
    self.a_sol_mean_3.setVisible(False)
    self.a_sol_sigma_3.setVisible(False) 

    self.label_K_minmax_74.setVisible(False)
    self.use_omega_dot_norm_pr_3.setVisible(False)
    self.omega_dot_mean_3.setVisible(False)
    self.omega_dot_sigma_3.setVisible(False) 


    self.label_K_minmax_73.setVisible(False)
    self.use_h_norm_pr_3.setVisible(False)
    self.h_mean_3.setVisible(False)
    self.h_sigma_3.setVisible(False) 

    self.label_K_minmax_72.setVisible(False)
    self.use_k_norm_pr_3.setVisible(False)
    self.k_mean_3.setVisible(False)
    self.k_sigma_3.setVisible(False) 


    self.label_K_minmax_75.setVisible(False) 
    self.use_lambda_norm_pr_3.setVisible(False)
    self.lambda_mean_3.setVisible(False)
    self.lambda_sigma_3.setVisible(False) 

    
    self.Planet1_41.setVisible(False) 
    self.Planet2_39.setVisible(False) 
    self.Planet1_42.setVisible(False) 
    self.Planet2_40.setVisible(False) 

    self.Planet1_43.setVisible(False) 
    self.Planet2_41.setVisible(False) 

    ############  Gaus End ###################

    ############  Jeff ###################

    self.label_ecc_minmax_19.setVisible(False)
    self.use_e_jeff_pr_1.setVisible(False)
    self.e_jeff_alpha_1.setVisible(False)
    self.e_jeff_beta_1.setVisible(False)

    self.label_omega_minmax_19.setVisible(False)
    self.use_om_jeff_pr_1.setVisible(False)
    self.om_jeff_alpha_1.setVisible(False)
    self.om_jeff_beta_1.setVisible(False)

    self.label_incl_minmax_19.setVisible(False)
    self.use_incl_jeff_pr_1.setVisible(False)
    self.incl_jeff_alpha_1.setVisible(False)
    self.incl_jeff_beta_1.setVisible(False)

    self.label_Omega_minmax_19.setVisible(False)
    self.use_Omega_jeff_pr_1.setVisible(False)
    self.Omega_jeff_alpha_1.setVisible(False)
    self.Omega_jeff_beta_1.setVisible(False)

    self.label_t0_minmax_19.setVisible(False)
    self.use_t0_jeff_pr_1.setVisible(False)
    self.t0_jeff_alpha_1.setVisible(False)
    self.t0_jeff_beta_1.setVisible(False) 

    self.label_pl_rad_minmax_19.setVisible(False)
    self.use_pl_rad_jeff_pr_1.setVisible(False)
    self.pl_rad_jeff_alpha_1.setVisible(False)
    self.pl_rad_jeff_beta_1.setVisible(False) 

    self.label_a_sol_minmax_19.setVisible(False)
    self.use_a_sol_jeff_pr_1.setVisible(False)
    self.a_sol_jeff_alpha_1.setVisible(False)
    self.a_sol_jeff_beta_1.setVisible(False) 

    self.label_K_minmax_102.setVisible(False)
    self.use_omega_dot_jeff_pr_1.setVisible(False)
    self.omega_dot_alpha_1.setVisible(False)
    self.omega_dot_beta_1.setVisible(False) 


    self.label_K_minmax_101.setVisible(False)
    self.use_h_jeff_pr_1.setVisible(False)
    self.h_alpha_1.setVisible(False)
    self.h_beta_1.setVisible(False) 

    self.label_K_minmax_100.setVisible(False)
    self.use_k_jeff_pr_1.setVisible(False)
    self.k_alpha_1.setVisible(False)
    self.k_beta_1.setVisible(False)  


    self.label_K_minmax_103.setVisible(False) 
    self.use_lambda_jeff_pr_1.setVisible(False)
    self.lambda_alpha_1.setVisible(False)
    self.lambda_beta_1.setVisible(False) 



    self.label_ecc_minmax_20.setVisible(False)
    self.use_e_jeff_pr_2.setVisible(False)
    self.e_jeff_alpha_2.setVisible(False)
    self.e_jeff_beta_2.setVisible(False)

    self.label_omega_minmax_20.setVisible(False)
    self.use_om_jeff_pr_2.setVisible(False)
    self.om_jeff_alpha_2.setVisible(False)
    self.om_jeff_beta_2.setVisible(False)

    self.label_incl_minmax_20.setVisible(False)
    self.use_incl_jeff_pr_2.setVisible(False)
    self.incl_jeff_alpha_2.setVisible(False)
    self.incl_jeff_beta_2.setVisible(False)

    self.label_Omega_minmax_20.setVisible(False)
    self.use_Omega_jeff_pr_2.setVisible(False)
    self.Omega_jeff_alpha_2.setVisible(False)
    self.Omega_jeff_beta_2.setVisible(False)

    self.label_t0_minmax_20.setVisible(False)
    self.use_t0_jeff_pr_2.setVisible(False)
    self.t0_jeff_alpha_2.setVisible(False)
    self.t0_jeff_beta_2.setVisible(False) 

    self.label_pl_rad_minmax_20.setVisible(False)
    self.use_pl_rad_jeff_pr_2.setVisible(False)
    self.pl_rad_jeff_alpha_2.setVisible(False)
    self.pl_rad_jeff_beta_2.setVisible(False) 

    self.label_a_sol_minmax_20.setVisible(False)
    self.use_a_sol_jeff_pr_2.setVisible(False)
    self.a_sol_jeff_alpha_2.setVisible(False)
    self.a_sol_jeff_beta_2.setVisible(False) 

    self.label_K_minmax_106.setVisible(False)
    self.use_omega_dot_jeff_pr_2.setVisible(False)
    self.omega_dot_alpha_2.setVisible(False)
    self.omega_dot_beta_2.setVisible(False) 


    self.label_K_minmax_105.setVisible(False)
    self.use_h_jeff_pr_2.setVisible(False)
    self.h_alpha_2.setVisible(False)
    self.h_beta_2.setVisible(False) 

    self.label_K_minmax_104.setVisible(False)
    self.use_k_jeff_pr_2.setVisible(False)
    self.k_alpha_2.setVisible(False)
    self.k_beta_2.setVisible(False)  


    self.label_K_minmax_107.setVisible(False) 
    self.use_lambda_jeff_pr_2.setVisible(False)
    self.lambda_alpha_2.setVisible(False)
    self.lambda_beta_2.setVisible(False) 


    self.label_ecc_minmax_21.setVisible(False)
    self.use_e_jeff_pr_3.setVisible(False)
    self.e_jeff_alpha_3.setVisible(False)
    self.e_jeff_beta_3.setVisible(False)

    self.label_omega_minmax_21.setVisible(False)
    self.use_om_jeff_pr_3.setVisible(False)
    self.om_jeff_alpha_3.setVisible(False)
    self.om_jeff_beta_3.setVisible(False)

    self.label_incl_minmax_21.setVisible(False)
    self.use_incl_jeff_pr_3.setVisible(False)
    self.incl_jeff_alpha_3.setVisible(False)
    self.incl_jeff_beta_3.setVisible(False)

    self.label_Omega_minmax_21.setVisible(False)
    self.use_Omega_jeff_pr_3.setVisible(False)
    self.Omega_jeff_alpha_3.setVisible(False)
    self.Omega_jeff_beta_3.setVisible(False)

    self.label_t0_minmax_21.setVisible(False)
    self.use_t0_jeff_pr_3.setVisible(False)
    self.t0_jeff_alpha_3.setVisible(False)
    self.t0_jeff_beta_3.setVisible(False) 

    self.label_pl_rad_minmax_21.setVisible(False)
    self.use_pl_rad_jeff_pr_3.setVisible(False)
    self.pl_rad_jeff_alpha_3.setVisible(False)
    self.pl_rad_jeff_beta_3.setVisible(False) 

    self.label_a_sol_minmax_21.setVisible(False)
    self.use_a_sol_jeff_pr_3.setVisible(False)
    self.a_sol_jeff_alpha_3.setVisible(False)
    self.a_sol_jeff_beta_3.setVisible(False) 

    self.label_K_minmax_110.setVisible(False)
    self.use_omega_dot_jeff_pr_3.setVisible(False)
    self.omega_dot_alpha_3.setVisible(False)
    self.omega_dot_beta_3.setVisible(False) 


    self.label_K_minmax_109.setVisible(False)
    self.use_h_jeff_pr_3.setVisible(False)
    self.h_alpha_3.setVisible(False)
    self.h_beta_3.setVisible(False) 

    self.label_K_minmax_108.setVisible(False)
    self.use_k_jeff_pr_3.setVisible(False)
    self.k_alpha_3.setVisible(False)
    self.k_beta_3.setVisible(False)  


    self.label_K_minmax_111.setVisible(False) 
    self.use_lambda_jeff_pr_3.setVisible(False)
    self.lambda_alpha_3.setVisible(False)
    self.lambda_beta_3.setVisible(False) 


    self.Planet2_48.setVisible(False) 
    self.Planet1_50.setVisible(False) 
    self.Planet1_51.setVisible(False) 
    self.Planet1_52.setVisible(False) 
    self.Planet2_50.setVisible(False) 
    self.Planet2_49.setVisible(False) 

    ############  Jeff End ###################

    self.label_Ma_minmax_1.setText("Phase [deg]") 
    self.label_Ma_minmax_2.setText("Phase [deg]") 
    self.label_Ma_minmax_3.setText("Phase [deg]")  

    self.label_Ma_minmax_10.setText("Phase [deg]") 
    self.label_Ma_minmax_11.setText("Phase [deg]") 
    self.label_Ma_minmax_12.setText("Phase [deg]")  

    self.label_Ma_minmax_19.setText("Phase [deg]") 
    self.label_Ma_minmax_20.setText("Phase [deg]") 
    self.label_Ma_minmax_21.setText("Phase [deg]")  





    self.label_109.setText("Test data symbol size global") 

    self.RV_plot_cross_hair.setText("Plot crosshair") 
    self.RV_o_c_plot_cross_hair.setText("o-c Plot crosshair") 
    self.plot_RV_GP_model.setText("subtract GP model from o-c") 

    self.label_18.setText("Sine minimizer algorithm")
    self.checkBox_first_RV_epoch.setText("First data [BJD]")
    self.label_kep_model_to_kill.setText("Fit interrupt [s]")





    ####### planet tabs ############

    self.Planet1.setText("Signal 1") 
    self.Planet2.setText("Signal 2") 
    self.Planet3.setText("Signal 3") 

    self.e1.setVisible(False) 
    self.e2.setVisible(False) 
    self.e3.setVisible(False) 

    self.use_e1.setVisible(False) 
    self.use_e2.setVisible(False) 
    self.use_e3.setVisible(False) 

    self.err_e1.setVisible(False) 
    self.err_e2.setVisible(False) 
    self.err_e3.setVisible(False)


    self.om1.setVisible(False) 
    self.om2.setVisible(False) 
    self.om3.setVisible(False) 

    self.err_om1.setVisible(False) 
    self.err_om2.setVisible(False) 
    self.err_om3.setVisible(False)

    self.use_om1.setVisible(False) 
    self.use_om2.setVisible(False) 
    self.use_om3.setVisible(False)


    self.incl1.setVisible(False) 
    self.incl2.setVisible(False) 
    self.incl3.setVisible(False) 

    self.err_incl1.setVisible(False) 
    self.err_incl2.setVisible(False) 
    self.err_incl3.setVisible(False) 

    self.use_incl1.setVisible(False) 
    self.use_incl2.setVisible(False) 
    self.use_incl3.setVisible(False) 

    self.Omega1.setVisible(False) 
    self.Omega2.setVisible(False) 
    self.Omega3.setVisible(False) 

    self.err_Omega1.setVisible(False) 
    self.err_Omega2.setVisible(False) 
    self.err_Omega3.setVisible(False) 

    self.use_Omega1.setVisible(False) 
    self.use_Omega2.setVisible(False) 
    self.use_Omega3.setVisible(False) 


    self.label_omega.setVisible(False)
    self.label_ecc.setVisible(False)
    self.label_incl.setVisible(False)

    self.label_Omega.setVisible(False)
    self.label_omega_2.setVisible(False)
    self.label_55.setVisible(False)


    self.label_a.setVisible(False)
    self.label_a1.setVisible(False)
    self.label_a2.setVisible(False)
    self.label_a3.setVisible(False)

    self.label_mp.setVisible(False)
    self.label_mass1.setVisible(False)
    self.label_mass2.setVisible(False)
    self.label_mass3.setVisible(False)

    self.label_tw.setVisible(False)
    self.label_t_peri1.setVisible(False)
    self.label_t_peri2.setVisible(False)
    self.label_t_peri3.setVisible(False)

    self.label_Rp.setVisible(False)
    self.label_ap.setVisible(False)





    self.label_19.setText("[mag./day]") 
    self.label_353.setText("[mag./day]") 
   # self.label_393.setText("<html><head/><body><p>Quad. add to &sigma; [mag.]</p></body></html>") 

    self.label_auto_fit_N_planets.setText("Maximum number of signals to look for") 


    self.actionGet_RV_model.setText("Get model")
    self.actionGet_RV_data.setText("Get data")
    self.label_322.setText("Data GLS and Test data GLS")

    self.label_382.setText("MLP")
    self.label_271.setText("GLS")
    self.label_312.setText("GLS o-c")
    self.label_381.setText("MLP")
    self.label_379.setText("MLP o-c")


    self.label_323.setVisible(False)
    self.label_316.setVisible(False)
    self.N_TLS_peak_to_point.setVisible(False)
    self.label_Rp.setVisible(False)
    self.label_Rp.setVisible(False)

    self.Planet3_4.setText("Offset [mag.]") 
    self.Planet3_6.setText("Jitter [mag.]") 

    self.label_stellar_mass.setVisible(False)
    self.value_stellar_mass.setVisible(False)

    self.label_K.setText("Ampl. [mag.]") 
    self.label_Ma.setText("Phase [deg]") 

    self.label_K_minmax_1.setText("Ampl. [mag.]") 
    self.label_K_minmax_2.setText("Ampl. [mag.]") 
    self.label_K_minmax_3.setText("Ampl. [mag.]") 

    self.label_K_minmax_10.setText("Ampl. [mag.]") 
    self.label_K_minmax_11.setText("Ampl. [mag.]") 
    self.label_K_minmax_12.setText("Ampl. [mag.]") 
    self.label_K_minmax_19.setText("Ampl. [mag.]") 
    self.label_K_minmax_20.setText("Ampl. [mag.]") 
    self.label_K_minmax_21.setText("Ampl. [mag.]") 


    self.button_orb_evol.setVisible(False)


    self.button_auto_fit.setText("Auto fit") 

    self.param_gui_tr      = gui_groups.param_gui_tr(self)
    self.use_param_gui_tr  = gui_groups.use_param_gui_tr(self)
    self.param_gui_wd      = gui_groups.param_gui_wd(self)
    self.use_param_gui_wd  = gui_groups.use_param_gui_wd(self)
    self.param_errors_gui_wd = gui_groups.param_errors_gui_wd(self)
    self.err_t0          = gui_groups.err_t0(self)
    self.err_pl_rad      = gui_groups.err_pl_rad(self)
    self.err_a_sol       = gui_groups.err_a_sol(self)

    for i in range(3*9):

        self.param_gui_tr[i].setVisible(False) 
        self.use_param_gui_tr[i].setVisible(False) 



    for i in range(9):
        self.param_gui_wd[i].setVisible(False) 
        self.use_param_gui_wd[i].setVisible(False) 
        self.param_errors_gui_wd[i].setVisible(False) 

        self.err_t0[i].setVisible(False) 
        self.err_pl_rad[i].setVisible(False)  
        self.err_a_sol[i].setVisible(False)  


    self.dataInspector_HARPS_RVBank.setVisible(False)


    self.label_370.setVisible(False)
    self.AMD_led.setVisible(False)



    self.label_52.setVisible(False)
    self.label_50.setVisible(False)
    self.label_46.setVisible(False)
    self.tls_period_min_use.setVisible(False)
    self.tls_period_max_use.setVisible(False)
    self.label_446.setVisible(False)
    self.label_447.setVisible(False)
    self.label_448.setVisible(False)
    self.tls_grid_step.setVisible(False)
    self.tls_ofac.setVisible(False)
    self.tls_min_period.setVisible(False)
    self.tls_max_period.setVisible(False)
    self.tls_fap_1.setVisible(False)
    self.tls_fap_2.setVisible(False)
    self.tls_fap_3.setVisible(False)
    self.line_9.setVisible(False)



    return 
################################################################################################
