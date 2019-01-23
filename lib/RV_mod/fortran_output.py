 #!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
sys.path.insert(0, '../lib')
import numpy as np
from Warning_log import Warning_log
from parameters import parameters,  parameter_errors
from functions import *
from kernel import kernel, summary



# We also need a class for the fortran output. Originally fortran code outputted everything in several files, now we want it all in console output, this causes some problems which will be solved by applying functions of this class

class fortran_output(object):
    # it's better to give initial values for these attributes, just in case the user calls functions in the wrong order
    best_par=[]
    RV_kep=[]
    keplerian_fit=[]
    params=parameters([],[],[],0.0,0.0)
    param_errors=parameter_errors([],[],[],0.0,0.0)
    ndata_str=''
    mfit_str=''
    rms_str=''
    chi_str=''
    epoch_str=''
    masses=[]
    semiM=[]
    JD_model=[]
    model=[]    
    jd=[]
    o_c=[]
    rv_obs=[]
    rv_error=[]
    data_set=[]
    mfit = 0
    loglik=0 
    reduced_chi2=1
    chi2=1
    rms=1
    
    def __init__(self,text, npl, ndset, stellar_mass):
        self.text=text
        self.npl=npl
        self.ndset=ndset
        self.stellar_mass=stellar_mass
        
    def dismantle_RV_kep(self): # save each column of RV_kep in a separate attribute 
        a=np.transpose(np.array(self.RV_kep)) # so columns are rows now
        if (len(a)==6): # otherwise something is wrong or keplerian_fit is empty, so we do nothing
            # we need to remove lines in which one of these entries was wrong

            indices=np.arange(min(list(map(len,a)))) # we will be deleting some of them as we go along
            self.jd,wrong_indices=convert_array_to_float(a[0],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.rv_obs,wrong_indices=convert_array_to_float(a[2],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.o_c,wrong_indices=convert_array_to_float(a[3],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.rv_error,wrong_indices=convert_array_to_float(a[4],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.data_set,wrong_indices=convert_array_to_int(a[5],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            # and now we filter out wrong lines
            self.jd=self.jd[indices]
            self.rv_obs=self.rv_obs[indices]
            self.o_c=self.o_c[indices]
            self.rv_error=self.rv_error[indices]
            self.data_set=self.data_set[indices]
            self.data_set=np.array([int(i)-1 for i in self.data_set]) # cause fortran names datasets starting from 1, and python from 0
        
    def dismantle_keplerian_fit(self): # save each column of keplerian_fit in a separate attribute 
        a=np.transpose(np.array(self.keplerian_fit, dtype=float)) # so columns are rows now
        if (len(a)==2): # otherwise something is wrong or keplerian_fit is empty, so we do nothing
            self.JD_model=a[0]
            self.model=a[1]
        return    
        
    def sort_out_text (self): # initially self.text will contain all that would normally be in fortran output files, we need to separate it into parts
        # similar action to read_file_as_array_of_arrays, except on a string
        T=self.text
        l=len(T)
        fortran_output_warnings=Warning_log([],'Reading fortran output')
        # now a loop which will do the sorting out
        i=0
        while i<l:
            if(len(T[i])==0): # skip empty lines
                i=i+1
            elif(T[i][0]=='ITMAX'): # error in amoeba
                raise RuntimeError ('Runtime error in amoeba (fortran subroutine). Change initial parameters.')
            elif(T[i][0]=='loglik,'): # these two lines always appear in the output
                self.loglik=float(T[i+1][0]) 
                self.reduced_chi2=float(T[i+1][1])
                self.chi2=float(T[i+1][2])
                self.rms=float(T[i+1][3])
                i=i+2
            elif(T[i][0]=='Best-fit'): # characteristic opening for the best_par.out file
                i0=i # save the index where file begins
                while 1:
                    i=i+1
                    if (i>l-1):
                        raise Exception ('Wrong fortran output')
                    if(len(T[i])==0): # empty lines
                        raise Exception ('Empty line in best_par part of fortran output in a place where it is not allowed!')
                    if (T[i][0]=='Jupiter'): # fourth last line of best_par.out begins this way
                        i=i+4 # so this will be the first line that doesn't belong to this file
                        self.best_par=T[i0:i] # now we can save the part of T which corresponds to best_par
                        break # and leave the while 1 loop

            elif(len(T[i])==6): # other two options are RV.out and keplerian_fit.out, the first one has 6 columns
                i0=i # save the index where file begins
                while 1:
                    i=i+1
                    if (i==l): # that means we reached the end of the file, so let's save what we have and leave
                        self.RV_kep=T[i0:i] # now we can save the part of T which corresponds to RV
                        break # and leave the while 1 loop                      
                    if (len(T[i])!=6 or not is_float(T[i][0])): # that means we went outside of the part of the output corresponding to RV_kep.out, second condition just in case
                        self.RV_kep=T[i0:i] # now we can save the part of T which corresponds to RV
                        break # and leave the while 1 loop     
                    
            elif(len(T[i])==2): # last option is keplerian_fit.out, which has 2 columns
                i0=i # save the index where file begins
                while 1:
                    i=i+1
                    if (i==l): # that means we reached the end of the file, so let's save what we have and leave
                        self.keplerian_fit=T[i0:i] # now we can save the part of T which corresponds to fit
                        break # and leave the while 1 loop                      
                    if (len(T[i])!=2 or not is_float(T[i][0])): # that means we went outside of the part of the output corresponding to RV_kep.out, second condition just in case
                        self.keplerian_fit=T[i0:i] # now we can save the part of T which corresponds to fit
                        break # and leave the while 1 loop                         
            else: # just in case there is something wrong with fortran output, skip this line
                i=i+1
                fortran_output_warnings.update_warning_list('Wrong data format in line %d of fortran output, line skipped.'%i) #i, not i+1, becase we increase it above    
        fortran_output_warnings.print_warning_log()
        return
    
    # this function is a bit of a work around, it puts together the part of the output which corresponds to best_par.out
    def print_stat_array(self):
        string_to_be_printed=''
        for line in self.best_par:
            for element in line:
                string_to_be_printed+='%s '%element
            string_to_be_printed+='\n'
        print(string_to_be_printed)
        return    

    # this one on the other hand reads bes_par.out and saves fitted parameters in object attributes
    def save_stat_array(self):
        planet_params=[]
        planet_params_errors=[]
        offsets=[]
        offset_errors=[]
        jitters=[]
        jitter_errors=[]
        linear_trend=0.0
        linear_trend_error=0.0
        i=0
        l=len(self.best_par)
        if(l==0): # to know that if there was no best_par in the output we do not really save anything here
            self.stat_array_saved=False
        else:
            self.stat_array_saved=True
        fortran_stat_warnings=Warning_log([],'Saving best fit parameters')
        while i<l:
            if (self.best_par[i][1]=='K'):# characteristic of beginning of planet parameters and errors information 
                for k in range(self.npl):
                    planet_params  = np.concatenate((planet_params,np.array(list(map(float,self.best_par[i+1+(k*2)][:7])))))
                    planet_params_errors  = np.concatenate((planet_params_errors,np.array(list(map(float,self.best_par[i+2+(k*2)][:7])))))
                i=i+2*self.npl+1
            elif (self.best_par[i][1]=='V0'):# characteristic of beginning of RV offsets and their errors information 
                for k in range(self.ndset):
                    offsets  = np.concatenate((offsets,np.array(list(map(float,self.best_par[i+1+(k*2)])))))
                    offset_errors  = np.concatenate((offset_errors,np.array(list(map(float,self.best_par[i+2+(k*2)])))))
                i=i+2*self.ndset+1
            elif (self.best_par[i][0]=='Jitters'):# characteristic of beginning of RV offsets and their errors information 
                for k in range(self.ndset):
                    jitters  = np.concatenate((jitters,np.array(list(map(float,self.best_par[i+1+(k*2)])))))
                    jitter_errors  = np.concatenate((jitter_errors,np.array(list(map(float,self.best_par[i+2+(k*2)])))))
                i=i+2*self.ndset+1
            elif (self.best_par[i][0]=='linear'): # characteristic of beginning of linear trend and it's error information
                linear_trend  = float(self.best_par[i+1][0])
                linear_trend_error  = float(self.best_par[i+2][0])
                i=i+3
            elif (self.best_par[i][0]=='ndata'): # rest of the information is below that
                self.ndata_str = self.best_par[i]
                self.mfit_str = self.best_par[i+1]
                self.rms_str   = self.best_par[i+2]
                self.chi_str   = self.best_par[i+3]
                self.epoch_str = self.best_par[i+4]
                self.masses  = list(map(float,self.best_par[i+6]))
                self.semiM = list(map(float,self.best_par[i+8]))
                self.mfit = int(self.mfit_str[2])
                i=i+9
            else:
                i=i+1
                fortran_stat_warnings.update_warning_list('Wrong data format in line %d of best fit parameter information in fortran output, line skipped. Please check if number of planets and number of datasets is specified correctly!'%i) #i, not i+1, becase we increase it above    
        self.params=parameters(offsets,jitters,planet_params,linear_trend,self.stellar_mass)
        self.param_errors=parameter_errors(offset_errors,jitter_errors,planet_params_errors,linear_trend_error,0.0)
        fortran_stat_warnings.print_warning_log()
        return
     
    def generate_summary(self):
        if(len(self.ndata_str)>2 and len(self.mfit_str)>2):
            dof = float(self.ndata_str[2]) - float(self.mfit_str[2]) # degrees of freedom for summary 
        else:
            dof= 1 # we have to give something y'know
        stat = summary(self.params,self.param_errors,dof=dof)
        return stat
       
    def modfit(self, print_stat=True):
        self.sort_out_text()
        if(print_stat): # if requested by user (passing print_stat=True to modfit) prints out the output file of the fortran code describing best fit parameters (by default True)
            self.print_stat_array()
        self.save_stat_array()
        self.dismantle_keplerian_fit()
        self.dismantle_RV_kep()
        results = kernel(self.generate_summary(), self.jd, self.rv_obs, self.rv_error,self.o_c, self.model, self.JD_model, self.npl,self.semiM,self.masses,self.data_set,self.stat_array_saved,self.reduced_chi2,self.chi2,self.rms,self.loglik, self.mfit) 
        return results
