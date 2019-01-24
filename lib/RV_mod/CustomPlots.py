#!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
#sys.path.insert(0, '../lib')
import numpy as np
from Warning_log import Warning_log
import matplotlib.pyplot as plt
plt.switch_backend('SVG') 


from matplotlib import gridspec





class Customplots(object):


    def plot_gls(self, sig=np.array([0.1,0.01,0.001]),ind_sig=2, filetosave='periodogram.png',linestyle = [':','--','-.','-']): # plot gls for RV and O-C, only if they exist
     
        plot_gls_warnings=Warning_log([],'Plotting gls')
        # first let's check if we have RV and O-C data
        if (len(self.rvs)==0):
            plot_gls_warnings.update_warning_list('No RV data to run gls, gls will not be run!')
        elif (len(self.o_c)==0):
            plot_gls_warnings.update_warning_list('No O-C data to run gls, gls will not be run!')
        else:
            # building objects to run gls on
            rv_data=signal_data(self.jd,self.rvs,self.rv_err,sig_for_gls=sig)
            oc_data=signal_data(self.jd,self.o_c,self.rv_err,sig_for_gls=sig)
            rv_data.gls(ind_sig=ind_sig)
            oc_data.gls(ind_sig=ind_sig)
            
            # preparing plot properties
            plt.rc('text',usetex= True)
            font = {'family' : 'serif','weight' : 'bold','size'   : 18,'serif':['Helvetica']}
            plt.rc('font', **font)    
            plt.figure(figsize=(8,6.5))
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1.5]) 
            gs.update(wspace=0.15)
            plt.subplots_adjust(hspace=0.005)
            ax1 = plt.subplot(gs[0,0]) # in ax1 we will plot the RV signal periodogram
            ax2 = plt.subplot(gs[1,0]) # in ax2 we will plot the O-C periodogram     

            ############### Ploting the best peaks ##########################         
            # loop to generate annotations next to highest peaks
            for l in range(min(3,rv_data.number_of_significant_peaks)): # highlight 3 highest peaks, unless fewer pass the significance criteria 
                ax1.annotate(r'%.2f d'%rv_data.best_per[l],(rv_data.best_per[l],rv_data.best_peaks[l]*0.95), fontsize=16, color = 'k')
  
            ax1.plot(rv_data.periods, rv_data.P_G, '-',color ='r', lw=1.3, label='generalized')
            ax2.plot(oc_data.periods, oc_data.P_G, '-',color ='r', lw=1.3, label='generalized')
 
            # limits for abssica axes, lower bound is always 1 day     
            ax1.set_xlim(1, rv_data.gls_range)
            ax2.set_xlim(1, rv_data.gls_range) 

            ############ plot the significance lines. ############

            j = 0 # iterator to go through the linestyle list

            xlim = (rv_data.periods[0], rv_data.periods[-1])

            for zi in rv_data.z:
                ax1.plot(xlim, (zi, zi), linestyle[j], color = '0.25', lw=1)
                j = j +1

            # logscale for abscissa axis only
            ax1.semilogx()
            ax2.semilogx() 
               
            ax2.set_xlabel('days', fontsize=22, color = 'k')
            ax1.set_ylabel(r'Power', fontsize=22, color = 'k')
          
            plt.savefig(filetosave, format='png', dpi=600 ) 
            plot_gls_warnings.print_warning_log()           
            return
        
                                     
              
    def plot_periodogram(self, filetosave='periodogram.png'):
        if(self.model_saved): # if we have the model we can run plot_gls on an appropriate kernel
            self.fit_results.rv_model.plot_gls(filetosave=filetosave)
        else: # otherwise we can only plot it for initial data
            object_for_gls=rvmodel(self.filelist.time,self.filelist.rvs,self.filelist.rv_err,[0.0]*len(self.filelist.time))
            object_for_gls.plot_gls(filetosave=filetosave)
        return
        
    def plot_time_series(self, filetosave='time_series.png', plot_title='RV time series', save = True, dpi=200):

        '''Generate a kernel object and plot it's time series'''
        
        if not (self.fit_performed): # to update in case we hadn't performed a fit, but added new datasets         
            self.fit_results=kernel(summary(parameters(0,0,0,0,0),parameter_errors([0],[0],[0],0,0)), self.filelist.time,self.filelist.rvs,self.filelist.rv_err,np.zeros(len(self.filelist.time)),np.zeros(len(self.filelist.time)),self.filelist.time,0,0,0,self.filelist.idset,0,0,0,0,0,0)
 
        self.fit_results.plot_ker(filetosave=filetosave, plot_title=plot_title, save=save, dpi=dpi)
        return
                   
    def plot_signal_one_planet(self,planet,filetosave='one_planet_signal.png', plot_title='One planet signal', save = True, dpi=200, warnings=None):
    
        '''Plots a signal for a chosen planet, understood as the overall signal minus contribution from all other planets'''     
        print_at_end=False

        if (warnings==None): # this maens we are running this function stand alone, rather than inside plot_signal_each_planet, so we need to create a new warnings object and print warnings at the end of this function
            warnings=Warning_log([],'Plotting one planet signal')
            print_at_end=True

        if (planet>self.npl):
            warnings.update_warning_list('Planet number %d higher than number of planets %d, cannot generate plot!'%(planet,self.npl))
        else:                 
            #special_task_copy=copy.deepcopy(self) # create new signal_fit object, which will be the same except semiamplitude for chosen planet will be zero
            #normal_copy=copy.deepcopy(self) # create new signal_fit object, which will be the same  
            special_task_copy=dill.copy(self)
            normal_copy=dill.copy(self)
       
            flag=normal_copy.fitting(fileinput=True, filename='temp_test', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,1],return_flag=True) # this should find signal from other planets
            if not (flag==1):
                warnings.update_warning_list('Failed to find signal from all planets!')
            else:
                special_task_copy.params.update_K(planet,0.0)    
                flag=special_task_copy.fitting(fileinput=True, filename='temp_test', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,1],return_flag=True) # this should find signal from other planets
                if not (flag==1):
                    warnings.update_warning_list('Failed to find signal from all planets except planet %d.'%planet)
                else:
                    if not (len(self.fit_results.rv_model.rvs)==len(special_task_copy.fit_results.rv_model.rvs)):
                        warnings.update_warning_list('Number of data points is %d, but after removing %d planet the number is different, i. e. %d. Something must have gone wrong.'%(len(normal_copy.fit_results.rv_model.rvs),planet,len(special_task_copy.fit_results.rv_model.rvs)))
                    rvs_to_pass=[]        
                    for i in range(min(len(normal_copy.fit_results.rv_model.rvs),len(special_task_copy.fit_results.rv_model.rvs))):
                        rvs_to_pass.append(normal_copy.fit_results.rv_model.rvs[i]-special_task_copy.fit_results.rv_model.rvs[i])
                    rvs_to_pass=np.array(rvs_to_pass)
                    if not (len(normal_copy.fit_results.model)==len(special_task_copy.fit_results.model)):
                        warnings.update_warning_list('Number of points in model is %d, but after removing %d planet the number is different, i. e. %d. Something must have gone wrong.'%(len(normal_copy.fit_results.model),planet,len(special_task_copy.fit_results.model)))
                    model_to_pass=[]  
                    for i in range(min(len(normal_copy.fit_results.model),len(special_task_copy.fit_results.model))):          
                        model_to_pass.append(float(normal_copy.fit_results.model[i])-float(special_task_copy.fit_results.model[i]))
                    model_to_pass=np.array(model_to_pass)
                    special_task_kernel=kernel(normal_copy.fit_results.stat, normal_copy.fit_results.rv_model.jd,rvs_to_pass,normal_copy.fit_results.rv_model.rv_err,normal_copy.fit_results.rv_model.o_c,model_to_pass,normal_copy.fit_results.model_jd,normal_copy.npl,0.0,0.0,normal_copy.fit_results.idset,normal_copy.fit_results.stat_array_saved,0.0,0.0,0.0,0.0,0)
                    special_task_kernel.plot_ker(filetosave=filetosave, plot_title=plot_title, save=save, dpi=dpi, plot_oc=False)
        
        if (print_at_end):
            warnings.print_warning_log()
        
        return
       
    def plot_signal_each_planet(self,filestosave=[],plot_titles=[],dpi=200):
        
        planetsignalplotwarnings=Warning_log([],'Plotting individual planet signals')       
        
        if (len(filestosave)<self.npl): # this means the user didn't provide custom filenames for all the plots
            if (len(filestosave)>0): # so he wanted to provide custom names, but didn't provide for all
                planetsignalplotwarnings.update_warning_list('Too few file names for the plots, will assume default for remaining!')
            for i in np.arange(len(filestosave),self.npl,1):
                filestosave.append('signal_plot_pl_%d.png'%(i+1))
            filestosave=np.array(filestosave)        
        if (len(filestosave)>self.npl):
            planetsignalplotwarnings.update_warning_list('Too many file names for the plots, additional will be ignored!')
            filestosave=filestosave[:self.npl]
        if (len(plot_titles)<self.npl): # this means the user didn't provide custom plottitles for all the plots
            if (len(plot_titles)>0): # so he wanted to provide custom names, but didn't provide for all
                planetsignalplotwarnings.update_warning_list('Too few plot titles, will assume default for remaining!')
            for i in np.arange(len(plot_titles),self.npl,1):
                plot_titles.append('signal_pl_%d'%i)
            plot_titles=np.array(plot_titles)        
        if (len(plot_titles)>self.npl):
            planetsignalplotwarnings.update_warning_list('Too many plot tiltes, additional will be ignored!')
            plot_titles=plot_titles[:self.npl]
        
        # Done with the warnings, now let's plot
        for i in range(self.npl):
            self.plot_signal_one_planet(i,filetosave=filestosave[i], plot_title=plot_titles[i], save = True, dpi=dpi, warnings=planetsignalplotwarnings)
                    
        planetsignalplotwarnings.print_warning_log()
        return         
        
        
