#!/usr/bin/python3


__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
import numpy as np
#import emcee
import emcee_ES as emcee


class CustomSampler(emcee.EnsembleSampler):
 
    
    def unique_rows(self):

        '''
        Given an array, remove identical rows and also sort it
        '''

        # Perform lex sort and get sorted data
        sorted_idx = np.lexsort(self.flatchain.T)
        sorted_data =  self.flatchain[sorted_idx,:]
        sorted_lnL  =  self.flatlnprobability[sorted_idx]
        # Get unique row mask
        row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))

        # Get unique rows 
        self.samples = sorted_data[row_mask] 
        self.lnL     = sorted_lnL[row_mask]
       # sorted_lnL =  lnL[sorted_idx]
        #self.lnL = sorted_lnL[row_mask]

        lnL_max_idx = np.argmax(self.lnL)
        #print(abs(lnL_min))

        # get samples at minimum Lnl
        self.maxlnL = self.samples[lnL_max_idx]
        self.lnL_max = self.lnL[lnL_max_idx]

        return


    def correct_rows(self,f,ndset,npl,hkl):

        '''Corrects angles and eccentricities for all samples'''

        i=0
        self.means=np.zeros(len(f))
        for k in range(len(f)):
            idx=f[k]
            if (idx<2*ndset):
                self.means[i]=np.mean(self.samples[:,i])
            elif (idx<2*ndset+7*npl):
                nr=idx-2*ndset
                #x=int(nr/7)
                if (np.mod(nr,7)<2):
                    self.means[i]=np.mean(self.samples[:,i]) 
                elif (np.mod(nr,7)==2): # correct eccentricities
                    #for j in range(len(self.samples)):
                       # if (self.samples[j,i]<0):
                       #     self.samples[j,i]=abs(self.samples[j,i])
                            #if(f[k+1]==i+1):
                            #    self.samples[j,i+1]=self.samples[j,i+1]+180.0
                            #if(f[k+2]==i+2):
                            #    self.samples[j,i+2]=self.samples[j,i+2]+-180.0
                    self.means[i]=np.mean(self.samples[:,i])
                elif (np.mod(nr,7)==3): # correct w to be in a 360 interval around mean value 
                    
                    if hkl == False:
                        self.samples[:, i] = self.samples[:,i]%360
                        meanw=self.circ_mean_np(self.samples[:,i])  
    
                        for j in range(len(self.samples)):
                            self.samples[j,i]=np.where(self.samples[j,i]<meanw-180.0,self.samples[j,i]+360.0,self.samples[j,i])
                            self.samples[j,i]=np.where(self.samples[j,i]>meanw+180.0,self.samples[j,i]-360.0,self.samples[j,i])
    
                        self.means[i]=meanw 
                        
                    else:
                        self.means[i]=np.mean(self.samples[:,i]) 
                        
                elif (np.mod(nr,7)==4):# correct M to be in a 360 interval around mean value
                    self.samples[:, i] = self.samples[:,i]%360
                    meanw=self.circ_mean_np(self.samples[:,i])  
                    for j in range(len(self.samples)):
                        self.samples[j,i]=np.where(self.samples[j,i]<meanw-180.0,self.samples[j,i]+360.0,self.samples[j,i])
                        self.samples[j,i]=np.where(self.samples[j,i]>meanw+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                        
                    self.means[i]=meanw 
            elif (idx<2*ndset+6*npl):# correct i to be in a 180 interval around mean value
                self.means[i]=np.mean(self.samples[:,i])
                meani=self.means[i]
 
            elif (idx<2*ndset+7*npl):# correct lineofnodes to be in a 360 interval around mean value 
                self.samples[:, i] = self.samples[:,i]%360
                meanw=self.circ_mean_np(self.samples[:,i])  
                for j in range(len(self.samples)):
                    self.samples[j,i]=np.where(self.samples[j,i]<meanw-180.0,self.samples[j,i]+360.0,self.samples[j,i])
                    self.samples[j,i]=np.where(self.samples[j,i]>meanw+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                    
                self.means[i]=meanw 
            else:
                self.means[i]=np.mean(self.samples[:,i])
            i=i+1
        return

    def circ_mean_np(self, angles,azimuth=True):  
        """find circular mean"""  
        rads = np.radians(angles)  
        av_sin = np.mean(np.sin(rads)) 
        av_cos = np.mean(np.cos(rads))  
        ang_rad = np.arctan2(av_sin,av_cos)  
        ang_deg = np.degrees(ang_rad)  
        if azimuth:  
            ang_deg = np.mod(ang_deg,360.)  
        return  ang_deg  

    def get_meadians(self,f):
        '''Gets median values in posterior sample'''
        self.median=np.zeros(len(f))
        for k in range(len(f)):
            self.median[k]=np.median(self.samples[:,k])
        return     

                
    def save_samples(self,f,ndset,npl,hkl):
        self.unique_rows()
        self.correct_rows(f,ndset,npl,hkl)
        self.get_meadians(f)
        return
