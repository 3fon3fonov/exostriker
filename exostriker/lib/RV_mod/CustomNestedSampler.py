#!/usr/bin/python


__author__ = 'Trifon Trifonov'

import sys 
#sys.path.insert(0, '../lib')
import numpy as np
#from dynesty import NestedSampler, DynamicNestedSampler
import dynesty


class CustomNestedSampler(dynesty.NestedSampler,dynesty.DynamicNestedSampler):

    
    def __init__(self, *args, **kwargs):
        super(CustomNestedSampler, self).__init__(*args, **kwargs)   
    
    def convert_to_samples(self):
        self.samples = self.results.samples
    
    def unique_rows(self):
    
        '''
        Given an array, remove identical rows and also sort it
        '''
    
        
    
    
        # Perform lex sort and get sorted data
        sorted_idx = np.lexsort(self.flatchain.T)
        sorted_data =  self.flatchain[sorted_idx,:]

        # Get unique row mask
        row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))

        # Get unique rows 
        out = sorted_data[row_mask] 
        self.results.samples = out
        
        
        # same for LnL
        lnL = np.hstack(self.lnprobability)       
        sorted_lnL =  lnL[sorted_idx]   
        self.lnL = sorted_lnL[row_mask]        
        
        
        lnL_min_idx = np.argmax(self.lnL)
        #print(abs(lnL_min))
        
        # get samples at minimum Lnl
        self.minlnL = out[lnL_min_idx]
        self.lnL_min = self.lnL[lnL_min_idx]

        return


    def correct_rows(self,f,ndset,npl):

        '''Corrects angles and eccentricities for all samples'''

        i=0
        self.means=np.zeros(len(f))
        for k in range(len(f)):
            idx=f[k]
            if (idx<2*ndset):
                self.means[i]=np.mean(self.results.samples[:,i])
            elif (idx<2*ndset+7*npl):
                nr=idx-2*ndset
                #x=int(nr/7)
                if (np.mod(nr,7)<2):
                    self.means[i]=np.mean(self.results.samples[:,i]) 
                elif (np.mod(nr,7)==2): # correct eccentricities
                    #for j in range(len(self.results.samples)):
                       # if (self.results.samples[j,i]<0):
                       #     self.results.samples[j,i]=abs(self.results.samples[j,i])
                            #if(f[k+1]==i+1):
                            #    self.results.samples[j,i+1]=self.results.samples[j,i+1]+180.0
                            #if(f[k+2]==i+2):
                            #    self.results.samples[j,i+2]=self.results.samples[j,i+2]+-180.0
                    self.means[i]=np.mean(self.results.samples[:,i])
                elif (np.mod(nr,7)==3): # correct w to be in a 360 interval around mean value 
                    self.results.samples[:, i] = self.results.samples[:,i]%360
                    meanw=self.circ_mean_np(self.results.samples[:,i])  

                    for j in range(len(self.results.samples)):
                        self.results.samples[j,i]=np.where(self.results.samples[j,i]<meanw-180.0,self.results.samples[j,i]+360.0,self.results.samples[j,i])
                        self.results.samples[j,i]=np.where(self.results.samples[j,i]>meanw+180.0,self.results.samples[j,i]-360.0,self.results.samples[j,i])

                    self.means[i]=meanw 
                elif (np.mod(nr,7)==4):# correct M to be in a 360 interval around mean value
                    self.results.samples[:, i] = self.results.samples[:,i]%360
                    meanw=self.circ_mean_np(self.results.samples[:,i])  
                    for j in range(len(self.results.samples)):
                        self.results.samples[j,i]=np.where(self.results.samples[j,i]<meanw-180.0,self.results.samples[j,i]+360.0,self.results.samples[j,i])
                        self.results.samples[j,i]=np.where(self.results.samples[j,i]>meanw+180.0,self.results.samples[j,i]-360.0,self.results.samples[j,i])
                        
                    self.means[i]=meanw 
            elif (idx<2*ndset+6*npl):# correct i to be in a 180 interval around mean value
                self.means[i]=np.mean(self.results.samples[:,i])
                meani=self.means[i]
 
            elif (idx<2*ndset+7*npl):# correct lineofnodes to be in a 360 interval around mean value 
                self.results.samples[:, i] = self.results.samples[:,i]%360
                meanw=self.circ_mean_np(self.results.samples[:,i])  
                for j in range(len(self.results.samples)):
                    self.results.samples[j,i]=np.where(self.results.samples[j,i]<meanw-180.0,self.results.samples[j,i]+360.0,self.results.samples[j,i])
                    self.results.samples[j,i]=np.where(self.results.samples[j,i]>meanw+180.0,self.results.samples[j,i]-360.0,self.results.samples[j,i])
                    
                self.means[i]=meanw 
            else:
                self.means[i]=np.mean(self.results.samples[:,i])
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
            self.median[k]=np.median(self.results.samples[:,k])
        return     

                
    def save_samples(self,f,ndset,npl):
        #self.unique_rows()
       # self.convert_to_samples()
        self.correct_rows(f,ndset,npl)
        self.get_meadians(f)
        return
