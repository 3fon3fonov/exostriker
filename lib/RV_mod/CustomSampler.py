#!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
#sys.path.insert(0, '../lib')
import numpy as np
import emcee



class CustomSampler(emcee.EnsembleSampler):

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
        self.samples = out

        # same for LnL
        lnL = np.hstack(self.lnprobability)
        sorted_lnL =  lnL[sorted_idx]
        self.lnL = sorted_lnL[row_mask]

        lnL_max_idx = np.argmax(self.lnL)
        #print(abs(lnL_min))

        # get samples at minimum Lnl
        self.maxlnL = out[lnL_max_idx]
        self.lnL_max = self.lnL[lnL_max_idx]

        return

    def correct_rows(self,f,ndset,npl):

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
                    self.means[i]=np.mean(self.samples[:,i]) 
                    meanw=self.means[i]
                   # for j in range(len(self.samples)):
                   #     self.samples[j,i]=np.where(self.samples[j,i]<meanw-180.0,self.samples[j,i]+360.0,self.samples[j,i])
                   #     self.samples[j,i]=np.where(self.samples[j,i]>meanw+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                    # now let's make sure meanw is between 0 and 360:
                    newmeanw=np.fmod(meanw,360.0)
                    delta=newmeanw-meanw
                    if not (delta==0):
                        for j in range(len(self.samples)):
                            self.samples[j,i]=self.samples[j,i]+delta
                elif (np.mod(nr,7)==4):# correct M to be in a 360 interval around mean value
                    self.means[i]=np.mean(self.samples[:,i]) 
                    meanM=self.means[i]
                   # for j in range(len(self.samples)):
                   #     self.samples[j,i]=np.where(self.samples[j,i]<meanM-180.0,self.samples[j,i]+360.0,self.samples[j,i])
                   #     self.samples[j,i]=np.where(self.samples[j,i]>meanM+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                    # now let's make sure meanw is between 0 and 360:
                    newmeanM=np.fmod(meanM,360.0)
                    delta=newmeanM-meanM
                    if not (delta==0):
                        for j in range(len(self.samples)):
                            self.samples[j,i]=self.samples[j,i]+delta
            elif (idx<2*ndset+6*npl):# correct i to be in a 180 interval around mean value
                self.means[i]=np.mean(self.samples[:,i])
                meani=self.means[i]
               # for j in range(len(self.samples)):
               #     self.samples[j,i]=np.where(self.samples[j,i]<meani-90.0,self.samples[j,i]+180.0,self.samples[j,i])
               #     self.samples[j,i]=np.where(self.samples[j,i]>meani+90.0,self.samples[j,i]-180.0,self.samples[j,i])
                # now let's make sure meani is between 0 and 180:
                newmeani=np.fmod(meani,180.0)
                delta=newmeani-meani
                if not (delta==0):
                    for j in range(len(self.samples)):    
                        self.samples[j,i]=self.samples[j,i]+delta
            elif (idx<2*ndset+7*npl):# correct lineofnodes to be in a 360 interval around mean value 
                self.means[i]=np.mean(self.samples[:,i])
                meancap=self.means[i]
              #  for j in range(len(self.samples)):
               #     self.samples[j,i]=np.where(self.samples[j,i]<meancap-180.0,self.samples[j,i]+360.0,self.samples[j,i])
               #     self.samples[j,i]=np.where(self.samples[j,i]>meancap+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                # now let's make sure meancap is between 0 and 360:
                newmeancap=np.fmod(meancap,360.0)
                delta=newmeancap-meancap
                if not (delta==0):
                    for j in range(len(self.samples)):
                        self.samples[j,i]=self.samples[j,i]+delta
            else:
                self.means[i]=np.mean(self.samples[:,i])
            i=i+1
        return


    def get_meadians(self,f):
        '''Gets median values in posterior sample'''
        self.median=np.zeros(len(f))
        for k in range(len(f)):
            self.median[k]=np.median(self.samples[:,k])
        return     

                
    def save_samples(self,f,ndset,npl):
        self.unique_rows()
        self.correct_rows(f,ndset,npl)
        self.get_meadians(f)
        return
