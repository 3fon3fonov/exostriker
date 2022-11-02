 #!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'


#import sys 
#sys.path.insert(0, '../lib')
import numpy as np

from .functions import *  
from .errors import Error, InputError, FittingError
from .Warning_log import Warning_log
     
# classes for processing RV files
NDSETMAX=20


class rvfile(object):
    
    def __init__(self,name,path):
        # let's make sure all the files are sorted, it won't take long and it will be needed for mean_value and first_observation functions
        #command = 'echo "$(sort %s)" > %s'%(path,path) # command to sort the file in place
        #text,flag=run_command_with_timeout(command, 3) # run the command  
        # now save the properties
        #path =  copy_file_to_datafiles(path)

        self.name=name
        self.path=path
        self.reading_in_progress=open(str(self.path),'r') # to be able to read file line by line using read_one_line function
        self.reading_in_progress.close() # for now let's keep it closed, though  

    # function to read just one line
    def read_one_line(self,offset=0.0): 
        if (self.reading_in_progress.closed): # if the file is currently closed we treat it the same way as end of file
            isline=False
            comment=False
            x=0
            y=0
            y_error=0        
        else:
            line=self.reading_in_progress.readline()
            if not line: # end of file, only isline value matters in this case but we need to declare all for output
                isline=False
                comment=False
                x=0
                y=0
                y_error=0
            else:
                isline=True
                if line[0] == '#': #skip comments, in this case x, y, y_error values are not important, just needed for function output
                    comment=True
                    x=0
                    y=0
                    y_error=0
                else:
                    line = line.strip()
                    #line = line.split()
                    line = [float(x) for x in line.split()]

                    if is_float(line[1]): #omit lines with invalid RV data, for example 'NaN', '-Inf' etc.
                        comment=False
                        x=float(line[0])
                        y=float(line[1])-offset # subtracting the offset at this point already
                        y_error=float(line[2])    
                    else: # lines with invalid RV data will also be treated as comments
                        comment=True
                        x=0
                        y=0
                        y_error=0
        return [x,y,y_error,isline,comment]

    # numerical integration of RV data to find the mean value as a first proxy for offset
    def mean_value(self):
        self.reading_in_progress=open(self.path,'r') # to be able to read file line by line using read_one_line function
        comment = True
        # In case there is only one data point we want to return it's abscissa, so we will add a tiny rectangle around this point to our integral, but we want it to be tiny, so it normally won't affect the results 
        dt=0.00001
        while(comment):
            x0,y0,y_error,isline,comment = self.read_one_line()
            if not (isline):
                raise InputError('RV file %s contains no data!'%self.name)
            else: # S will be our integral, t0 is going to be the beginning of the interval over which we calculate our integral
                S=y0*dt
                t0=x0-dt              
        while 1: # now we have the starting point, let's go through the rest of the file
            x,y,y_error,isline,comment = self.read_one_line()
            if not (isline):
                break # the file is over, we finish
            elif (comment):
                continue
            else:
                S=S+(x-x0)*(y+y0)/2 # add a trapezoid to the integral 
                x0=x
                y0=y              
        self.reading_in_progress.close()
        return S/(x0-t0) # divide integral by interval length

    # time of first observation in a given dataset
    def first_datapoint(self):
        self.reading_in_progress=open(self.path,'r') # to be able to read file line by line using read_one_line function
        comment = True
        while(comment):
            x0,y0,y_error,isline,comment = self.read_one_line()
            if not (isline):
                raise InputError('RV file %s contains no data!'%self.name)
        self.reading_in_progress.close()
        return x0            

class rvfile_list(object): # this will store a list of rvfile objects
    def __init__(self,ndset,names,paths): # in initialization we use default offsets and jitters

        if (ndset>NDSETMAX):
            raise InputError('Too many data sets! Maximum number is %d'%NDSETMAX)
        self.ndset=ndset
        self.files=[]
        self.time=[]
        self.rvs=[]
        self.rv_err=[]
        self.idset=[]
        for i in range(self.ndset):
            self.files.append(rvfile(names[i],paths[i])) # creating a rvfile object for each data set
        
         
    def add_datafile(self,name,path): # add another data file to the list
        #path =  copy_file_to_datafiles(path)

        if (os.path.isfile(path)):
            self.files.append(rvfile(name,path))
            self.ndset=self.ndset+1    
            flag=1     
        else:
            warnings=Warning_log(['Path %s does not correspond to a valid file! Request ignored'%path], 'Adding a dataset %s'%name)
            warnings.print_warning_log()
            flag=0
        return flag
         
    def remove_datafile(self,number): # add another data file to the list

        if not (number<self.ndset):
            warnings=Warning_log(['Dataset index outside of range'], 'Removing a dataset %d'%number)
            warnings.print_warning_log()
            flag=0
        else:
            self.ndset-=1
            self.files.pop(number)
            self.idset  = self.idset[np.where(self.idset != number)]
            self.time   = self.time[np.where(self.idset != number)]
            self.rvs    = self.rvs[np.where(self.idset != number)]
            self.rv_err = self.rv_err[np.where(self.idset != number)] 
            self.idset[self.idset > number] -= 1  
            flag=1     
 
        return flag     
         
    def first_observation (self): # search for the earliest observation
        if (len(self.files)==0): # shouldn't happen, but just in case
            raise Exception ('No RV files provided!')
        else: 
            # each file is by now sorted, so the first line of the file is it's earliest Julian Day
            firstone = self.files[0].first_datapoint()
            for i in range(1,len(self.files)): # search all datafiles for earliest JD of all
                anotherone = self.files[i].first_datapoint() # to compare with current minimum
                if (anotherone<firstone): 
                    firstone=anotherone
        return firstone
        
   # def read_rvfiles(self,offsets,justthenewone=False):
        
        
    #    x =
   #     y = 
    #    y_error = 
    #    data_set =       
   


    def read_rvfile(self,path,  ndset,justthenewone=False):

        x       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
        y       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
        y_error = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
 
        data_set = np.array([ndset-1]*len(x))
       

        # saving sorted arrays to appropriate object attributes        
        self.time = np.concatenate((self.time,x))
        self.rvs = np.concatenate((self.rvs,y))
        self.rv_err = np.concatenate((self.rv_err,y_error))
        #self.idset = data_set
        
        dataset_ind = np.concatenate((self.idset,data_set)) #there is a problem with np.concatenate() it always returns floats not integers
        self.idset = dataset_ind.astype(int) # this is a quick fix
        #print(self.idset)

        #sorting data by time
        ind_sort = np.argsort(self.time)
        self.time = self.time[ind_sort]
        self.rvs = self.rvs[ind_sort]
        self.rv_err = self.rv_err[ind_sort]
        self.idset = self.idset[ind_sort]
        return 

    def read_rvfiles(self,offsets,justthenewone=False):
        i = 0
        x = []
        y = []
        y_error = []
        data_set = []
        if (justthenewone): # this means we had already read rvfiles, but we're calling a function to read a new freshly added one
            i=0#self.ndset-1 #not working!!!!
        else:
            i=0
        while i<self.ndset:
            self.files[i].reading_in_progress=open(self.files[i].path,'r')
            while 1:
                xx, yy, yy_error, isline, comment = self.files[i].read_one_line(offsets[i])
                if not (isline):
                    break
                elif (comment):
                    continue
                else:
                    x.append(xx)
                    y.append(yy)
                    y_error.append(yy_error)
                    data_set.append(i)       
            self.files[i].reading_in_progress.close()   
                   
            i += 1
        #casting all lists on arrays which work faster 
        x = np.array(x)
        y = np.array(y)
        y_error = np.array(y_error)
        data_set = np.array(data_set)
       

        # saving sorted arrays to appropriate object attributes        
        self.time = x #np.concatenate((self.time,x))
        self.rvs = y#np.concatenate((self.rvs,y))
        self.rv_err = y_error#np.concatenate((self.rv_err,y_error))
        self.idset = data_set
        
        dataset_ind = np.concatenate((self.idset,data_set)) #there is a problem with np.concatenate() it always returns floats not integers
        self.idset = dataset_ind.astype(int) # this is a quick fix
        #print(self.idset)

        #sorting data by time
        ind_sort = np.argsort(self.time)
        self.time = self.time[ind_sort]
        self.rvs = self.rvs[ind_sort]
        self.rv_err = self.rv_err[ind_sort]
        self.idset = self.idset[ind_sort]
        return 
        
