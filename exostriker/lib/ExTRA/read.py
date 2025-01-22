
import numpy as np
from .hipparcos import hip_JD


#reading out RV data

def RV_read(path):
    with open(str(path),"r") as f:
        rawdata=f.readlines()
    data=[]
    for line in rawdata:
        line=line.strip()
        line=line.split(" ")
        line=' '.join(line).split()
        line=line[:3]
        line=np.array(line)
        data.append(line)
    data=np.transpose(data)
    data=data.astype("float64")
    return data

#ordering RV data the way i want: t,data,err

def RV_order(RV_data):
    
    t=[]
    dat=[]
    err=[]
    #putting data in order:
    t_ind=list(range(0,len(RV_data)))
    for i in t_ind:
        t.append(RV_data[i][0])
        dat.append(RV_data[i][1])
        err.append(RV_data[i][2])
    return t,dat,err
################

#Reading out HIP data:

def hip_read(path):
    """Returns HIP astrometric data and time for hip measurements in a array"""
    with open(str(path),"r") as g:
        next(g)
        rows=[]
        lines=g.readlines()
        for line in lines:
            rows.append(line.split())
    i=0
    while i<len(rows):
        rows[i]=rows[i][1:7]
        rows[i][0]=str(float(rows[i][0]))#+2400000.5)
        i=i+1
    HIP2=np.transpose(rows)
    HIP=HIP2.astype("float64")
    HIP_epochs,A_5,A_3,A_4,A_8,A_9=HIP #A3=cos--x , #A4=sin ---y


    ################
    #HIP_epochs,A_5,A_3,A_4,A_8,A_9=HIP2 
    ################



    A_6=A_3*(HIP_epochs)
    A_7=A_4*(HIP_epochs)


    hip_ad=np.array([A_3,A_4,A_5,A_6,A_7,A_8,A_9])

    t_HIP=hip_JD(hip_ad)

    return hip_ad,t_HIP

