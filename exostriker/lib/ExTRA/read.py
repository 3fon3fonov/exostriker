
import numpy as np
from astropy.table import Table
import pandas as pd
from .hipparcos import hip_JD

__all__ = ["RV_read", "RV_order", "hip_data", "hip_read","hip_init"]


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
def hip_data(location,hip_id,format=""):
    
    
    if hip_id[:3]=="HIP":
        number=hip_id[3:]
    else:
        number=hip_id
    

    if len(number)==6:
        id_full="H"+str(number)
    if len(number)==5:
        id_full="H0"+str(number)
    if len(number)==4:
        id_full="H00"+str(number)
    if len(number)==3:
        id_full="H000"+str(number)
    if len(number)==2:
        id_full="H0000"+str(number)
    if len(number)==1:
        id_full="H00000"+str(number)




    destination=location+id_full[:4]+"/"+id_full+".d"




    cols = ["IORB", "EPOCH", "PARF", "CPSI", "SPSI", "RES", "SRES"]

    if format=="pandas":
        df = pd.read_csv(
        destination,
        comment="#",
        sep='\s+',
        names=cols)

    else:
        df = Table.read(
        destination,
        format="csv",
        comment="#",
        delimiter=" ",
        names=cols)

    

    header_lines = []

    with open(destination) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line.strip())
            else:
                break  # stop when actual data starts

        header=[]
        for line in header_lines[5:-2]:
            line=line.split()[1:]
            header.append(line)


    return df,header


def hip_init(table_raw,timeformat="jd"):
    """Returns HIP astrometric data and time for hip measurements in a array"""
    number,HIP_epochs,A_5,A_3,A_4,A_8,A_9=table_raw #A3=cos--x , #A4=sin ---y
    A_6=A_3*(HIP_epochs)
    A_7=A_4*(HIP_epochs)


    hip_ad=np.array([A_3,A_4,A_5,A_6,A_7,A_8,A_9])

    t_HIP=hip_JD(hip_ad,format=timeformat)

    return hip_ad,t_HIP


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

