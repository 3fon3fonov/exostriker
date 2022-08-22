#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import matplotlib.pyplot as plt


# In[2]:


#calculating thiele constants
def thiele(a,omega,Omega,i):
    A=a*(np.cos(omega)*np.cos(Omega)-np.sin(omega)*np.sin(Omega)*np.cos(i))
    B=a*(np.cos(omega)*np.sin(Omega)+np.sin(omega)*np.cos(Omega)*np.cos(i))
    F=a*(-np.sin(omega)*np.cos(Omega)-np.cos(omega)*np.sin(Omega)*np.cos(i))
    G=a*(-np.sin(omega)*np.sin(Omega)+np.cos(omega)*np.cos(Omega)*np.cos(i))
    return A,B,F,G


# in the following i constructed 2 methods for calculating E, the first one is better we can talk about it

# In[10]:


def calc_E(e,M,n=30):
    final=[]
    for j in M:
        if e<=0.8:
            E=j
        if e>0.8:
            E=np.pi
        values=np.ones(n)
        values[0]=E
        i=0
        while i<=n:
            i+=1
            values[i]=values[i-1]-(values[i-1]-e*np.sin(values[i-1])-j)/(1-e*np.cos(values[i-1]))
            if abs(values[i]-values[i-1])<1e-8:
                values[-1]=values[i]
                break
            if i==n:
                return print("error in calculating E")
        final.append(values[-1])
    
    final=np.array(final)
        
        
    return final

#Second method, u can ignore it for now:

# #calculating E via newton method, input M array, output E array
# def calc_E2(e,M,n=10):
#     if e<=0.8:
#         E=M
#     if e>0.8:
#         E=np.pi*np.ones(len(M))
#     #The next line constructs an 2d array with starting values E for every M and length according to n.
#     values=np.vstack([E,np.ones((n,len(M)))])
#     i=1
#     while i<n+1:
#         values[i]=values[i-1]-(values[i-1]-e*np.sin(values[i-1])-M)/(1-e*np.cos(values[i-1]))
#         i=i+1
#     return values[n]

# In[4]:


def coords(P,e,om,i,Om,T0,a,t):
    
    #calculating thiele constants
    const=thiele(a,om,Om,i)
    
    #calculating E
    M=2*np.pi*((t-T0)%P)/P
    E=calc_E(e,M)

    #calculating elliptical rectangular coords from E and e
    X=np.cos(E)-e
    Y=(1-e**2)**0.5 *np.sin(E)

    #with thiele, X,Y, we can calculate the final position x,y for a time t
    x=const[0]*X+const[2]*Y
    y=const[1]*X+const[3]*Y
    
    return x,y


# In[11]:


#data given for X or Y
def residual(model,data):
    #residuals
    res=model-data
    return res


# In[6]:


def loglikelihood_am(model,data,err,s): #s=jitter
    #first term just depends on number of obs:
    term_1=-np.log(2*np.pi)*len(data)/2
    
    #term 2 and term 3 are summation
    
    term_2=-0.5*sum(np.log(err**2 +s**2))
    #residual function from before just calculates the diff between model and obs
    res=residual(model,data)
    term_3=-0.5*sum((res**2) /(err**2 +s**2))
        
        
    ln_L=term_1-0.5*(term_2+term_3)
    
    return ln_L
        
        
        
        
        


# In[7]:


#reminder: the data,errors and the times must be given in arrays
def final(x_data,x_err,s_x,y_data,y_err,s_y,P,e,om,i,Om,T0,a,t):
    
    x_mod,y_mod=coords(P,e,om,i,Om,T0,a,t)
    #calculating the likelihoods for both dimensions
    L_x=loglikelihood_am(x_mod,x_data,x_err,s_x)
    L_y=loglikelihood_am(y_mod,y_data,y_err,s_y)
    #adding them togther because its allowed
    L=L_x+L_y
    
    return L,x_mod,y_mod


# #very simple case example
x_data=np.array([1,-0.5,-0.5,1])
x_err=np.ones(4)*0.1
y_data=np.array([0,0.86,-0.86,0])
y_err=np.ones(4)*0.1
t=np.linspace(0,365,4)
data_am = final(x_data,x_err,0,y_data,y_err,0,365,0,0,math.radians(45),0,0,1,t)

# #model example
t=np.linspace(0,365,10000)
model_am = coords(365,0,0,math.radians(45),0,0,1,t)

plt.plot(model_am[0],model_am[1], 'r-')
plt.plot(data_am[1],data_am[2], '*')
plt.show()
# In[ ]:




