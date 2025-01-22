import numpy as np
from .RVcomb import *
from .RVsolo import *
from .hipparcos import *
from .astrometry import *


def loglikelihood(res,err,s): #s=jitter
    '''
    Computes the negative loglikelihood
    
    Parameters:
    ----------
    res: array
        residual, observed minus computed
    err: array
        error of data
    s: float
        jitter
    '''
    #first term just depends on number of obs:
    term_1=-np.log(2*np.pi)*len(res)/2
    
    #term 2 and term 3 are summation
    
    term_2=-0.5*sum(np.log(err**2 +s**2))
    #residual function from before just calculates the diff between model and obs
    term_3=-0.5*sum((res**2) /(err**2 +s**2))
        
        
    ln_L=term_1+(term_2+term_3)
    
    return -1*ln_L #NEGAITVE loglikeli---> MINIMIZE this
        



##########
#RVs#
##########


#v_array is an array  of the shifts done to each RV set

#for multiple RVs, these compute the likelihood
def L_RVs(t,data, err, s, v0, K, P, e, om, T0):
    """
    Calculates the loglikelihood of RV data.
        
    Parameters
    ----------
    t: array
        time of measurements
    
    data: array
        RV data
    v0:
        offset
    K,P,e,om,T0 : floats
        orbital model parameters
    s: float
        jitter
    
    Returns
    -------
    L : float
        The summed loglikelihood  *-1 -----> needs to be minimized
     """
    v_mod=RV_solo(v0,K,P,e,om,T0,t)
    RV_res=data-v_mod
    L=loglikelihood(RV_res,err,s)
    return L


def L_RVs_comb(t,data,err,s,v0,P,e,om,i,T0,a,parallax):
    """
    Calculates the loglikelihood of RV data when combined with astrometry
        
    Parameters
    ----------
    t: array
        time of measurements
    
    data: array
        RV data
    v0:
        offset
    P,e,om,i,T0,a,parallax : floats
        orbital model parameters
    s: float
        jitter
    
        
    Returns
    -------
    L : float
        The summed loglikelihood  *-1 -----> needs to be minimized
     """

    v_mod=RV_comb(v0,P,e,om,i,T0,a,parallax,t)
    RV_res=data-v_mod
    L=loglikelihood(RV_res,err,s)
    return L


###########
##HIP##
###########


#This function only works for gaia models for J2016 and Hipparcos data!
#hip_ad stands for hip "astrometric data" and hip_stand for hip "standard model solution"
def L_hipold(hip_ad,hip_stand,gaia,correction,P,e,om,i,Om,T0,a,s_hip=0):
        """
        Calculates the loglikelihood of a corrected gaia model with hipparcos data.
            

        Parameters
        ----------

        hip_ad : array
            astrometric hipparcos data
            A3,A4,A5,A6,A7,A8,A9=hip_ad
            
        hip_stand : array
            hipparcos standard model solution (contained in header of file mostly)
            
        gaia: array
            corrected gaia standard model solution
            one can basically put any standard model solution in here, that has standard epoch J2016.

        correction:
            changes made to the gaia standard model, these are 5 fit parameters

        P,e,om,i,Om,T0,a : floats
            orbital model parameters


        s_hip: float
            jitter for hipparcos data
            
        Returns
        -------
        L_hip : float
            The summed loglikelihood  *-1 -----> needs to be maximized
         """
        




    
        
        




        #Hipparchos
        #order data
        A3,A4,A5,A6,A7,A8,A9=hip_ad
        #hip_ad=np.array([A3,A4,A5,A6,A7])
        hip_stand=np.array(hip_stand)
        
        
        #calculating JD for every hip data
        t_HIP=hip_JD(hip_ad)
        
        #We are clculating corrections to the gaia model, but
        #Since gaia has the values for standardepoch J2016, we need to recompute asc_star and dec for year
        #J1991 to calculate the new residuals

        #gaia standard model solution

        c_gaia=stand_correct(gaia,correction)
        asc,dec,parallax,mu_a_star,mu_d=c_gaia
        
        #standardepochs for hip (1991) and gaia (2016) in JD
        t_1991=2448349.0625
        t_2016=2457389.0
        
        #pos recalc computes the new asc and dec position, we shift from 2016 to 1991
        gaia1991_asc,gaia1991_dec=pos_recalc(c_gaia,t_2016,t_1991)
                                        
        
        #The derivations in the hip data compute the change for the abscissa with given !!asc_star!! 
        #we need to multiply the catalogue values with the cos of their respective declination.
        gaia1991_asc_star=gaia1991_asc*np.cos(np.radians(gaia1991_dec))
        hip_stand[0]=hip_stand[0]*np.cos(np.radians(hip_stand[1]))

        #gaia with shifted asc_star and shifted dec to 1991
        gaia1991=(gaia1991_asc_star,gaia1991_dec,parallax,mu_a_star,mu_d)
        
        
        
                                         
                                         
                                         
        #the new hip residuals for gaia catalogue standard solution:

        #A8 is the abs residual
        #A8=np.array(A8)
        #new residual due to gaia correction:
        c_res_hip=abs_res(A8,gaia1991,hip_stand,hip_ad)
        

        
        
 
        #Now we have the remaining residuals, where the orbital motion is still contained
        #now we need to subtract the orbit, but in hipparcos manner
        x_O,y_O=orbit(P,e,om,i,Om,T0,a,t_HIP)
        
        res_hip_final=c_res_hip-(A3*x_O+A4*y_O)
        #we multiply the orbit positions with the respective hipparcos derivation and subtract them from
        #the remaining residual
        
        
        
       
        #error of hip residual
        A9=np.array(A9)

        
        L_hip=loglikelihood(res_hip_final,A9,s_hip)

        return L_hip


#This function only works for standard models with given standard epoch and Hipparcos data!
#hip_ad stands for hip "astrometric data" and hip_stand for hip "standard model solution"
def L_hip(hip_ad,hip_stand,standard_model,correction,P,e,om,i,Om,T0,a,Sepoch=2457389.0,s_hip=0):
        """
        Calculates the loglikelihood of a corrected gaia model with hipparcos data.
            

        Parameters
        ----------

        hip_ad : array
            astrometric hipparcos data
            A3,A4,A5,A6,A7,A8,A9=hip_ad
            
        hip_stand : array
            hipparcos standard model solution (contained in header of file mostly)
            
        standard_model: array
            standard model solution
            one can basically put any standard model solution in here, with given standard epoch "Sepoch"

        correction:
            changes made to the gaia standard model, these are 5 fit parameters

        P,e,om,i,Om,T0,a : floats
            orbital model parameters

        Sepoch: float
            standard epoch for standard model solution


        s_hip: float
            jitter for hipparcos data
            
        Returns
        -------
        L_hip : float
            The summed loglikelihood  *-1 -----> needs to be minimized
         """
        




    
        
        




        #Hipparchos
        #order data
        A3,A4,A5,A6,A7,A8,A9=hip_ad

        hip_stand=np.array(hip_stand)

        #calculating JD for every hip data

        t_HIP=hip_JD(hip_ad)

        
        
        
       
        
        #We are clculating corrections to the gaia model, but
        #Since gaia has the values for standardepoch J2016, we need to recompute asc_star and dec for year
        #J1991 to calculate the new residuals

        #gaia standard model solution

        c_stand=stand_correct(standard_model,correction)
        
        asc,dec,parallax,mu_a_star,mu_d=c_stand
        
        #standardepochs for hip (1991) and gaia DR3 (2016) in JD
        t_1991=2448349.0625
        
       
        
        #pos recalc computes the new asc and dec position, we shift from 2016 to 1991
        stand1991_asc,stand1991_dec=pos_recalc(c_stand,Sepoch,t_1991)
        #print(stand1991_asc)
                                        
        
        #The derivations in the hip data compute the change for the abscissa with given !!asc_star!! 
        #we need to multiply the catalogue values with the cos of their respective declination.
        stand1991_asc_star=stand1991_asc*np.cos(np.radians(stand1991_dec))
        hip_stand[0]=hip_stand[0]*np.cos(np.radians(hip_stand[1]))

        
        
        
        #gaia with shifted asc_star and shifted dec to 1991
        stand1991=np.array((stand1991_asc_star,stand1991_dec,parallax,mu_a_star,mu_d))


        #print(stand1991)

        
        
        
        
        
                                         
                                         
                                         
        

        #A8 is the abs residual
        #A8=np.array(A8)
        #new residual due to standard model correction:
        
        c_res_hip=abs_res(A8,stand1991,hip_stand,hip_ad)

        #print(c_res_hip)
    
    
        

        
        
 
        #Now we have the remaining residuals, where the orbital motion is still contained
        #now we need to subtract the orbit, but in hipparcos manner
        x_O,y_O=orbit(P,e,om,i,Om,T0,a,t_HIP)
        
        res_hip_final=c_res_hip-(A3*x_O+A4*y_O)
        #we multiply the orbit positions with the respective hipparcos derivation and subtract them from
        #the remaining residual
        
        
        
        
        
        #error of hip residual
        A9=np.array(A9)

        #print(res_hip_final)
        
        L_hip=loglikelihood(res_hip_final,A9,s_hip)

        return L_hip




###########
##HST##
###########


#def L_relative(hst_ad,hst_offsets,hst_proplax,c_gaia,P,e,om,i,Om,T0,a,earth,s_x_hst=0,s_y_hst=0):
#        
#
#        """
#        Calculates the loglikelihood of a corrected gaia model with relative astrometric data.
#            
#
#        Parameters
#        ----------
#        hst_ad : array
#            hipparcos standard model solution (contained in header of file mostly)
#            
#        hst_offsets : array,
#            offsets
#
#        hst_propxlax : array,
#            propermotion in RA,DEC and parallax, HST solution
#
#        c_gaia: array
#            corrected gaia standard model solution
#            one can basically put any standard model solution in here, that has standard epoch J2016.
#
#        P,e,om,i,Om,T0,a : floats
#            orbital model parameters
#
#        earth: 3D-array
#            position of earth for every time stamp of 2D measurements,
#            calculated preemtivly with earth_position(t)!
#
#
#        s_x_hst: float
#            jitter for hst data in RA
#
#        s_y_hst : float
#            hitter for hst data in dec
#            
#        Returns
#        -------
#        L_hst : float
#            The summed loglikelihood  *-1 -----> needs to be MAXIMIZED
#         """
#        #HST offset:
#        
#        #measured data:
#        t_HST,x_HST,x_err_HST,y_HST,y_err_HST=hst_ad
#
#        HST_asc,HST_dec=hst_offsets
#        
#        
#        #benedicts relative model
#        mu_a_star,mu_d,parallax=hst_proplax
#        
#        
#        
#        #gaia, benedict difference
#    
#        
#        
#        #calculate changes in tangential plane, need to add D_asc and D_dec at the end
#        
#        D_x0,D_y0=standard_model(c_gaia[0]+((HST_asc)/3.6e6)
#                                 ,c_gaia[1]+((HST_dec)/3.6e6)
#                                 ,hst_proplax[0]-c_gaia[2]
#                                 ,hst_proplax[1]-c_gaia[3]
#                                 ,hst_proplax[2]-c_gaia[4]
#                                 ,t_HST,earth,tangential=1)
#    
#        
#        #since we correct asc, but are right now in the tangential plane, cos(dec) is needed.
#        D_x=D_x0+HST_asc*np.cos(np.radians(c_gaia[1]+((HST_dec)/3.6e6)))
#        D_y=D_y0+HST_dec
#                         
#                                        
#        #
#        #calculate new residuals:
#        c_res_x=x_HST+D_x
#        c_res_y=y_HST+D_y
#        
#
#        #calculate residuals after we subtract the orbit:
#        xO,yO=orbit(P,e,om,i,Om,T0,a,t_HST)
#        
#        res_xO=c_res_x-xO
#        res_yO=c_res_y-yO
#        
#        
#                      
#                      
#    
#        L_hst_x=loglikelihood(res_xO,x_err_HST,s_x_hst)
#        L_hst_y=loglikelihood(res_yO,y_err_HST,s_y_hst)
#     
#     
#        return L_hst_x+L_hst_y #return is negative loglikeli ==> max this



def L_combined(RV_data,RV_err,t_RVs,#RV data
               hip_ad,hip_stand, #HIP data
               gaia,
               #Parameters:
               v_array,#RV offsets
               correction,P,e,om,i,Om,T0,a,#gaia solution and orbital parameters
               s_RVs=0,s_hip=0):
        """
        Calculates the loglikelihood of a corrected gaia model with hipparcos and RV data.
            

        Parameters
        ----------

        RV_data : list of n arrays 
        RV data [km/s]
        n being the amount of different RV datasets

        RV_err : list of n arrays
        RV error [km/s]
        n being the amount of different RV datasets

        t_RVs: list of n arrays
        RV times of measurement [JD]
        n being the amount of different RV datasets


        hip_ad : array
            astrometric hipparcos data
            
        hip_stand : array
            hipparcos standard model solution (contained in header of file mostly)
            
        gaia: array
            corrected gaia standard model solution
            one can basically put any standard model solution in here, that has standard epoch J2016.

        v_array: array
        Each offset for each RV dataset in ONE array

        correction: array
            changes made to the gaia standard model, these are 5 fit parameters

        P,e,om,i,Om,T0,a : floats
            orbital model parameters

        s_RVs: array
        jitter for each dataset

        s_hip: float
            jitter for hipparcos data
            
        Returns
        -------
        L_hip : float
            The summed loglikelihood  *-1 -----> needs to be maximized
         """
        

        
        parallax=gaia[2]+correction[2]
        
     
        Lrv=L_RVs_comb(RV_data,RV_err,s_RVs,v_array,P,e,om,i,T0,a,parallax,t_RVs)

        
        Lhip=L_hip(hip_ad,hip_stand,gaia,correction,P,e,om,i,Om,T0,a,s_hip)
        
        final=Lrv+Lhip
        

        return final


     
