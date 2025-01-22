
import numpy as np
from numba import jit

from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric


#calculating thiele constants
#Output in mas
def thiele(a,omega,Omega,i):
    """Thiele Innes Constants"""
    
    
    #omega=omega+np.pi
    
    A=a*(np.cos(omega)*np.cos(Omega)-np.sin(omega)*np.sin(Omega)*np.cos(i)) 
    B=a*(np.cos(omega)*np.sin(Omega)+np.sin(omega)*np.cos(Omega)*np.cos(i)) 
    F=a*(-np.sin(omega)*np.cos(Omega)-np.cos(omega)*np.sin(Omega)*np.cos(i)) 
    G=a*(-np.sin(omega)*np.sin(Omega)+np.cos(omega)*np.cos(Omega)*np.cos(i)) 
    return A,B,F,G


@jit(nopython=True)
def calc_E(e,M,n=30):
    """Newton Raphson Algorythm to calculate the eccentric anomaly given the mean anomaly and the eccentricity"""
    
   
    final=np.ones(len(M))
    for index, j in enumerate(M):
        if e<=0.8:
            E=j
        if e>0.8:
            E=np.pi
        values=np.ones(n)
        values[0]=E
        i=1
        while i<n:
            values[i]=values[i-1]-(values[i-1]-e*np.sin(values[i-1])-j)/(1-e*np.cos(values[i-1]))
            if abs(values[i]-values[i-1])<1e-8:
                values[-1]=values[i]
                break
            i=i+1
            if i==n:
                return print("error in calculating E")
        final[index] = values[-1]
    
    # final=np.array(final)
            
    return final

def orbit(P,e,om,i,Om,T0,a,t):
    
    
    #calculating thiele constants
    const=thiele(a,om,Om,i)
    
    #calculating E
    M=2*np.pi*((t-T0)%P)/P
    E=calc_E(e,M)

    #calculating elliptical rectangular coords from E and e
    X=np.cos(E)-e
    Y=(1-e**2)**0.5 *np.sin(E)
    
    #const=A,B,F,G
    #with thiele, X,Y, we can calculate the final position x,y for a time t
    x=const[1]*X+const[3]*Y
    y=const[0]*X+const[2]*Y
    
    
    
    
    #to provide a right handed system and fulfill conventions, x is NORTH and y is EAST
    
    return x,y



###Standard model

#position of earth, important for parallax factors
def earth_position(t): 
    """Computes the earths barycentric position in AU given a time in JD"""
    t=Time(t,format="jd")
    with solar_system_ephemeris.set("de440"):
        loc=get_body_barycentric("earth",t) #earth position from barycentre
        
        loc=loc.xyz.value/1.4959787e+8 #convert from km to AU
    #the aequinox is at -1 0 0, north summer is at 0 -1 -0.4!
    return loc


def parallax_factors(asc,dec,pos_earth): #earth pos x,y,z, in AU #asc,dec in degree
    """Calculates the parallax factors depending on earths position and the objects position at the sky"""
    
    asc=np.radians(asc)
    dec=np.radians(dec)
    
    p_a=(-pos_earth[0]*np.sin(asc)+pos_earth[1]*np.cos(asc))/np.cos(dec)
    
    p_d=(-pos_earth[0]*np.cos(asc)-pos_earth[1]*np.sin(asc))*np.sin(dec)+pos_earth[2]*np.cos(dec)
    return p_a,p_d #parallax factors



def standard_model(asc,dec,parallax,mu_a_star,mu_d,t,earth,Sepoch=2457389.0,tangential=True):

    """
    Calculates the position of an object in respect to a standard epoch

    Parameters:
    ---------
    asc,dec,mu_a_star,mu_d : floats
        the standard model solution for Sepoch
    t: array
        timestamps

    earth: array
        position of earth at given timestamps,(use ExTRA.earth_position(t) to compute these)

    Sepoch: float
        standard epoch of input standard model

    tangential: bool
        gives tangential position if True and absolute position if false, both in [mas]
    
    Returns:
    ----------
    asc_final,dec_final : Tuple,floats
        new coordinates for Epoch1 in [mas]
    """

    
    
    #give asc and deg in degree
    
    #tangential=False means total position
    #tangential=True leaves out asc and dec and returns the position in the tangential plane
    
    t0=Sepoch
    #t0=time of asc and dec measurement, so a standard epoch, its mostly J2000, or 2451545.0JD, but
    #for cases like Hipparchos its J1991.25, or 2448349.0625JD
    #for gaia it is J2016, or 2456389.0
    
    
    asc_star=asc*np.cos(np.radians(dec)) #we need the RA* to calculate the final position
    
    p_a,p_d=parallax_factors(asc,dec,earth)
    if tangential:
        
        
        a=mu_a_star*((t-t0)/365.25)+p_a*parallax
        d=mu_d*((t-t0)/365.25)+p_d*parallax
    
    if not tangential:
        asc_star=asc_star*(3.6e6) #convert to mas
        dec=dec*(3.6e6) #convert to mas

    #t0 is the time for given asc and dec.
    #plugging in t lets us calculate p_a and p_d and therefore the new star position.


        a=asc_star+mu_a_star*((t-t0)/365.25)+p_a*parallax
        d=dec+mu_d*((t-t0)/365.25)+p_d*parallax
    
    return a,d #in mas

#this function recalculates asc and dec between STANDARD EPOCHS meaning for example from J2000 to J2016.
#epoch0 is the old and epoch1 the new epoch
def pos_recalc(standmodel,Epoch0,Epoch1):
    """
    Recalculates the position of an object given a different standard epoch

    Parameters:
    ---------
    standmodel: array
        the standard model solution for Epoch0
        given as asc,dec,parallax,mu_a_star,mu_d
    Epoch0: float [JD]
        the standard epoch for given solution, in [JD]
    Epoch1: float [JD]
        the new standard epoch for given solution
    
    Returns:
    ----------
    asc_final,dec_final : Tuple,floats
        new coordinates for Epoch1 in [mas]
    """
    asc,dec,parallax,mu_a_star,mu_d=standmodel
    asc_star=asc*np.cos(np.radians(dec))
    
    asc_star=asc_star*(3.6e6) #convert to mas
    dec=dec*(3.6e6) #convert to mas


    a=asc_star+mu_a_star*((Epoch1-Epoch0)/365.25)
    d=dec+mu_d*((Epoch1-Epoch0)/365.25)
        
    
    dec_final=d/(3.6e6)
    asc_final=a/(3.6e6*np.cos(np.radians(dec_final)))
        
    return asc_final,dec_final








def secondary_mass(M_prime,parallax,P,e,i,a,unit="jup"):
    """computes Mass of the secondary"""



    M_jup=1.89813*1e27
    M_sun=1.989*1e30
    
    K1=2*np.pi*a*np.sin(i)*1.496e11/(parallax)
    K2=(P*86400*(1-e**2)**0.5)
    K=K1/K2
    G=6.6743*1e-11
    print("K:")
    print(K)
    #print((2*np.pi*G)*(M_s**2))
    M_s=((K**3) *(P*24*60**2)/(2*np.pi*G*np.sin(i)**3)*(M_prime**2))**(1/3)

    if unit=="jup":
        return M_s/M_jup
    if unit=="sun":
        return M_s/M_jup


def stand_correct(stand,correction):
    """Corrects a standard model with a correction by adding the correction onto it.
        Parameters:
    ---------
    stand: array
    the standard model (asc,dec,parallax,mu_a,mu_d) in [deg,deg,mas,mas/y,mas/y]

    correction: array
    the correction, same format as standard model in [mas!,mas!,mas,mas/y,mas/yr]

    
    Returns: array
    corrected standard solution in [deg,deg,mas,mas/y,mas/y]
    ----------
    
    """

    new=np.zeros(5)


    

    
    new[1]=stand[1]+correction[1]/(3.6e6)
    new[2]=stand[2]+correction[2]
    new[3]=stand[3]+correction[3]
    new[4]=stand[4]+correction[4]

    stand_star=stand[0]*np.cos(np.radians(new[1]))

    stand_star_shifted=stand_star+correction[0]/3.6e6

    new[0]=stand_star_shifted/(np.cos(np.radians(new[1])))


    

    return new
