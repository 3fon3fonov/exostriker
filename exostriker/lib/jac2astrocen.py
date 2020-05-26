#!/usr/bin/python
 


#from array import *
#from numpy import *
#from pylab import *
 
import numpy as np

G  = 6.67384e-11 
pi = np.pi

twopi = 2.0 * pi
piby2 = .50 * pi
AU_m = 149597870700.00000



 
##################################################################################################################

 

#**********************************************************************
#                    ORBEL_ZGET.F
#**********************************************************************
#     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
#*          given Q (Fitz. notation.)
#*
#*             Input:
#*                           q ==>  parabola mean anomaly. (real scalar)
#*             Returns:
#*                  orbel_zget ==>  eccentric anomaly. (real scalar)
#*
#*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
#*     REMARKS: For a parabola we can solve analytically.
#*     AUTHOR: M. Duncan 
#*     DATE WRITTEN: May 11, 1992.
#*     REVISIONS: May 27 - corrected it for negative Q and use power
#*              series for small Q.
#***********************************************************************

def orbel_zget(q):

    iflag = 0

    if q < 0.0:
        iflag = 1.
        q = -q
 

    if  q < 1.e-3:
        orbel_zget = q*(1.0 - (q*q/3.0)*(1.0 -q*q))
    else:
        x = 0.50*(3.0*q + np.sqrt(9.0*(q**2) +4.0))
        tmp = x**(1.0/3.0)
        orbel_zget = tmp - 1.0/tmp
 

    if iflag == 1: 
        orbel_zget = -orbel_zget
        q = -q
 
    return orbel_zget
 
# Not tested!
#------------------------------------------------------------------------------
#

#------------------------------------------------------------------
#
#**********************************************************************
#                    ORBEL_FHYBRID.F
#**********************************************************************
#     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
#
#             Input:
#                           e ==> eccentricity anomaly. (real scalar)
#                           n ==> hyperbola mean anomaly. (real scalar)
#             Returns:
#               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
#
#     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
#                 For larger N, uses FGET
#     REMARKS: 
#     AUTHOR: M. Duncan 
#     DATE WRITTEN: May 26,1992.
#     REVISIONS: 
#     REVISIONS: 2/26/93 hfl
#**********************************************************************

def orbel_fhybrid(e,n):

    abn = n
    if n > 0.0:
        abn = -abn

    if abn < 0.6360*e -0.60:
        orbel_fhybrid = orbel_flon(e,n)
    else:
        orbel_fhybrid = orbel_fget(e,n)

    return orbel_fhybrid

#not implemented yet
#**********************************************************************

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      MCO_KEP.FOR    (ErikSoft  7 July 1999)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Author: John E. Chambers 
# translated to Py: T.T. 
#c
# Solves Kepler's equation for eccentricities less than one.
# Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
#
#  e = eccentricity
#  l = mean anomaly      (radians)
#  u = eccentric anomaly (   "   )
#
#------------------------------------------------------------------------------
#
      
def mco_kep(e,oldl):
 

#
# Reduce mean anomaly to lie in the range 0 < l < pi
    if oldl >= 0:
        l = oldl%twopi

    else:
        l = (oldl%twopi) + twopi

    sign = 1.0

    if l > pi:
        l = twopi - l
        sign = -1.0
    ome = 1.0 - e

    if l >= 0.450 or e < 0.550:

    #
    # Regions A,B or C in Nijenhuis
    # -----------------------------
    #
    # Rough starting value for eccentric anomaly
        if  l > ome:
            u1 = ome
        else:
            if  l > (pi-1.0-e):
                u1 = (l+e*pi)/(1.0+e)
            elif l <= (pi-1.0-e):
                u1 = l + e
 
    #
    # Improved value using Halley's method
 
        if u1 > piby2:
            x = pi - u1
        else:
            x = u1
 
        x2 = x*x
        sn = x*(1.0 + x2*(-0.16605 + x2*0.00761) )
        dsn = 1.0 + x2*(-0.49815 + x2*0.03805)

        if u1 > piby2: 
            dsn = -dsn

        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1.0 - e*dsn
        u2 = u1 - f0/(f1 - 0.50*f0*f2/f1)

    else:
    #
     #Region D in Nijenhuis
    # ---------------------
    #
    # Rough starting value for eccentric anomaly
        z1 = 4.0*e + 0.50
        p = ome / z1
        q = 0.50 * l / z1
        p2 = p*p
        z2 = np.exp(np.log(np.sqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.0*q / ( z2 + p + p2/z2 )
    #
    # Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - 0.0750*u1*z3 / (ome + z1*z2 + 0.3750*z3)
        u2 = l + e*u2*( 3.0 - 4.0*u2*u2 )
 
#
# Accurate value using 3rd-order version of Newton's method
# N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!
#
# First get accurate values for u2 - sin(u2) and 1 - cos(u2)
 
    if  u2 > piby2:
        z3 = pi - u2

    else:
        z3 = u2
 
    if z3 > (0.50*piby2):
        x = piby2 - z3
    else:
        x = z3
 
    x2 = x*x
    ss = 1.0
    cc = 1.0
    
    ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
    cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. - x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. - x2/306.))))))))
 
    if z3 > (0.50*piby2):
        z1 = cc + z3 - 1.0
        z2 = ss + z3 + 1.0 - piby2
    else:
        z1 = ss
        z2 = cc
 
    if u2 > piby2:
        z1 = 2.0*u2 + z1 - pi
        z2 = 2.0 - z2
 
    f0 = l - u2*ome - e*z1
    f1 = ome + e*z2
    f2 = .50*e*(u2-z1)
    f3 = e/6.0*(1.0-z2)
    z1 = f0/f1
    z2 = f0/(f2*z1+f1)
    temp = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
# 
#------------------------------------------------------------------------------
#

# works OK!

    return temp

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      MCO_SIN.FOR    (ErikSoft  12 June 1998)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Calculates sin and cos of an angle X (in radians)
#
#------------------------------------------------------------------------------
#

def mco_sine (x):
 
    if x > 0:
        x = x%twopi
    else:
        x =  x%twopi + twopi
 
    cx = np.cos(x)
 
    if  x > pi:
        sx = -np.sqrt(1.0 - cx*cx)
    else:
        sx =  np.sqrt(1.0 - cx*cx)
    return sx,cx

# # works OK!
#------------------------------------------------------------------------------
#



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      MCO_SINH.FOR    (ErikSoft  12 June 1998)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Calculates sinh and cosh of an angle X (in radians)
#
#------------------------------------------------------------------------------
#

def mco_sinh(x):
 
    sx = np.sinh(x)
    cx = np.sqrt(1.0 + sx*sx)
 
    return sx,cx
 
# # works OK!
#------------------------------------------------------------------------------
#
 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      MCO_EL2X.FOR    (ErikSoft  7 July 1999)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Author: John E. Chambers 
# translated to Py: T.T. 
#
# Calculates Cartesian coordinates and velocities given Keplerian orbital
# elements (for elliptical, parabolic or hyperbolic orbits).
#
# Based on a routine from Levison and Duncan's SWIFT integrator.
#
#  gm = grav const * (central + secondary mass)
#  q = perihelion distance
#  e = eccentricity
#  i = inclination                 )
#  p = longitude of perihelion !!! )   in
#  n = longitude of ascending node ) radians
#  l = mean anomaly                )
#
#  x,y,z = Cartesian positions  ( units the same as a )
#  u,v,w =     "     velocities ( units the same as sqrt(gm/a) )
#
#------------------------------------------------------------------------------

def mco_el2x(gm,q,e,i,p,n,l):
 
 
    g = p - n

# Rotation factors
    si,ci = mco_sine(i)
    sg,cg = mco_sine(g)
    sn,cn = mco_sine(n)
    z1 = cg * cn
    z2 = cg * sn
    z3 = sg * cn
    z4 = sg * sn
    d11 =  z1 - z4*ci
    d12 =  z2 + z3*ci
    d13 = sg * si
    d21 = -z3 - z2*ci
    d22 = -z4 + z1*ci
    d23 = cg * si

# Semi-major axis
    a = q / (1.0 - e)

# Ellipse
    if  e < 1.0:
        romes = np.sqrt(1.0 - e*e)
        temp = mco_kep(e,l)
        se,ce = mco_sine(temp)
        z1 = a * (ce - e)
        z2 = a * romes * se
        temp = np.sqrt(gm/a) / (1.0 - e*ce)
        z3 = -se * temp
        z4 = romes * ce * temp
    else:
# Parabola
        if e == 1.0:  
            ce = orbel_zget(l)
            z1 = q * (1.0 - ce*ce)
            z2 = 2.0 * q * ce
            z4 = np.sqrt(2.0*gm/q) / (1.0 + ce*ce)
            z3 = -ce * z4
        else:
# Hyperbola
            romes = np.sqrt(e*e - 1.0)
            temp = orbel_fhybrid(e,l) #not implemented?
            se,ce = mco_sinh(temp)
            z1 = a * (ce - e)
            z2 = -a * romes * se
            temp = np.sqrt(gm/abs(a)) / (e*ce - 1.0)
            z3 = -se * temp
            z4 = romes * ce * temp

    x = d11 * z1  +  d21 * z2
    y = d12 * z1  +  d22 * z2
    z = d13 * z1  +  d23 * z2
    u = d11 * z3  +  d21 * z4
    v = d12 * z3  +  d22 * z4
    w = d13 * z3  +  d23 * z4
 
    return x,y,z,u,v,w
 

 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      MCO_J2H.FOR    (ErikSoft   2 March 2001)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Author: John E. Chambers 
# translated to Py: T.T. 
#
# Converts Jacobi coordinates to coordinates with respect to the central
# body.
# N.B. The Jacobi coordinates of the small bodies are assumed to be equal
# ===  to their coordinates with respect to the central body.
#
#------------------------------------------------------------------------------
#

def mco_j2h(m,x,v):


    xh = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    vh = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]

    nbig = len(m) -1

    mtot = m[1]
    temp = m[1] / (mtot + m[0])
    mx = temp * x[0][1]
    my = temp * x[1][1]
    mz = temp * x[2][1]
    mu = temp * v[0][1]
    mv = temp * v[1][1]
    mw = temp * v[2][1]

    xh[0][1] = x[0][1] + mx
    xh[1][1] = x[1][1] + my
    xh[2][1] = x[2][1] + mz
    vh[0][1] = v[0][1] + mu
    vh[1][1] = v[1][1] + mv
    vh[2][1] = v[2][1] + mw


    for jj in range(2,3):
        j =  jj
 
        xh[0][j] = x[0][j] + mx
        xh[1][j] = x[1][j] + my
        xh[2][j] = x[2][j] + mz
        vh[0][j] = v[0][j] + mu
        vh[1][j] = v[1][j] + mv
        vh[2][j] = v[2][j] + mw
        mtot = mtot + m[j]
        temp = m[j] / (mtot + m[0])
        mx = mx  +  temp * x[0][j]
        my = my  +  temp * x[1][j]
        mz = mz  +  temp * x[2][j]
        mu = mu  +  temp * v[0][j]
        mv = mv  +  temp * v[1][j]
        mw = mw  +  temp * v[2][j]

    if  nbig > 2:
        xh[0][nbig] = x[0][nbig] + mx
        xh[1][nbig] = x[1][nbig] + my
        xh[2][nbig] = x[2][nbig] + mz
        vh[0][nbig] = v[0][nbig] + mu
        vh[1][nbig] = v[1][nbig] + mv
        vh[2][nbig] = v[2][nbig] + mw
 
 
    return xh,vh


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      MCO_H2B.FOR    (ErikSoft   2 March 2001)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Author: John E. Chambers
#
# Converts coordinates with respect to the central body to barycentric
# coordinates.
#
#------------------------------------------------------------------------------
#

def mco_h2b(m,xh,vh): 


 
# Calculate coordinates and velocities of the central body
        
    xb = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    vb = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    mtot2 = 0.0
 
#Calculate coordinates and velocities of the central body
    for jj in range(1,3):
        j = 3-jj
        mtot2 = mtot2  +  m[j]
        xb[0][0] = xb[0][0]  +  m[j] * xh[0][j]
        xb[1][0] = xb[1][0]  +  m[j] * xh[1][j]
        xb[2][0] = xb[2][0]  +  m[j] * xh[2][j]
        vb[0][0] = vb[0][0]  +  m[j] * vh[0][j]
        vb[1][0] = vb[1][0]  +  m[j] * vh[1][j]
        vb[2][0] = vb[2][0]  +  m[j] * vh[2][j]
 

    temp = -1.0 / (mtot2 + m[0])
    xb[0][0] = temp * xb[0][0]
    xb[1][0] = temp * xb[1][0]
    xb[2][0] = temp * xb[2][0]
    vb[0][0] = temp * vb[0][0]
    vb[1][0] = temp * vb[1][0]
    vb[2][0] = temp * vb[2][0]
 

#Calculate the barycentric coordinates and velocities
    for jj in range(1,3):
        j = 3-jj
        xb[0][j] = xh[0][j] + xb[0][0]
        xb[1][j] = xh[1][j] + xb[1][0]
        xb[2][j] = xh[2][j] + xb[2][0]
        vb[0][j] = vh[0][j] + vb[0][0]
        vb[1][j] = vh[1][j] + vb[1][0]
        vb[2][j] = vh[2][j] + vb[2][0]
 

    return xb,vb
#
#------------------------------------------------------------------------------
 
 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      MCO_X2EL.FOR    (ErikSoft  23 January 2001)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Author: John E. Chambers 
# translated to Py: T.T. 
#
# Calculates Keplerian orbital elements given relative coordinates and
# velocities, and GM = G times the sum of the masses.
#
# The elements are: q = perihelion distance
#                   e = eccentricity
#                   i = inclination
#                   p = longitude of perihelion (NOT argument of perihelion!!)
#                   n = longitude of ascending node
#                   l = mean anomaly (or mean longitude if e < 1.e-8)
#
 
def mco_x2el(gm,x,y,z,u,v,w):
 
    hx = y * w  -  z * v
    hy = z * u  -  x * w
    hz = x * v  -  y * u
    h2 = hx*hx + hy*hy + hz*hz
    v2 = u * u  +  v * v  +  w * w
    rv = x * u  +  y * v  +  z * w
    r = np.sqrt(x*x + y*y + z*z)
    h = np.sqrt(h2)
    s = h2 / gm
 
# Inclination and node
    ci = hz / h
    if abs(ci) < 1.0:
        i = np.acos(ci)
        n = np.atan2(-hy,hx) #(hx,-hy) 
        if  n < 0.0: 
            n = n + twopi
    else:
        if ci > 0.0:  
            i = 0.0
        if ci < 0.0: 
            i = pi
        n = 0.0
 
#
# Eccentricity and perihelion distance
    temp = 1.0  +  s * (v2 / gm  -  2.0 / r)
    if temp <= 0:
        e = 0.0
    else:
        e = np.sqrt(temp)
 
    q = s / (1.0 + e)
#
# True longitude
    if  hy != 0.0:
        to = -hx/hy
        temp = (1.0 - ci) * to
        tmp2 = to * to
        true = np.atan2((y*(1.0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp))
    else:
        true = np.atan2(y * ci, x)
 
    if ci < 0: 
        true = true + pi
 
    if e < 3.e-8:
        p = 0.0
        l = true
    else:
        ce = (v2*r - gm) / (e*gm)
#
# Mean anomaly for ellipse
        if e < 1:
            if abs(ce) > 1.0: 
                ce = np.copysign(1.0,ce)
            bige = np.acos(ce)
            if rv < 0: 
                bige = twopi - bige
            l = bige - e*np.sin(bige)
        else:
 
# Mean anomaly for hyperbola
            if ce < 1.0: 
                ce = 1.0
            bige = np.log( ce + np.sqrt(ce*ce-1.0) )
            if rv < 0.0: 
                bige = -bige
            l = e*np.sinh(bige) - bige
 
# Longitude of perihelion
        cf = (s - r) / (e*r)
        if  abs(cf) > 1.0: 
            cf = np.copysign(1.0,cf)
        f = np.acos(cf)
        if  rv < 0.0: 
            f = twopi - f
        p = true - f
        p = (p + twopi + twopi)%twopi
  
    if l < 0.0:
        l = l + twopi
    if l > twopi:
        l = l%twopi
 
    return q,e,i,p,n,l
 
 
 


