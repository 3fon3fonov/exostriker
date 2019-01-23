c*************************************************************************
c                        RMVS.INC
c*************************************************************************
C Include file for the encounter subroutines
c
c Author:  Hal Levison  
c Date:    3/7/93
c Last revision: 

c...	scale factor for hill's sphere to take shorter time step
       real*8 RHSCALE
       parameter (RHSCALE=3.5)
c       parameter (RHSCALE=5.0)

c...	scale factor for hill's sphere to go to planet centric coord
       real*8 RHPSCALE
       parameter (RHPSCALE=1.0)
c       parameter (RHPSCALE=1.5)

c..    ratio of the number of time steps in the 
c...        encounter (helocentric) vs normal
       integer NTENC
c       parameter (NTENC=100)
       parameter (NTENC=10)

c..    ratio of the number of time steps in the 
c...        planetcentric encounter  vs heliocentric 
       integer NTPHENC
c       parameter (NTPHENC=5)
       parameter (NTPHENC=3)

c..    ratio of the number of time steps in the 
c...        encounter (planetcentric) vs normal
       integer NTPENC
       parameter (NTPENC=NTENC*NTPHENC)
c-------------------------------------------------------------------------
