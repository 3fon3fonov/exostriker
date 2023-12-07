c*************************************************************************
c                        BS.INC
c*************************************************************************
C Include file for the bs integrator
c
c Author:  Hal Levison  
c Date:    5/17/93
c Last revision: 

c...   Size of internal arrays
       integer N6DBS
       parameter(N6DBS=(6*(NTPMAX+NPLMAX)))
       integer NTEMP
       parameter(NTEMP=(72*(NTPMAX+NPLMAX)))

c...   maximum numbers of internations
       integer NTRYS
       parameter(NTRYS=30)
  
c------------------------------------------------------------------------