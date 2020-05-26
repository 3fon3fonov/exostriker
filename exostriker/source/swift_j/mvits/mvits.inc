c*************************************************************************
c                        MVITS.INC
c*************************************************************************
C Include file for the individual time step code
c
c Author:  Hal Levison  
c Date:    12/16/93
c Last revision: 

c..   Maximum  number of different timesteps
       integer NUMDT
       parameter (NUMDT=5)


c...	Maximum ratio of time steps
       integer NDTMAX
       parameter (NDTMAX=2**NUMDT)     ! 32

c-------------------------------------------------------------------------
