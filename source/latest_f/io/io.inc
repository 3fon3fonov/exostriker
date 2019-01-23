c*************************************************************************
c                        IO.INC
c*************************************************************************
C Include file for io routines
c
c Author:  Hal Levison  
c Date:    2/21/94
c Last revision: 2/24/94

c...   number of bytes in iflgchk
       integer IO_NBITS
       parameter(IO_NBITS=6)
c                          bit 0 set ==>  write big binary data file
c                          bit 1 set ==>  write real*4 binary file rather 
c                                         than int*2: ignored if bit0=F
c                          bit 2 set ==>  calc energy of system wrt time
c                          bit 3 set ==>  calc jacobi of the test particles
c                          bit 4 set ==>  check if particles are removed
c                          bit 5 set ==>  include J2 and J4 terms
c-------------------------------------------------------------------------
