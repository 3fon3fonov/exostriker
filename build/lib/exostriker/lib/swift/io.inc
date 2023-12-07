!*************************************************************************
!                        IO.INC
!*************************************************************************
! Include file for io routines
!
! Author:  Hal Levison  
! Date:    2/21/94
! Last revision: 2/24/94

!...   number of bytes in iflgchk
    integer IO_NBITS
    parameter(IO_NBITS=6)
!                          bit 0 set ==>  write big binary data file
!                          bit 1 set ==>  write real*4 binary file rather 
!                                         than int*2: ignored if bit0=F
!                          bit 2 set ==>  calc energy of system wrt time
!                          bit 3 set ==>  calc jacobi of the test particles
!                          bit 4 set ==>  check if particles are removed
!                          bit 5 set ==>  include J2 and J4 terms
!-------------------------------------------------------------------------
