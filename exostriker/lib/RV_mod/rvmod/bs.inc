!*************************************************************************
!                        BS.INC
!*************************************************************************
! Include file for the bs integrator
!
! Author:  Hal Levison  
! Date:    5/17/93
! Last revision: 

!...   Size of internal arrays
integer N6DBS
parameter(N6DBS = (6 * (NTPMAX + NPLMAX)))
integer NTEMP
parameter(NTEMP = (72 * (NTPMAX + NPLMAX)))

!...   maximum numbers of internations
integer NTRYS
parameter(NTRYS = 300)

!------------------------------------------------------------------------