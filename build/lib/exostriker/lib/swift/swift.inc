!*************************************************************************
!                        SWIFT.INC
!*************************************************************************
! Include file for SWIFT
!
! Author:  Hal Levison  
! Date:    2/2/93
! Last revision: 3/7/93

    implicit NONE    ! you got it baby

!...   Version of Swift
    real(8) VER_NUM
    parameter(VER_NUM=2.0d0)

!...   Maximum array size
    integer  NPLMAX, NTPMAX
!       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun 
    parameter  (NPLMAX = 51)   ! max number of planets, including the Sun
    parameter  (NTPMAX = 4096) ! max number of test particles

!...   Size of the test particle integer status flag
    integer NSTATP            ! Number of status parameters
    parameter  (NSTATP = 3)
    integer NSTAT            ! Number of status parameters
    parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet

!...   Size of the test particle integer status flag
    integer NSTATR
    parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

!...   convergence criteria for danby
    real(8) DANBYAC , DANBYB
    parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

!...    loop limits in the Laguerre attempts
    integer NLAG1, NLAG2
    parameter(NLAG1 = 50, NLAG2 = 400)

!...    A small number
    real(8) TINY
    PARAMETER(TINY=4.D-15)

!...    trig stuff
    real(8) PI,TWOPI,PIBY2,DEGRAD
    parameter (PI = 3.14159265358979D0)
    parameter (TWOPI = 2.0D0 * PI)
    parameter (PIBY2 = PI/2.0D0)
    parameter (DEGRAD = 180.0D0 / PI)

!-------------------------------------------------------------------------
