c*************************************************************************
c                            RMVS_STEP_IN_MVPL.F
c*************************************************************************
c Rearrange the order of the planets for planocentric encounter. Also
c moves the arraw from xpltb to xplte
c
c             Input:
c                 ipl            ==>  The planet in the center
c                 nbod           ==>  number of massive bodies (int scalar)
c                 mass           ==>  Heliocentric mass of bodies (real array)
c                 xpl,ypl,zpl    ==>  Heliocentric position of planet wrt time
c                                       (real arrays)
c             Output:
c                 masst          ==>  rearranged masses
c              xpltb,ypltb,zpltb ==>  rearranged pl x's at beginning of step
c                                       (1d real arrays)
c              xplte,yplte,zplte ==>  rearranged pl x's at end of step
c                                       (1d real arrays)
c
c
c Remarks: Adopted from hal's wiscl_fk.f
c Authors:  Hal Levison 
c Date:    2/19/93
c Last revision: 

      subroutine rmvs_step_in_mvpl(ipl,nbod,mass,xpl,ypl,zpl,
     &             masst,xpltb,ypltb,zpltb,xplte,yplte,zplte)


      include '../swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer ipl,nbod
      real*8 mass(NPLMAX)
      real*8 xpl(NPLMAX),ypl(NPLMAX)
      real*8 zpl(NPLMAX)

c...  Outputs:
      real*8 xpltb(NPLMAX),ypltb(NPLMAX),zpltb(NPLMAX)
      real*8 xplte(NPLMAX),yplte(NPLMAX),zplte(NPLMAX)
      real*8 masst(NPLMAX)

c...  Internals
      integer i

c----
c...  Executable code 

c...  Move things over
      do i = 1,nbod
         xpltb(i) = xplte(i) 
         ypltb(i) = yplte(i) 
         zpltb(i) = zplte(i) 
      enddo

c...  first just move the planets
      do i = 2,nbod
         masst(i) = mass(i)
         xplte(i) = xpl(i) - xpl(ipl)
         yplte(i) = ypl(i) - ypl(ipl)
         zplte(i) = zpl(i) - zpl(ipl)
      enddo

c...  now switch things
      masst(1) = mass(ipl)
      xplte(1) = 0.0
      yplte(1) = 0.0
      zplte(1) = 0.0

      masst(ipl) = mass(1)
      xplte(ipl) = -xpl(ipl)
      yplte(ipl) = -ypl(ipl)
      zplte(ipl) = -zpl(ipl)

      return
      end   ! rmvs_step_in_mvpl
c---------------------------------------------------------------------------


