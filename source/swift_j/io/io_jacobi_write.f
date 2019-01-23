c*************************************************************************
c                          IO_JACOBI_WRITE.F
c*************************************************************************
c Does the write for anal_jacobi_write
c
c      Input:
c            i1st           ==>  =0 if first write, =1 if not (int scalar)
c            t              ==>  current time (real scalar)
c            jac0           ==>  Initial values of the jacobio constants 
c                                  (real array)
c            dj             ==>  change in Jacobi const. for part (real arrays)
c            nw             ==>  number of bodies (int scalar)
c            iu             ==>  unit to write to
c            fopenstat      ==>  The status flag for the open 
c                                statements of the output files.  
c                                          (character*80)
c
c Remarks: If the particle is not active a value of -100 is written out 
c Authors:  Hal Levison 
c Date:    8/12/93
c Last revision: 10/3/96

      subroutine io_jacobi_write(i1st,t,jac0,dj,nw,iu,fopenstat)


      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer nw,iu,i1st
      real*8 t,dj(NTPMAX),jac0(NTPMAX)
      character*(*) fopenstat

c...  Internals
      integer i,ierr

c----
c...  Executable code 

      if(i1st.eq.0) then

         call io_open(iu,'jacobi.out',fopenstat,'FORMATTED',ierr)

         if(ierr.ne.0) then
           write(*,*) ' SWIFT ERROR: in anal_jacobi_write '
           write(*,*) '     Could not open jacobi.out '
           call util_exit(1)
         endif

         write(*,*) ' Initial values of the Jacobi Constant: '
         write(*,2) (jac0(i),i=1,nw)

         write(iu,2) t,(dj(i),i=1,nw)
 2       format(1x,11(2x,1p1e12.5))

      else

         call io_open(iu,'jacobi.out','append','FORMATTED',ierr)
         write(iu,2) t,(dj(i),i=1,nw)
         
      endif

      close(iu)

      return
      end                       ! io_jacobi_write
c-------------------------------------------------------------------------
