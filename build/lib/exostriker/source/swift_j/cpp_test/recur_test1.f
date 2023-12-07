      i = 1
      call count(i)
      stop
      end

c-------------------------------------------
      subroutine count(i)
      write(*,*) i
      i = i + 1
      if(i.le.10) then
         call count(i)
      endif
      return
      end
