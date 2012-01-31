      subroutine advance(f,neqn,y,t,tout,relerr,abserr,iflag)
Cf2py intent(in,hide) neqn
Cf2py intent(in) t,tout,relerr,abserr 
Cf2py intente(out) iflag
Cf2py intent(in,out) y
      integer neqn,iflag
      double precision y(neqn),t,tout,relerr,abserr
      external f

      integer iwork(5)
      double precision work(3+6*neqn)
      logical finished 
      integer triedcounter
      character*80 msg

      triedcounter = 0
      finished = .false.
      iflag = 1
      msg = ''
      do while (.not. finished)
         call rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
         triedcounter = triedcounter + 1
c     succeed
         if (iflag .lt. 3) go to 300
c     Tried more than two times. Return with failure.
         if (triedcounter .gt. 2) go to 300
c     Tried again when iflag = 3 or 4.
         if (iflag .lt. 5) go to 500
c     Invalid inputs
         if (iflag .eq. 8) then
            msg = 'Invalid input parameters'
            go to 300
         end if
c     One-step integrator is altered as suggestion.
         if (iflag .eq. 7) then
            msg = 'One-step integrator is altered as suggestion'
            iflag = -1
            go to 400
         end if
c     Too much accuracy requirement.
         if (iflag .eq. 6) then
            msg = 'Too much accuracy requirement.'
            go to 300
         end if
c     atol should be positive to continue.
         msg = 'atol should be positive to continue.'     

 300     finished = .true.
 400     if (iflag .gt. 2) write (*,*) msg
 500  continue
      end do
      return
      end
