      subroutine solve(spcrad,neq,f,yin,yout,trr,nt,rtol,atol,info,
     1     idid, termin, itermin, nstop)
c--------------------------------------------------------------------------
c
c  This is a wrapper for RKC.f to solve ODE on a sequence of time points.
c  Source code and details of RKC.f can be found in website of netlib.
c---------------------------------------------------------------------------
      external spcrad, f
      integer neq, nt, info, idid, nstop
      double precision yin, yout, trr, rtol, atol
      dimension yin(neq), yout(nt,neq), trr(nt), atol(*), info(4)
c------------------------------------------------------------------
c  RKCDID is a labelled common block that communicates statistics
c         about the integration process:
c         common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm
c
c         The integer counters are:
c
c        NFE      number of evaluations of F used
c                 to integrate the initial value problem
c        NSTEPS   number of integration steps
c        NACCPT   number of accepted steps
c        NREJCT   number of rejected steps
c        NFESIG   number of evaluations of F used 
c                 to estimate the spectral radius
c        MAXM     maximum number of stages used
c
c        This data can be used to monitor the work and terminate a run
c        that proves to be unacceptably expensive.  Also, if MAXM should
c        be far beyond 100, the problem is too expensive for RKC and
c        alternative methods should be considered.
      integer           nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm
c------------------------------------------------------------------
      double precision t, tend, ycurrent, spcrad
      double precision, allocatable ::  work(:)
      logical finish, repeated
      integer n, i, resterm, ptr
      dimension ycurrent(neq)
c------Initial status----------------------------------------------

      if (info(2) .eq. 0) then
         allocate(work(8 + 5*neq))
      else
         allocate(work(8 + 4*neq))
      end if
      finish = .false.
      repeated = .false.
      n = 1
      resterm = 0
      do 10 i = 1, neq
 10      yout(1, i) = yin(i)
      do while ((n .lt. nt) .and. (.not. finish))
         t = trr(n)
         tend = trr(n + 1)
         idid = 0
         call rkc(spcrad,neq,f,yin,t,tend,rtol,atol,info,
     &          work,idid)
         if (maxm .ge. 1000) then
            write (*,*) 'The problem is too expensive for RKC.'
            write (*,*) 'Alternative methods should be considered.'
         end if
         if (idid .eq. 3) go to 200
         if (idid .ne. 1) go to 300
c------ Successful iteration --------------------------------------
         n = n + 1
         repeated = .false.
         do 20 i = 1, neq
 20         yout(n, i) = yin(i)
         if (nsteps .ge. 5000) then
            write (*,*) ' Quit because of too much work.'
            go to 300
         end if
C*** Function TERMINATE returns with 1. 
C*** Stop iteration and return with current values.
         if (itermin .eq. 1) then
            call termin(yout, trr, n - 1, nt, neq, resterm)
            if (resterm .eq. 1) then
               idid = 0
               finish = .true.
            end if
         end if
         go to 400
C***  idid return with abnormal status
C***  When idid is 3, it indicates improper error control
C***  For some j, ATOL(j) = 0 and Y(j) = 0.
 200  write (*,*) 'Improper error control.'
      write (*,*) 'For some j, ATOL(j) = 0 and Y(j) = 0.'
      if (.not. repeated) then
         ptr = nint(work(7))
         do 30 i = 1, neq
            ycurrent(i) = work(ptr + i - 1) 
            if (ycurrent(i) .eq. 0.0) then
               if (info(4) .eq. 0) then 
                  write (*,*) 'ATOL is reset to 1e-8.'
                  write (*,*) 'Iteration is restarted.'
                  atol(1) = 1d-8
               else if (atol(i) .eq. 0.0) then
                  write (*,*) 'ATOL(', j, ')is reset to 1e-8.'
                  write (*,*) 'Iteration is restarted.'
                  atol(i) = 1d-8
               end if
            end if
 30      continue
         repeated = .true.
      else
         write (*,*) 'Quit because of repeated improper error-control.'
         go to 300
      end if
      go to 400
C*** Stop iteration and return     
 300     finish = .true.
 400  continue
      end do
      nstop = n
      return
      end subroutine

