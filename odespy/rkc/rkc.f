      subroutine rkc(spcrad,neqn,f,y,t,tend,rtol,atol,info,work,idid)
c--------------------------------------------------------------------------
c
c  ABSTRACT:  RKC integrates initial value problems for systems of first
c  order ordinary differential equations.  It is based on a family of
c  explicit Runge-Kutta-Chebyshev formulas of order two.  The stability
c  of members of the family increases quadratically in the number of
c  stages m. An estimate of the spectral radius is used at each step to
c  select the smallest m resulting in a stable integration. RKC is
c  appropriate for the solution to modest accuracy of mildly stiff problems
c  with eigenvalues of Jacobians that are close to the negative real axis.
c  For such problems it has the advantages of explicit one-step methods and
c  very low storage. If it should turn out that RKC is using m far beyond
c  100, the problem is not mildly stiff and alternative methods should be
c  considered.  Answers can be obtained cheaply anywhere in the interval
c  of integration by means of a continuous extension evaluated in the
c  subroutine RKCINT.
c
c  The initial value problems arising from semi-discretization of
c  diffusion-dominated parabolic partial differential equations and of
c  reaction-diffusion equations, especially in two and three spatial
c  variables, exemplify the problems for which RKC was designed.  Two
c  example programs, ExA and ExB, are provided that show how to use RKC.
c
c---------------------------------------------------------------------------
c  USAGE:  RKC integrates a system of NEQN first order ordinary differential
c  equations specified by a subroutine F from T to TEND.  The initial values
c  at T are input in Y(*).  On all returns from RKC, Y(*) is an approximate
c  solution at T.  In the computation of Y(*), the local error has been
c  controlled at each step to satisfy a relative error tolerance RTOL and
c  absolute error tolerances ATOL(*).  The array INFO(*) specifies the way
c  the problem is to be solved.  WORK(*) is a work array. IDID reports
c  success or the reason the computation has been terminated.
c
c  FIRST CALL TO RKC
c
c  You must provide storage in your calling program for the arrays in the
c  call list -- Y(NEQN), INFO(4), WORK(8+5*NEQN).  If INFO(2) = 0, you can
c  reduce the storage for the work array to WORK(8+4*NEQN).  ATOL may be
c  a scalar or an array.  If it is an array, you must provide storage for
c  ATOL(NEQN).  You must declare F in an external statement, supply the
c  subroutine F and the function SPCRAD, and initialize the following
c  quantities:
c
c    NEQN:  The number of differential equations. Integer.
c
c    T:     The initial point of the integration. Double precision.
c           Must be a variable.
c
c    TEND:  The end of the interval of integration.  Double precision.
c           TEND may be less than T.
c
c    Y(*):  The initial value of the solution.  Double precision array
c           of length NEQN.
c
c    F:     The name of a subroutine for evaluating the differential
c           equation.  It must have the form
c
c             subroutine f(neqn,t,y,dy)
c             integer          neqn
c             double precision t,y(neqn),dy(neqn)
c             dy(1)    = ...
c             ...
c             dy(neqn) = ...
c             return
c             end
c
c  RTOL,
C  ATOL(*):  At each step of the integration the local error is controlled
c            so that its RMS norm is no larger than tolerances RTOL, ATOL(*).
c            RTOL is a double precision scalar. ATOL(*) is either a double
c            precision scalar or a double precision array of length NEQN.
c            RKC is designed for the solution of problems to modest accuracy.
c            Because it is based on a method of order 2, it is relatively
c            expensive to achieve high accuracy.
c
c            RTOL is a relative error tolerance.  You must ask for some
c            relative accuracy, but you cannot ask for too much for the
c            precision available.  Accordingly, it is required that
c            0.1 >= RTOL >= 10*uround. (See below for the machine and
c            precision dependent quantity uround.)
c
c            ATOL is an absolute error tolerance that can be either a
c            scalar or an array.  When it is an array, the tolerances are
c            applied to corresponding components of the solution and when
c            it is a scalar, it is applied to all components.  A scalar
c            tolerance is reasonable only when all solution components are
c            scaled to be of comparable size.  A scalar tolerance saves a
c            useful amount of storage and is convenient.  Use INFO(*) to
c            tell RKC whether ATOL is a scalar or an array.
c
c            The absolute error tolerances ATOL(*) must satisfy ATOL(i) >= 0
c            for i = 1,...,NEQN.  ATOL(j)= 0 specifies a pure relative error
c            test on component j of the solution, so it is an error if this
c            component vanishes in the course of the integration.
c
c            If all is going well, reducing the tolerances by a factor of
c            0.1 will reduce the error in the computed solution by a factor
c            of roughly 0.2.
c
c  INFO(*)   Integer array of length 4 that specifies how the problem
c            is to be solved.
c
c  INFO(1):  RKC integrates the initial value problem from T to TEND.
c            This is done by computing approximate solutions at points
c            chosen automatically throughout [T, TEND].  Ordinarily RKC
c            returns at each step with an approximate solution. These
c            approximations show how y behaves throughout the interval.
c            The subroutine RKCINT can be used to obtain answers anywhere
c            in the span of a step very inexpensively. This makes it
c            possible to obtain answers at specific points in [T, TEND]
c            and to obtain many answers very cheaply when attempting to
c            locating where some function of the solution has a zero
c            (event location).  Sometimes you will be interested only in
c            a solution at TEND, so you can suppress the returns at each
c            step along the way if you wish.
c
c  INFO(1)  = 0 Return after each step on the way to TEND with a
c               solution Y(*) at the output value of T.
c
c           = 1 Compute a solution Y(*) at TEND only.
c
c  INFO(2):  RKC needs an estimate of the spectral radius of the Jacobian.
c            You must provide a function that must be called SPCRAD and
c            have the form
c
c              double precision function spcrad(neqn,t,y)
c              integer          neqn
c              double precision t,y(neqn)
c
c              spcrad = < expression depending on info(2) >
c
c              return
c              end
c
c            You can provide a dummy function and let RKC compute the
c            estimate. Sometimes it is convenient for you to compute in
c            SPCRAD a reasonably close upper bound on the spectral radius,
c            using, e.g., Gershgorin's theorem.  This may be faster and/or
c            more reliable than having RKC compute one.
c
c  INFO(2)  = 0  RKC is to compute the estimate internally.
c                Assign any value to SPCRAD.
c
c           = 1  SPCRAD returns an upper bound on the spectral
c                radius of the Jacobian of f at (t,y).
c
c  INFO(3):  If you know that the Jacobian is constant, you should say so.
c
c  INFO(3)  = 0  The Jacobian may not be constant.
c
c           = 1  The Jacobian is constant.
c
c  INFO(4):  You must tell RKC whether ATOL is a scalar or an array.
c
c  INFO(4)  = 0  ATOL is a double precision scalar.
c
c           = 1  ATOL is a double precision array of length NEQN.
c
c  WORK(*):  Work array.  Double precision array of length at least
c            8 + 5*NEQN if INFO(2) = 0 and otherwise, 8 + 4*NEQN.
c
c  IDID:     Set IDID = 0 to initialize the integration.
c
c
c
c  RETURNS FROM RKC
c
c  T:        The integration has advanced to T.
c
c  Y(*):     The solution at T.
c
c  IDID:     The value of IDID reports what happened.
c
c                          SUCCESS
c
c    IDID     = 1 T = TEND, so the integration is complete.
c
c             = 2 Took a step to the output value of T.  To continue on
c                 towards TEND, just call RKC again.   WARNING:  Do not
c                 alter any argument between calls.
c
c                 The last step, HLAST, is returned as WORK(1). RKCINT
c                 can be used to approximate the solution anywhere in
c                 [T-HLAST, T] very inexpensively using data in WORK(*).
c
c                 The work can be monitored by inspecting data in RKCDID.
c
c                          FAILURE
c
c             = 3 Improper error control: For some j, ATOL(j) = 0
c                 and Y(j) = 0.
c
c             = 4 Unable to achieve the desired accuracy with the
c                 precision available.  A severe lack of smoothness in
c                 the solution y(t) or the function f(t,y) is likely.
c
c             = 5 Invalid input parameters:  NEQN <= 0, RTOL > 0.1,
c                 RTOL < 10*UROUND, or ATOL(i) < 0 for some i.
c
c             = 6 The method used by RKC to estimate the spectral
c                 radius of the Jacobian failed to converge.
c
c  RKCDID is a labelled common block that communicates statistics
c         about the integration process:
c         common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm
c
c         The integer counters are:
c
c        NFE      number of evaluations of F used
c                   to integrate the initial value problem
c        NSTEPS   number of integration steps
c        NACCPT   number of accepted steps
c        NREJCT   number of rejected steps
c        NFESIG   number of evaluations of F used
c                   to estimate the spectral radius
c        MAXM     maximum number of stages used
c
c        This data can be used to monitor the work and terminate a run
c        that proves to be unacceptably expensive.  Also, if MAXM should
c        be far beyond 100, the problem is too expensive for RKC and
c        alternative methods should be considered.
c
c--------------------------------------------------------------------------
c
c   CAUTION: MACHINE/PRECISION ISSUES
c
c     UROUND (the machine precision) is the smallest number such that
c     1 + UROUND > 1, where 1 is a floating point number in the working
c     precision. UROUND is set in a parameter statement in RKC. Its
c     value depends on both the precision and the machine used, so it
c     must be set appropriately.  UROUND is the only constant in RKC
c     that depends on the precision.
c
c     This version of RKC is written in double precision. It can be changed
c     to single precision by replacing DOUBLE PRECISION in the declarations
c     by REAL and changing the type of the floating point constants set in
c     PARAMETER statements from double precision to real.
c
c--------------------------------------------------------------------------
c
c  Authors: B.P. Sommeijer and J.G. Verwer
c           Centre for Mathematics and Computer Science (CWI)
c           Kruislaan 413
c           1098 SJ  Amsterdam
c           The Netherlands
c           e-mail: bsom@cwi.nl
c
c           L.F. Shampine
c           Mathematics Department
c           Southern Methodist University
c           Dallas, Texas 75275-0156
c           USA
c           e-mail: lshampin@mail.smu.edu
c
c  Details of the methods used and the performance of RKC can be
c  found in
c
c         B.P. Sommeijer, L.F. Shampine and J.G. Verwer
c         RKC: an Explicit Solver for Parabolic PDEs.
c              Technical Report MAS-R9715, CWI, Amsterdam, 1997
c
c  This source code for RKC and some examples, as well as the
c  reference solution to the second example can also be obtained
c  by anonymous ftp from the address ftp://ftp.cwi.nl/pub/bsom/rkc
c------------------------------------------------------------------
      integer          neqn,info(*),idid
      double precision y(neqn),t,tend,rtol,atol(*),work(*)
c
c*********************************************************************
c  uround is set here for IEEE double precision arithmetic.
      double precision uround
      parameter       (uround=2.22d-16)
c*********************************************************************
c
      double precision zero,rmax,rmin
      parameter       (zero=0d0,rmax=0.1d0,rmin=10d0*uround)
      integer          i,ptr1,ptr2,ptr3,ptr4
      logical          array,valid
      save
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm
      external         f
c
      if(idid .eq. 0) then
c----------------------
c  Test the input data.
c----------------------
         array = info(4) .eq. 1
         valid = neqn .gt. 0
         if((rtol .gt. rmax) .or. (rtol .lt. rmin)) valid = .false.
         if(atol(1) .lt. zero) valid = .false.
         if(array) then
           do 10 i = 2, neqn
             if(atol(i) .lt. zero) valid = .false.
10         continue
         endif
         if(.not. valid) then
           idid = 5
           return
         endif
c-----------------------------------
c  Initialize counters and pointers.
c-----------------------------------
         nfe = 0
         nsteps = 0
         naccpt = 0
         nrejct = 0
         nfesig = 0
         maxm = 0
c-----------------------------------------------------------
c  work(*) contains information needed for interpolation,
c  continuation after a return, and working storage. Items
c  relevant here are:
c
c  The last step taken, hlast, is work(1).
c  The current t is work(2).
c  The number of equations, neqn, is work(3).
c  The unit roundoff, uround, is work(4).
c  The square root of uround, sqrtu, is work(5).
c  The maximum step size, hmax, is work(6).
c  The base address for the solution is ptr1 = nint(work(7)).
c  The solution at t starts at ptr1.
c  The derivative of the solution at t starts at ptr2.
c  The solution at t-hlast starts at ptr3.
c  The derivative of the solution at t-hlast starts at ptr4.
c  The estimated dominant eigenvector starts at ptr4 + neqn.
c------------------------------------------------------------
        work(2) = t
        work(3) = neqn
        work(4) = uround
        work(5) = sqrt(uround)
        ptr1 = 8
        work(7) = ptr1
        ptr2 = ptr1 + neqn
        ptr3 = ptr2 + neqn
        ptr4 = ptr3 + neqn
      elseif(idid .ne. 2) then
        write(*,*) ' RKC was called with an illegal value of IDID.'
        stop
      endif
c
      call rkclow(spcrad,neqn,t,tend,y,f,info,rtol,atol,work,
     &            work(ptr1),work(ptr2),work(ptr3),work(ptr4),idid)
      return
      end

      subroutine rkclow(spcrad,neqn,t,tend,y,f,info,rtol,atol,work,
     &                  yn,fn,vtemp1,vtemp2,idid)
c----------------------------------------------------------------------
c  RKC is an interface to RKCLOW where the actual solution takes place.
c----------------------------------------------------------------------
      integer          neqn,info(*),idid
      double precision t,tend,y(*),rtol,atol(*),work(*),
     &                 yn(*),fn(*),vtemp1(*),vtemp2(*)
      external         f, spcrad
c
      double precision one,onep1,onep54,p1,p4,p8,
     &                 ten,zero,one3rd,two3rd
      parameter       (one=1d0,onep1=1.1d0,onep54=1.54d0,
     &                 p1=0.1d0,p4=0.4d0,p8=0.8d0,ten=10d0,
     &                 zero=0d0,one3rd=1d0/3d0,two3rd=2d0/3d0)
      integer          i,m,mmax,nstsig
      double precision absh,est,err,errold,fac,h,hmax,hmin,hold,
     &                 spcrad,sprad,tdir,temp1,temp2,
     &                 uround,wt,ylast,yplast,at
      logical          array,last,newspc,jacatt
      save
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm
c
c---------------------------------
c    Initialize on the first call.
c---------------------------------
      if(idid .eq. 0) then
        array = info(4) .eq. 1
        uround = work(4)
        mmax = nint(sqrt(rtol/(10d0*uround)))
        mmax = max(mmax,2)
        newspc = .true.
        jacatt = .false.
        nstsig = 0
        do 10 i = 1, neqn
          yn(i) = y(i)
10      continue
        call f(neqn,t,yn,fn)
        nfe = nfe + 1
        tdir = sign(one,tend - t)
        hmax = abs(tend - t)
        work(6) = hmax
        hmin = ten*uround*max(abs(t),hmax)
      endif
c------------------------------------
c  Start of loop for taking one step.
c------------------------------------
20    continue
c----------------------------------------------
c  Estimate the spectral radius of the Jacobian
c  when newspc = .true..  A convergence failure
c  in rkcrho is reported by idid = 6.
c----------------------------------------------
      if(newspc) then
        if(info(2) .eq. 1) then
          sprad = spcrad(neqn,t,yn)
        else
          call rkcrho(neqn,t,f,yn,fn,vtemp1,vtemp2,work,sprad,idid)
          if(idid .eq. 6) return
        endif
        jacatt = .true.
      endif

c-------------------------------
c  Compute an initial step size.
c-------------------------------
      if(nsteps .eq. 0) then
        absh = hmax
        if(sprad*absh .gt. one) absh = one/sprad
        absh = max(absh,hmin)
        do 30 i = 1,neqn
          vtemp1(i) = yn(i) + absh*fn(i)
30      continue
        call f(neqn,t+absh,vtemp1,vtemp2)
        nfe = nfe + 1
        est = zero
        at = atol(1)
        do 40 i = 1,neqn
          if(array) at = atol(i)
          wt = at + rtol*abs(yn(i))
          if(wt .eq. zero) then
            idid = 3
            return
          endif
          est = est + ((vtemp2(i) - fn(i))/wt)**2
40      continue
        est = absh*sqrt(est/neqn)
        if(p1*absh .lt. hmax*sqrt(est)) then
          absh = max(p1*absh/sqrt(est), hmin)
        else
          absh = hmax
        endif
      endif
c------------------------------------------------------------
c  Adjust the step size and determine the number of stages m.
c------------------------------------------------------------
      last = .false.
      if(onep1*absh .ge. abs(tend - t)) then
        absh = abs(tend - t)
        last = .true.
      endif
      m = 1 + int(sqrt(onep54*absh*sprad + one))
c----------------------------------------------------------
c  Limit m to mmax to control the growth of roundoff error.
c----------------------------------------------------------
      if(m .gt. mmax) then
        m = mmax
        absh = (m**2 - 1)/(onep54*sprad)
        last = .false.
      endif
      maxm = max(m,maxm)
c--------------------------------------------
c  A tentative solution at t+h is returned in
c  y and its slope is evaluated in vtemp1(*).
c--------------------------------------------
      h = tdir*absh
      hmin = ten*uround*max(abs(t),abs(t + h))
      call step(neqn,f,t,yn,fn,h,m,y,vtemp1,vtemp2)
      call f(neqn,t+h,y,vtemp1)
      nfe = nfe + m
      nsteps = nsteps + 1
c-------------------------------------------------------------
c  Estimate the local error and compute its weighted RMS norm.
c-------------------------------------------------------------
      err = zero
      at = atol(1)
      do 50 i = 1, neqn
        if(array) at = atol(i)
        wt = at + rtol*max(abs(y(i)),abs(yn(i)))
        if(wt .eq. zero) then
          idid = 3
          return
        endif
        est = p8*(yn(i) - y(i)) + p4*h*(fn(i) + vtemp1(i))
        err = err + (est/wt)**2
50    continue
      err = sqrt(err/neqn)
c
      if(err .gt. one) then
c-------------------
c  Step is rejected.
c-------------------
        nrejct = nrejct + 1
        absh = p8*absh/(err**one3rd)
        if(absh .lt. hmin) then
          idid = 4
          return
        else
          newspc = .not. jacatt
          goto 20
        endif
      endif
c-------------------
c  Step is accepted.
c-------------------
      naccpt = naccpt + 1
      t = t + h
      jacatt = info(3) .eq. 1
      nstsig = mod(nstsig+1,25)
      newspc = .false.
      if((info(2) .eq. 1) .or. (nstsig .eq. 0)) newspc = .not. jacatt
c------------------------------------------------------
c  Update the data for interpolation stored in work(*).
c------------------------------------------------------
      work(1) = h
      work(2) = t
      do 60 i = 1, neqn
         ylast = yn(i)
         yplast = fn(i)
         yn(i) = y(i)
         fn(i) = vtemp1(i)
         vtemp1(i) = ylast
         vtemp2(i) = yplast
60    continue
      fac = ten
      if(naccpt .eq. 1) then
        temp2 = err**one3rd
        if(p8 .lt. fac*temp2) fac = p8/temp2
      else
        temp1 = p8*absh*errold**one3rd
        temp2 = abs(hold)*err**two3rd
        if(temp1 .lt. fac*temp2) fac = temp1/temp2
      endif
      absh = max(p1,fac)*absh
      absh = max(hmin,min(hmax,absh))
      errold = err
      hold = h
      h = tdir*absh
      if(last) then
        idid = 1
        return
      elseif(info(1) .eq. 0) then
        idid = 2
        return
      else
        goto 20
      endif
      end

      subroutine step(neqn,f,t,yn,fn,h,m,y,yjm1,yjm2)
c--------------------------------------------------
c  Take a step of size H from T to T+H to get Y(*).
c--------------------------------------------------
      integer          neqn,m
      double precision t,yn(neqn),fn(neqn),h,
     &                 y(neqn),yjm1(neqn),yjm2(neqn)
      external         f
c
      double precision one,two,four,c13,zero
      parameter       (one=1d0,two=2d0,four=4d0,c13=13d0,zero=0d0)
      integer          i,j
      double precision ajm1,arg,bj,bjm1,bjm2,dzj,dzjm1,dzjm2,
     &                 d2zj,d2zjm1,d2zjm2,mu,mus,nu,
     &                 temp1,temp2,thj,thjm1,thjm2,w0,w1,
     &                 zj,zjm1,zjm2
c
      w0 = one + two/(c13*m**2)
      temp1 = w0**2 - one
      temp2 = sqrt(temp1)
      arg = m*log(w0 + temp2)
      w1 = sinh(arg)*temp1 / (cosh(arg)*m*temp2 - w0*sinh(arg))
      bjm1 = one/(two*w0)**2
      bjm2 = bjm1
c---------------------------
c  Evaluate the first stage.
c---------------------------
      do 10 i = 1, neqn
        yjm2(i) = yn(i)
10    continue
      mus = w1*bjm1
      do 20 i = 1, neqn
        yjm1(i) = yn(i) + h*mus*fn(i)
20    continue
      thjm2  = zero
      thjm1  = mus
      zjm1   = w0
      zjm2   = one
      dzjm1  = one
      dzjm2  = zero
      d2zjm1 = zero
      d2zjm2 = zero
c------------------------------
c  Evaluate stages j = 2,...,m.
c------------------------------
      do 50 j = 2, m
        zj   =   two*w0*zjm1 - zjm2
        dzj  =   two*w0*dzjm1 - dzjm2 + two*zjm1
        d2zj =   two*w0*d2zjm1 - d2zjm2 + four*dzjm1
        bj   =   d2zj/dzj**2
        ajm1 =   one - zjm1*bjm1
        mu   =   two*w0*bj/bjm1
        nu   = - bj/bjm2
        mus  =   mu*w1/w0
c---------------------------------------------
c  Use the y array for temporary storage here.
c---------------------------------------------
        call f(neqn,t + h*thjm1,yjm1,y)
        do 30 i = 1, neqn
          y(i) = mu*yjm1(i) + nu*yjm2(i) + (one - mu - nu)*yn(i) +
     &           h*mus*(y(i) - ajm1*fn(i))
30      continue
        thj = mu*thjm1 + nu*thjm2 + mus*(one - ajm1)
c------------------------------------
c  Shift the data for the next stage.
c------------------------------------
        if(j .lt. m) then
          do 40 i = 1, neqn
            yjm2(i) = yjm1(i)
            yjm1(i) = y(i)
40        continue
          thjm2  = thjm1
          thjm1  = thj
          bjm2   = bjm1
          bjm1   = bj
          zjm2   = zjm1
          zjm1   = zj
          dzjm2  = dzjm1
          dzjm1  = dzj
          d2zjm2 = d2zjm1
          d2zjm1 = d2zj
        endif
50    continue
      return
      end

      subroutine rkcint(work,arg,yarg)
c-------------------------------------------------------------------------
c  RKCINT is used to compute approximate solutions at specific t and to
c  compute cheaply the large number of approximations that may be needed
c  for plotting or locating when events occur.
c
c  After a step to T, RKC provides HLAST, the step just taken, in WORK(1).
c  In other entries of WORK(*) it provides the data needed to interpolate
c  anywhere in [T-HLAST, T]. YARG(*), the approximate solution at t = ARG
c  computed by interpolation in RKCINT has the same order of accuracy as
c  the Y(*) computed directly by RKC.
c
c  INPUT:
c
c    WORK(*)   Double precision array returned by RKC.
c
c    ARG       The point at which a solution is desired. Double precision.
c
c  OUTPUT:
c
c    YARG(*)   The approximate solution at t = ARG.  Double precision
c              array of length neqn.
c--------------------------------------------------------------------------
      double precision work(*),arg,yarg(*)
c
      double precision one,two,three
      parameter       (one=1d0,two=2d0,three=3d0)
      integer          i,neqn,ptr1,ptr2,ptr3,ptr4
      double precision a1,a2,b1,b2,s,hlast,t,tlast
c
c---------------------------------------------------------------------
c  The data needed for interpolation are stored in work(*) as follows:
c
c  The last step taken, hlast, is work(1).
c  The current t is work(2).
c  The number of equations, neqn, is work(3).
c  The base address for the solution is ptr1 = nint(work(7))
c  The solution at t starts at ptr1.
c  The derivative of the solution at t starts at ptr2.
c  The solution at t-hlast starts at ptr3.
c  The derivative of the solution at t-hlast starts at ptr4.
c---------------------------------------------------------------------
      hlast = work(1)
      t = work(2)
      tlast = t - hlast
      neqn = nint(work(3))
      ptr1 = nint(work(7))
      ptr2 = ptr1 + neqn
      ptr3 = ptr2 + neqn
      ptr4 = ptr3 + neqn
c
      s  = (arg - tlast)/hlast
      a1 = (one + two*s)*(s - one)**2
      a2 = (three - two*s)*s**2
      b1 = hlast*s*(s - one)**2
      b2 = hlast*(s - one)*s**2
c
      do 10 i = 1, neqn
        yarg(i) = a1*work(ptr3+i-1) + a2*work(ptr1+i-1) +
     &            b1*work(ptr4+i-1) + b2*work(ptr2+i-1)
10    continue
      return
      end

      subroutine rkcrho(neqn,t,f,yn,fn,v,fv,work,sprad,idid)
c---------------------------------------------------------------
c  RKCRHO attempts to compute a close upper bound, SPRAD, on
c  the spectral radius of the Jacobian matrix using a nonlinear
c  power method.  A convergence failure is reported by IDID = 6.
c---------------------------------------------------------------
      integer          neqn,idid
      double precision t,yn(neqn),fn(neqn),v(neqn),fv(neqn),work(*),
     &                 sprad
      external         f
c
      integer          itmax
      parameter       (itmax=50)
      double precision zero,one,onep2,p01
      parameter       (zero=0d0,one=1d0,onep2=1.2d0,p01=0.01d0)
      integer          i,iter,index,ptr5
      double precision uround,sqrtu,ynrm,sigma,sigmal,
     &                 dynrm,dfnrm,vnrm,small
      integer          nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/  nfe,nsteps,naccpt,nrejct,nfesig,maxm
c
      uround = work(4)
      sqrtu = work(5)
c------------------------------------------------------------
c  hmax = work(6).  sprad smaller than small = 1/hmax are not
c  interesting because they do not constrain the step size.
c------------------------------------------------------------
      small = one/work(6)
c---------------------------------------------------------
c  The initial slope is used as guess when nsteps = 0 and
c  thereafter the last computed eigenvector.  Some care
c  is needed to deal with special cases. Approximations to
c  the eigenvector are normalized so that their Euclidean
c  norm has the constant value dynrm.
c---------------------------------------------------------
      ptr5 = nint(work(7)) + 4*neqn
      if(nsteps .eq. 0) then
        do 10 i = 1,neqn
          v(i) = fn(i)
10      continue
      else
        do 20 i = 1,neqn
          v(i) = work(ptr5+i-1)
20      continue
      endif
      ynrm = zero
      vnrm = zero
      do 30 i = 1,neqn
        ynrm = ynrm + yn(i)**2
        vnrm = vnrm + v(i)**2
30    continue
      ynrm = sqrt(ynrm)
      vnrm = sqrt(vnrm)
      if(ynrm .ne. zero .and. vnrm .ne. zero) then
        dynrm = ynrm*sqrtu
        do 40 i = 1,neqn
          v(i) = yn(i) + v(i)*(dynrm/vnrm)
40      continue
      elseif(ynrm .ne. zero) then
        dynrm = ynrm*sqrtu
        do 50 i = 1, neqn
          v(i) = yn(i) + yn(i)*sqrtu
50      continue
      elseif(vnrm .ne. zero) then
        dynrm = uround
        do 60 i = 1,neqn
          v(i) = v(i)*(dynrm/vnrm)
60      continue
      else
        dynrm = uround
        do 70 i = 1,neqn
          v(i) = dynrm
70      continue
      endif
c--------------------------------------------
c  Now iterate with a nonlinear power method.
c--------------------------------------------
      sigma = zero
      do 110 iter = 1, itmax
        call f(neqn,t,v,fv)
        nfesig = nfesig + 1
        dfnrm = zero
        do 80 i = 1, neqn
          dfnrm = dfnrm + (fv(i) - fn(i))**2
80      continue
        dfnrm = sqrt(dfnrm)
        sigmal = sigma
        sigma = dfnrm/dynrm
c----------------------------------------------------------
c  sprad is a little bigger than the estimate sigma of the
c  spectral radius, so is more likely to be an upper bound.
c----------------------------------------------------------
        sprad = onep2*sigma
        if(iter .ge. 2 .and.
     &     abs(sigma - sigmal) .le. max(sigma,small)*p01) then
          do 90 i = 1,neqn
            work(ptr5+i-1) = v(i) - yn(i)
90        continue
          return
        endif
c--------------------------------------
c  The next v(*) is the change in f
c  scaled so that norm(v - yn) = dynrm.
c--------------------------------------
        if(dfnrm .ne. zero) then
          do 100 i = 1,neqn
            v(i) = yn(i) + (fv(i) - fn(i))*(dynrm/dfnrm)
100       continue
        else
c-------------------------------------------------------
c  The new v(*) degenerated to yn(*)--"randomly" perturb
c  current approximation to the eigenvector by changing
c  the sign of one component.
c-------------------------------------------------------
          index = 1 + mod(iter,neqn)
          v(index) = yn(index) - (v(index) - yn(index))
        endif
110   continue
c-------------------------------------------
c  Set flag to report a convergence failure.
c-------------------------------------------
      idid = 6
      return
      end
