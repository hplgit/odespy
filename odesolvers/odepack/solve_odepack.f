      subroutine solve(termin, itermin, nt, nstart, nstop, f,  
     1                 neq, yin, trr, itol, rtol, atol,   
     2                 istate, iopt, rworkin, lrwin, lrw, iworkin,
     3                 liwin, liw, jac, jaccolumn, mf, g, 
     6                 ng, rinfo, iinfo, subroutinename)
C***Input parameters:
C
C     external functions: 
C
C       termin    -- function to detect stop event.
c                 termin(yrr,trr,step_no,nt,neq,resterm) --> resterm = 0/1
C       f         -- function to define dy/dt in ODE system.
c                 f(neq,t,y,dy) --> dy
C       jac       -- function to define jacobian matrix in Lsode, Lsoda, Lsodar
C                 jac(neq,t,y,ml,mu,pd,nrowpd)        -->  df/dy
C       jaccolumn -- a specified column of Jacobian (df/du) matrix with 
C                    arbitrary sparse structure in Lsodes.
C                 jac(neq,t,y,j,ia,ja,pdj)  --> pdj = jth column in Jacobian
C       g         -- defines the equations for root-finding.
C
C     integer inputs: 
C       itermin   -- 0/1, whether TERMIN is defined.
C       iopt      -- 0/1, whether optional inputs are involved.
C       itol      -- 1/2/3/4, indicates the type of error control.
C       istate    -- 0/1/2/3, indicates the state of calculation.
C       nstart    -- the step number to start iteration.
C       lrw       -- Estimated length requirement for real work array.
C       liw       -- Estimated length requirement for integer work array.
C       ng        -- how many equations are involved in root-finding.
C       iwork_in  -- array to hold optional integer inputs.
C       mf        -- indicate method choice
C
C     Real inputs:
C       yin       -- complete array to hold y before current iteration.
C                    dimension(nstart, neq)
C       trr       -- time-points. 
C       rtol      -- relative tolerance.
C       atol      -- absolute tolerance.
C       rwork_in  -- array to hold optional float inputs.

C     Other    :   
C       solver_name-- String for subroutine name in origianl ODEPACK package.

C***Output parameters:
C     nstop    : Step number for current successful iteration.
C     yin      : 2d-array for Y for current successful iteration. 
C     istate   : current status.
C     rinfo    : Optional float outputs.
C     iinfo    : Optional integer outputs.
 
      external f, jac, jaccolumn, termin, g
      integer itermin, neq, itol, istate, iopt, lrw, liw, lrwin, liwin, 
     1        iworkin, mf, nt, nstart, nstop, iinfo, ng, yoti, nnz
      double precision yin, trr, rtol, atol, rworkin, rinfo
      dimension yin(nt, neq), trr(nt), rtol(*), atol(*), 
     1        rworkin(lrwin), iworkin(liwin), iinfo(17 + ng), rinfo(5)
      character*(*) subroutinename
Cf2py intent(in, out) yrr, istate
Cf2py intent(out) nstop, rinfo, iinfo
Cf2py intent(in) termin, f, itol, rtol, atol, itermin, trr, itask, iopt, rwork
Cf2py intent(in) iwork, jac, mf, nstart, subroutinename, jaccolumn, g, ng
Cf2py intent(hide) neq, nt, lrw, liw

      external dlsode, dlsodes, dlsoda, dlsodar
      integer i, j, resterm, itask, neqnzz
      integer, allocatable :: jroot(:)
      integer, allocatable :: iwork(:)
      double precision, allocatable :: rwork(:)

      double precision y, t, tout
      dimension y(neq), neqnzz(2)
      logical finish

C***Initialization
      n = nstart
      resterm = 0
      do 10 i = 1, neq
 10      y(i) = yin(n,i)
      finish = .FALSE.
      itask = 1
      if (ng .gt. 0) allocate(jroot(ng))
C***Work arrays need to be initialized to required length
      allocate(iwork(liw))
      allocate(rwork(lrw))
      do 20 i = 1, liwin
 20      iwork(i) = iworkin(i)
      do 30 i = 1, lrwin 
 30      rwork(i) = rworkin(i)
C***Start while loop
      do while ((n .lt. nt) .and. (.not. finish))
         t = trr(n)
         tout = trr(n + 1)
         if (subroutinename .eq. 'dlsode') then
            call dlsode(f, neq, y, t, tout, itol, rtol, atol, itask,  
     1                  istate, iopt, rwork, lrw, iwork, liw, jac, mf)
         else if (subroutinename .eq. 'dlsodes') then
            neqnzz(1) = neq
            if (liwin .gt. 30) then
               neqnzz(2) = iworkin(31 + neq) - 1
            else
               neqnzz(2) = neq*neq/2
            end if
            call dlsodes(f, neq, y, t, tout, itol, rtol, atol, itask,  
     1                  istate, iopt, rwork, lrw, iwork, liw, jaccolumn,
     2                  mf)
         else if (subroutinename .eq. 'dlsoda') then
            call dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask,  
     1                  istate, iopt, rwork, lrw, iwork, liw, jac, mf)
         else if (subroutinename .eq. 'dlsodar') then
            if (ng .gt. 0) then
               call dlsodar(f, neq, y, t, tout, itol, rtol, atol, itask,  
     1                      istate, iopt, rwork, lrw, iwork, liw, jac,  
     2                      mf, g, ng, jroot)
            else 
               call dlsodar(f, neq, y, t, tout, itol, rtol, atol, itask,  
     1                      istate, iopt, rwork, lrw, iwork, liw, jac, 
     2                      mf, g, ng, jrootdummy)
            end if
         end if

C*** Error occurs. Pause iteration and return to python part.
         if (istate .lt. 0) go to 200
C*** Continue iteration after successful iteration.
         n = n + 1
         do 60 i = 1, neq
 60         yin(n, i) = y(i)

C*** Function TERMINATE returns with 1. 
C*** Stop iteration and return with current values.
         if (itermin .eq. 1) then
            call termin(yin, trr, n - 1, nt, neq, resterm)
            if (resterm .eq. 1) then
               istate = 0
               finish = .TRUE.
            end if
         end if
         if (istate .eq. 3) then
            if (subroutinename .eq. 'dlsodar') go to 300
            istate = -10
            go to 200
         end if 
         go to 500 
C***istate is a negative number, -- abnormal status
C***Stop iteration and return with extra output arrays
C***rinfo: optional output floats.
C***iinfo: optional output integers + jroot.
 200     do 70 i = 1, 5
 70         rinfo(i) = rwork(10 + i)
         do 80 i = 1, 17
 80         if (liw .ge. i + 9) iinfo(i) = iwork(9 + i)
         if (ng .eq. 0) go to 400
 300     do 90 i = 1, ng
 90         iinfo(17 + i) = jroot(i)
C***istate = 3, which indicates roots founded
C***Pause iteration and return with current values T & Y
 400     finish = .TRUE.
 500  continue
      end do
      nstop = n
      return
      end




