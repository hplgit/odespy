c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        auxiliary routines for all drivers
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/auxil/report.f
c
c     This is revision
c     $Id: report.F,v 1.20 2006/10/31 09:00:44 testset Exp $
c
c-----------------------------------------------------------------------
c
c     This source file contains the common subprograms used by
c     all Test Set driver programs.
c
c     To get CPU times you should configure gettim:
c     use your favorite editor to search for: ConfigureTiming
c     this is the only thing you need to configure.
c
c     However, if you wish to add problems / solvers / drivers
c     you have to configure more:
c     in that case use your favorite editor to search for:
c
c        ConfigureProblem
c        ConfigureSolver
c        ConfigureDriver
c
c-----------------------------------------------------------------------

ConfigureTiming: let gettim return CPU time
      real function gettim()
c
c gettim is supposed to deliver the cpu time (user time)
c beware of the clock resolution: typically 1/50 - 1/100 s
c
c in Fortran 77 timing is highly machine/compiler specific
c on some machines you might even have to use a C function instead
c

c dummy; returns 0
      write(*,*) 'WARNING: report.f: ',
     +           'currently no timing; activate timing by configuring',
     +           'the function gettim in the file report.f'
      gettim = 0
c
c if the compiler supports the intrinsic cpu_time
c     call cpu_time(gettim)
c
c best bet for UNIX
c
c     real*4 etime
c     real*4 total, tarray(2)
c     total = etime(tarray)
c     gettim = tarray(1)

c Cray
c
c     real second
c     gettim = second()

      return
      end

      subroutine chkdat(driver,solver,drivd)
      character*(*) driver, solver
      integer drivd
c
c     check date of problem interface
c
c     to be full proof we should check all release dates
c     i.e. of the driver, solver, problem, and of this file: report.f
c
c     however, the most important date: the release date of the solver
c     is not available in the solver
c     so we only do a sanity check on the problem interface date
c
      external lnblnk
      integer  pidate
      integer probd
      probd=pidate()

      probd=probd
      if (probd.ne.drivd) then
         write(*,*) 'ERROR: report.f: ',
     +              'unknown Test Set problem interface date', probd
         if (probd.gt.drivd) then
            write(*,*) 'you should probably get and ',
     +                 'install a more recent ',
     +                 driver(:lnblnk(driver)),'.f'
         else
            write(*,*) 'you should probably get and ',
     +                 'install a more recent problem'
         endif
         write(*,*) 'from http://www.dm.uniba.it/~testset/'
         stop
      endif

      return
      end

      subroutine getinp(driver,problm,solver,fullnm,
     +                  tolvec,rtol,atol,h0,solmax)
      character*(*) driver, problm, solver, fullnm
      double precision rtol(*), atol(*), h0, solmax
      logical tolvec
c
c     get input
c     check whether driver, problem, and solver are supported
c     read rtol, atol, h0, solmax
c
      integer i
      double precision tol
      logical error

      logical batch
      external batch,lnblnk

      character*3 FF, LF, PROMPT
      call getfmt(FF,LF,PROMPT)

      error  = .false.

ConfigureProblem: straightforward
      if (.not.(
     +   problm.eq.'chemakzo'  .or.
     +   problm.eq.'hires'     .or.
     +   problm.eq.'vdpol'     .or.
     +   problm.eq.'vdpolm'     .or.
     +   problm.eq.'rober'     .or.
     +   problm.eq.'orego'     .or.
     +   problm.eq.'e5'        .or.
     +   problm.eq.'pollu'     .or.
     +   problm.eq.'ringmod'   .or.
     +   problm.eq.'andrews'   .or.
     +   problm.eq.'transamp'  .or.
     +   problm.eq.'medakzo'   .or.
     +   problm.eq.'emep'      .or.
     +   problm.eq.'nand'      .or.
     +   problm.eq.'pump'      .or.
     +   problm.eq.'wheel'     .or.
     +   problm.eq.'tba'       .or.
     +   problm.eq.'caraxis'   .or.
     +   problm.eq.'fekete'    .or.
     +   problm.eq.'plei'      .or.
     +   problm.eq.'water'     .or.
     +   problm.eq.'crank'     .or.
     +   problm.eq.'beam'
     +   ))
     +then
         write(*,*) 'ERROR: report.f: unknown testset problem: ', problm
         error = .true.
      endif

ConfigureDriver: straightforward
      if (.not.(
     +   driver.eq.'bimdd'     .or.
     +   driver.eq.'dassld'    .or.
     +   driver.eq.'dopri5d'   .or.
     +   driver.eq.'gamdd'     .or.
     +   driver.eq.'mebdfd'    .or.
     +   driver.eq.'mebdfid'   .or.
     +   driver.eq.'psided'    .or.
     +   driver.eq.'psoded'    .or.
     +   driver.eq.'radaud'    .or.
     +   driver.eq.'radau5d'   .or.
     +   driver.eq.'voded'     
     +   ))
     +then
         write(*,*) 'WARNING: report.f: unknown testset driver: ', driver
         write(*,*)
      endif

ConfigureSolver: straightforward
      if (.not.(
     +   solver.eq.'BIMD'     .or.
     +   solver.eq.'DDASSL'    .or.
     +   solver.eq.'DOPRI5'   .or.
     +   solver.eq.'GAMD'     .or.
     +   solver.eq.'MEBDFDAE' .or.
     +   solver.eq.'MEBDFI'   .or.
     +   solver.eq.'PSIDE'    .or.
     +   solver.eq.'PSODE'    .or.
     +   solver.eq.'RADAU'    .or.
     +   solver.eq.'RADAU5'   .or.
     +   solver.eq.'VODE'
     +   ))
     +then
         write(*,*) 'ERROR: report.f: unknown testset solver: ', solver
         error = .true.
      endif

      if (error) stop

      write(*,'('//LF//',a,//,'//LF//',a,1x,a,1x,a,1x,a,/)')
     +        ' Test Set for IVP Solvers (release 2.3)',
     +        ' Solving',
     +        fullnm(:lnblnk(fullnm)),
     +        'using',
     +        solver(:lnblnk(solver))

      tolvec  = .false.
      rtol(1) = 0d0
      atol(1) = 0d0
      h0      = 0d0
      solmax  = 0d0

      if (.not. batch()) then

ConfigureSolver: read rtol/atol
ConfigureProblem: read special tolerance
         write(*,'('//LF//',a)') 'User input:'
         write(*,*)
         if (problm .eq. 'pump') then
            write(*,'('//LF//',a'//PROMPT//')')
     +              'give (pump) error tolerance: '
            read *, tol
            rtol(1) = tol
            atol(1) = tol
         elseif  (solver .eq. 'PSODE') then
            write(*,'('//LF//',a'//PROMPT//')')
     +              'give (mixed) error tolerance: '
            read *, tol
            rtol(1) = tol
            atol(1) = tol
         else
            write(*,'('//LF//',a'//PROMPT//')')
     +              'give relative error tolerance: '
            read *, rtol(1)
            write(*,'('//LF//',a'//PROMPT//')')
     +              'give absolute error tolerance: '
            read *, atol(1)
         endif
ConfigureSolver: read h0, only if used by solver
         if (
     +      solver.eq.'BIMD'     .or.
     +      solver.eq.'GAMD'     .or.
     +      solver.eq.'MEBDFDAE' .or.
     +      solver.eq.'MEBDFI'   .or.
     +      solver.eq.'PSODE'    .or.
     +      solver.eq.'RADAU'    .or.
     +      solver.eq.'RADAU5'
     +      )
     +   then
            write(*,'('//LF//',a'//PROMPT//')')
     +              'give initial stepsize: '
            read *, h0
         endif
ConfigureSolver: read solmax, only if used by solver
         if (solver.eq.'PSODE') then
            write(*,'('//LF//',a'//PROMPT//')')
     +              'give solmax: '
            read *, solmax
         endif

      else

ConfigureSolver: generally CWI only (i.e. normally batch() .eq. .false.)
ConfigureSolver: read rtol, atol, h0, solmax; same as above
         if (solver .eq. 'PSODE') then
            read *, rtol(1), atol(1), h0, solmax
         elseif (
     +      solver.eq.'BIMD'     .or.
     +      solver.eq.'GAMD'     .or.
     +      solver.eq.'MEBDFDAE' .or.
     +      solver.eq.'MEBDFI'   .or.
     +      solver.eq.'RADAU'    .or.
     +      solver.eq.'RADAU5'
     +   ) then
            read *, rtol(1), atol(1), h0
         else
            read *, rtol(1), atol(1)
         endif

      endif

c     a hack to allow input of e.g. rtol = 10d0 ** (-1.5d0) exactly
c     negative input is assumed to be the log10 value

      if (atol(1).lt.0d0) atol(1)=10d0**atol(1)
      if (rtol(1).lt.0d0) rtol(1)=10d0**rtol(1)
      if (h0     .lt.0d0) h0     =10d0**h0
      if (solmax .lt.0d0) solmax =10d0**solmax

      if (batch()) then

ConfigureSolver: generally CWI only (i.e. normally batch() .eq. .false.)
ConfigureSolver: print rtol, atol, h0, solmax; see read above
         if (solver .eq. 'PSODE') then
            write(*,'('//LF//',a,e12.3)')
     +         'with (mixed) error tolerance:  ', rtol(1)
         else
            write(*,'('//LF//',a,e12.3)')
     +         'with relative error tolerance: ', rtol(1)
            write(*,'('//LF//',a,e12.3)')
     +         '     absolute error tolerance: ', atol(1)
         endif
         if (
     +      solver.eq.'BIMD'     .or.
     +      solver.eq.'GAMD'     .or.
     +      solver.eq.'MEBDFDAE' .or.
     +      solver.eq.'MEBDFI'   .or.
     +      solver.eq.'PSODE'    .or.
     +      solver.eq.'RADAU'    .or.
     +      solver.eq.'RADAU5'
     +   ) then
            write(*,'('//LF//',a,e12.3)')
     +              '     initial stepsize:         ', h0
         endif
         if (solver.eq.'PSODE') then
            write(*,'('//LF//',a,e12.3)')
     +              '     solmax:                   ',
     +              solmax
         endif
      endif

ConfigureProblem: special (vector) tolerances
c This is now done in the function settolerances in the 
c definition of each problem
      

      return
      end

      subroutine getscd(mescd,scd,neqn,yref,y,problm,tolvec,atol,rtol,
     +                  printout)
      integer neqn
      double precision mescd, scd, yref(neqn), y(neqn)
      character*(*) problm     
      logical tolvec, printout
      double precision rtol(neqn), atol(neqn)
c
c     compute the accuracy of the numerical solution
c
      integer i
      double precision aerr, aerrmx, aerrl2, rerr, rerrmx, rerrl2
      double precision mixed_rerr, mixed_rerrmx, min_value
      double precision scda,scdr
      logical skipi, skipime
      character*1 skipc, skipcme
      integer numa, numr

      character*3 FF, LF, PROMPT

c
c     the double precision ieee rounding unit
c
      double precision uround
      parameter (uround=1.D-16)

      call getfmt(FF,LF,PROMPT)

      numa   = 0
      numr   = 0
      numme  = 0
      aerrl2 = 0d0
      rerrl2 = 0d0 
      aerrmx = 0d0
      rerrmx = 0d0
      mixed_rerrmx = 0d0
      min_value = 1.0d-200

cwe need the tolerances to computes the mixed error
      
      if (.not.tolvec) then
         do i=2,neqn
            atol(i)=atol(1)
            rtol(i)=rtol(1)
         enddo
      endif 
      
      if (printout) then
      write(*,'('//LF//',a)') 'Numerical solution:'
      write(*,*)
      write(*,'('//LF//',a,t41,a)')
     +' '                                  ,
     + '            scd                    ',
     +'       solution component'        ,
     + '---------------------------     ignore      ',
     +' '                                ,
     + 'mixed       abs        rel    mix - abs,rel',
     +'----------------------------------           ', 
     +'-----      -----      -----   -------------   '

      endif
      do 10 i=1,neqn

ConfigureProblem: skipping components in scd computation
         skipi =
     +           problm.eq.'emep'   .and. (i.eq.36 .or. i.eq.38)
     +      .or. problm.eq.'crank'   .and. (i .ge. 8 .and. i .le. 21)
     +      .or. problm.eq.'andrews'.and. i.gt.7
     +      .or. problm.eq.'nand'   .and. i.ne.5
     +      .or. problm.eq.'pump'   .and. i.eq.9
     +      .or. problm.eq.'tba'    .and.
     +              (i.ne.49+175 .and. i.ne.130+175 .and. i.ne.148+175)
     +      .or. problm.eq.'fekete' .and. i.gt.60
     +      .or. problm.eq.'plei'   .and. i.gt.14
     +      .or. problm.eq.'beam'   .and. i.gt.40


         skipime = .false.
 
    
         if (skipi) then
            skipc = 'y'
         else
            skipc = ' '
         endif


         if (skipime) then
            skipcme = 'y'
         else
            skipcme = ' '
         endif

         aerr = abs(yref(i)-y(i))
         if (.not. skipi) then
            aerrmx = max(aerrmx,aerr)
            aerrl2 = aerrl2 + aerr**2
            numa   = numa + 1
         endif
        
         
         if (abs(yref(i)).ge. min_value) then 
            rerr = abs((yref(i)-y(i))/(yref(i)))
            if (.not. skipi) then      
               rerrmx = max(rerrmx,rerr)
               rerrl2 = rerrl2 + rerr**2
               numr   = numr + 1
            endif
         endif
c      mixed relative-absolute error:
c      each code tries to put it below rtol(i)  
         mixed_rerr = abs((yref(i)-y(i))
     +                /((atol(i)/rtol(i))+abs(yref(i))))
         if (.not. skipime) then
            mixed_rerrmx = max(mixed_rerrmx,mixed_rerr)
            numme = numme+1
         end if
        
         if (printout) then 
         if (aerr.eq.0d0) then
            write(*,
     +         '('//LF// ', ''y('',i3,'') = '','//
     +         'e24.16e3,                  t72,a,t80,a)'
     +      ) i,y(i),             skipcme,skipc
         elseif (abs(yref(i)).le. min_value ) then
            write(*,
     +         '('//LF// ', ''y('',i3,'') = '','//
     +         'e24.16e3,t36,f10.2,1x,f10.2,         t72,a,t80,a)'
     +      ) i,y(i),-log10(mixed_rerr),-log10(aerr),     skipcme,
     +         skipc
         else
            write(*,
     +         '('//LF// ', ''y('',i3,'') = '','//
     +         'e24.16e3,t36,f10.2,1x,f10.2,1x,f10.2,t72,a,t80,a)'
     +     ) i,y(i),-log10(mixed_rerr),-log10(aerr),-log10(rerr),
     +       skipcme,  skipc
         endif
         end if


   10 continue

      aerrl2 = sqrt(aerrl2)
      rerrl2 = sqrt(rerrl2)

      if (aerrmx .eq.0d0)  aerrmx = uround
      if (rerrmx .eq.0d0)  rerrmx = uround
      if (mixed_rerrmx.eq.0d0) mixed_rerrmx = uround

      scda  = -log10(aerrmx)
      scdr  = -log10(rerrmx)
      mescd = -log10(mixed_rerrmx)

      if (printout) then

      write(*,*)
      write(*,'('//LF//',a,t36,i10,1x,i10,1x,i10)')
     +    'used components for  scd', numme,numa, numr
      write(*,'('//LF//',a,t36,f10.2,1x,f10.2,1x,f10.2)')
     +   'scd of Y (maximum norm)',
     +    mescd, scda, scdr
c      if (rerrmx .gt. 0d0) then
c         write(*,'('//LF//',a,t36,f10.2,1x,f10.2,1x,f10.2)')
c     +   'scd of Y (maximum norm)',
c     +   -log10(mixed_rerrmx), -log10(aerrmx), -log10(rerrmx)
c      else
c         write(*,'('//LF//',a,t36,f10.2,1x,f10.2)')
c     +   'scd of Y (maximum norm)',
c     +   -log10(mixed_rerrmx), -log10(aerrmx)
c      end if
c     write(*,'('//LF//',a,t36,f10.2,1x,f10.2)')
c    +   'scd of Y (Euclidian norm)',
c    +   -log10(aerrl2), -log10(rerrl2)
      write(*,*)
      write(*,*)

      end if
ConfigureProblem: use relative (default) or absolute error
ConfigureProblem: for the scd computation

      if (printout) then 
         write(*,'('//LF//',a,t36,f10.2)')
     +      'using mixed error yields mescd', mescd
      end if
      if (problm .eq. 'medakzo') then
         scd = scda
         if (printout) then 
         write(*,'('//LF//',a,t47,f10.2)')
     +           'using absolute error yields scd', scd
         endif
      else
c        if (rerrmx .gt. 0d0) then
         scd = scdr
         if (printout) then 
         write(*,'('//LF//',a,t58,f10.2)')
     +           'using relative error yields scd', scd
         end if
c        else
c         scd = 0
c        end if
      endif

      return
      end

      subroutine printsol(neqn,y,problm)
      integer neqn
      double precision  y(neqn)
      character*(*) problm

      integer i  

      character*3 FF, LF, PROMPT
      call getfmt(FF,LF,PROMPT)
   
c      
c  print of the computed solution
c  useful for computing the reference solution
c
      write(*,'('//LF//',a)') 'Numerical solution:'
      write(*,*)

      do i=1,neqn
        write(*,
     +         '('//LF// ', ''       y('',i3,'') = '','//
     +         'e24.16e3)'
     +     ) i,y(i)
      end do


      return
      end


      subroutine report(
     +   driver, problm, solver,
     +   rtol, atol, h0, solmax,
     +   iwork, cputim, scd,mescd
     +)
      character*(*)    driver, problm, solver
      integer          iwork(*)
      double precision rtol, atol, h0, solmax
      real             cputim
      double precision scd,mescd
c
c     print integration characteristics including CPU time
c
      integer          nsteps, naccpt, nfe, njac, nlu
      logical          batch
      external         batch

      character*(*)     FS, RS
c
c for LaTeX tabular
c
c     parameter        (FS = '''&''', RS = '''\\\\''')
c
      parameter        (FS = ''' ''', RS = '''''')

      character*3 FF, LF, PROMPT
      call getfmt(FF,LF,PROMPT)

ConfigureSolver: get integration characteristics from iwork
      if (solver.eq.'BIMD') then
         nsteps = iwork(21)
         naccpt = iwork(26)
         do i=22,25
           nsteps = nsteps + iwork(i)
           naccpt = naccpt + iwork(i+5)
         end do
         nfe    = iwork(12)
         njac   = iwork(13)
         nlu    = iwork(14)
      elseif (solver.eq.'DDASSL') then
         nsteps = iwork(11)+iwork(14)+iwork(15)
         naccpt = iwork(11)
         nfe    = iwork(12)
         njac   = iwork(13)
         nlu    = 0
      elseif (solver.eq.'DOPRI5') then
         nsteps = iwork(18)
         naccpt = iwork(19)
         nfe    = iwork(17)
         njac   = 0
         nlu    = 0
      elseif (solver.eq.'GAMD' ) then
         nsteps = 0
         do 12 i=12,23
            nsteps = nsteps + iwork(i)
   12    continue
         naccpt = iwork(12)+iwork(13)+iwork(14)+iwork(15)
         nfe    = iwork(10)
         njac   = iwork(11)
         nlu    = iwork(24)
      elseif (solver.eq.'MEBDFDAE') then
         nsteps = iwork(5)+iwork(6)
         naccpt = iwork(5)
         nfe    = iwork(7)
         njac   = iwork(8)
         nlu    = iwork(9)
      elseif (solver.eq.'MEBDFI') then
         nsteps = iwork(5)+iwork(6)
         naccpt = iwork(5)
         nfe    = iwork(7)
         njac   = iwork(8)
         nlu    = iwork(9)
      elseif (solver.eq.'PSIDE') then
         nsteps = iwork(15)
         naccpt = iwork(15)-iwork(16)-iwork(17)-iwork(18)-iwork(19)
         nfe    = iwork(11)
         njac   = iwork(12)
         nlu    = iwork(13)
      elseif (solver.eq.'PSODE') then
         nsteps = iwork(1)
         naccpt = iwork(2)
         nfe    = iwork(3)
         njac   = iwork(4)
         nlu    = iwork(5)
      elseif (solver.eq.'RADAU' .or. solver.eq.'RADAU5') then
         nsteps = iwork(16)
         naccpt = iwork(17)
         nfe    = iwork(14)
         njac   = iwork(15)
         nlu    = iwork(19)
      elseif (solver.eq.'VODE') then
         nsteps = iwork(11)+iwork(21)+iwork(22)
         naccpt = iwork(11)
         nfe    = iwork(12)
         njac   = iwork(13)
         nlu    = iwork(19)
      else
         write(*,*) 'ERROR: report.f: does not support ',
     +              solver, ' (yet)'
         stop
      endif

ConfigureSolver: not all solvers provide all integration characteristics
      write(*,*)
      write(*,'('//LF//',a)')
     +        'Integration characteristics:'
      write(*,*)
      write(*,'('//LF//',3x,a,t35,i8)')
     +        'number of integration steps', nsteps
      write(*,'('//LF//',3x,a,t35,i8)')
     +        'number of accepted steps', naccpt
      write(*,'('//LF//',3x,a,t35,i8)')
     +        'number of f evaluations', nfe
      write(*,'('//LF//',3x,a,t35,i8)')
     +        'number of Jacobian evaluations', njac
      if (solver.ne.'DASSL') then
         write(*,'('//LF//',3x,a,t35,i8)')
     +           'number of LU decompositions', nlu
      endif
      write(*,*)
      write(*,'('//LF//',a,t38,f8.4,a)')
     +        'CPU-time used:', cputim, ' sec'

      if (batch()) then
         write(*,*)
         write(*,'(a,'//
     +           '3(a,'//FS//'),'//
     +           '4(e12.3,'//FS//'),'//
     +           'f8.2,'//FS//','//
     +           '5(i8,'//FS//'),'//
     +           'f10.4,'//FS//','//
     +           'f8.2,'//FS//')')
     +      'TestRun: ',
     +      driver, problm, solver,
     +      rtol, atol, h0, solmax,
     +      scd,
     +      nsteps, naccpt, nfe, njac, nlu,
     +      cputim, 
     +      mescd
      endif

      return
      end

      logical function batch()
      batch = .false.
      return
      end

      subroutine getfmt(LF, FF, PROMPT)
      character*(*) LF, FF, PROMPT
c
c standard Fortran carriage control
c
      LF = '''1'''
      FF = ''' '''
      PROMPT = '   '
c
c if supported, you might prefer
c
c     LF = ''''''
c     FF = ''''''
c     PROMPT = ',$'
c
      return
      end

      integer function lnblnk(string)
      character*(*) string
c
c     return index of last non blank character in string
c
      integer i

      do 10 i=len(string),1,-1
         if (string(i:i) .ne. ' ') then
            lnblnk = i
            return
         endif
   10 continue

      lnblnk = 0
      return

      end


      integer function lrnblnk(string)
      character*(*) string
c
c     return index of last non blank character in string
c
      integer i

      do 10 i=1,len(string)
         if (string(i:i) .ne. ' ') then
            lrnblnk = i
            return
         endif
   10 continue

      lrnblnk = 1
      return

      end


      subroutine mfileout(problm,solver,filepath,nindsol,indsol)
c utility that generate a MATLAB and a SCILAB function to print the plots
c of the solution
      character problm*8, solver*8
      character  filemout*140,fileindout*140
      character  filepath*100, filesciout*140
      integer nindsol, indsol(*)
      external lnblnk

      integer i
     
      write(filemout,'(a,a,a,a)') filepath(1:lnblnk(filepath)),
     + problm(1:lnblnk(problm)), solver(1:lnblnk(solver)),'.m'

      write(filesciout,'(a,a,a,a)') filepath(1:lnblnk(filepath)),
     + problm(1:lnblnk(problm)), solver(1:lnblnk(solver)),'.sci'


      write(fileindout,'(a,a,a,a)') filepath(1:lnblnk(filepath)),
     +   problm(1:lnblnk(problm)), solver(1:lnblnk(solver)),'.ind'

       
      open(UNIT=95,FILE=filemout)
      write(95,*)'% This file is automatically generated by '
      write(95,*)'% the IVPtestset file report.f '
      write(95,*) 'nindsol = ',nindsol, ';'
      write(95,*)  'indsol = ['
      do i=1,nindsol 
           write(95,*) indsol(i)
      end do
      write(95,*) '];'
      write(95,*)'cd ',filepath(1:lnblnk(filepath))
      write(95,*) 'addpath(''', filepath(1:lnblnk(filepath)),''')'
      write(95,*)  'fl=IVPfplot(''', problm(1:lnblnk(problm)),
     +       ''',''',   solver(1:lnblnk(solver)),''',indsol,1,...'
      write(95,*)
     +       '''',filepath(1:lnblnk(filepath)),''');'
      close(95)

      open(unit=95,FILE=fileindout)
      write(95,'(i3)') nindsol
      close(95)

       open(UNIT=95,FILE=filesciout)
       write(95,*)'// This file is automatically generated by '
       write(95,*)'// the IVPtestset file report.f '
       write(95,*) 'nindsol = ',nindsol, ';'
       write(95,*)  'indsol = ['
       do i=1,nindsol 
           write(95,*) indsol(i)
       end do
       write(95,*) '];'
       write(95,*) 'getf(strcat([''',filepath(1:lnblnk(filepath)),''',  
     +             ''IVPfplot.sci'']))'
       write(95,*)  'fl=IVPfplot(''', problm(1:lnblnk(problm)),
     +       ''',''',   solver(1:lnblnk(solver)),''',indsol,1,...'
       write(95,*)
     +       '''',filepath(1:lnblnk(filepath)),''');'
       close(95)
       open(unit=95,FILE=fileindout)
       write(95,'(i3)') nindsol
       close(95)


      return
      end
