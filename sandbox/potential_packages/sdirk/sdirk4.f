      SUBROUTINE SDIRK4(N,FCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC ,IJAC,MLJAC,MUJAC,
     &                  MAS ,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,LRCONT,IDID)
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A STIFF 
C     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  BY'=F(X,Y).
C     TH METHOD USED IS A SINGLY DIAGONALLY IMPLICIT RUNGE-KUTTA METHOD
C     OF ORDER 4 (WITH STEP SIZE CONTROL).
C     C.F. SECTION IV.6
C
C     AUTHORS: E. HAIRER AND G. WANNER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND 
C              E-MAIL:  HAIRER@CGEUGE51.BITNET,  WANNER@CGEUGE51.BITNET
C     
C     THIS CODE IS PART OF THE BOOK:
C         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
C         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
C         SPRINGER-VERLAG (1990)               
C      
C     VERSION OF JUNE 22, 1995
C
C     INPUT PARAMETERS  
C     ----------------  
C     N           DIMENSION OF THE SYSTEM 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F)
C                    REAL*8 X,Y(N),F(N)
C                    F(1)=...   ETC.
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, 
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
C                 ADAPTS ITS STEP SIZE. STUDY THE CHOSEN VALUES FOR A FEW
C                 STEPS IN SUBROUTINE "SOLOUT", WHEN YOU ARE NOT SURE.
C                 (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
C                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
C                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
C                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(N,X,Y,DFY,LDFY)
C                    REAL*8 X,Y(N),DFY(LDFY,N)
C                    DFY(1,1)= ...
C                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS
C                 FURNISHED BY THE CALLING PROGRAM.
C                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
C                    BE FULL AND THE PARTIAL DERIVATIVES ARE
C                    STORED IN DFY AS
C                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
C                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
C                    THE PARTIAL DERIVATIVES ARE STORED
C                    DIAGONAL-WISE AS
C                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
C
C     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
C                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
C                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
C
C     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
C                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN 
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C
C     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLJAC=N.
C
C     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
C     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
C
C     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
C                 MATRIX M.
C                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
C                 MATRIX AND NEEDS NOT TO BE DEFINED;
C                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
C                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
C                    SUBROUTINE MAS(N,AM,LMAS)
C                    REAL*8 AM(LMAS,N)
C                    AM(1,1)= ....
C                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
C                    AS FULL MATRIX LIKE
C                         AM(I,J) = M(I,J)
C                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
C                    DIAGONAL-WISE AS
C                         AM(I-J+MUMAS+1,J) = M(I,J).
C
C     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
C                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
C                       MATRIX, MAS IS NEVER CALLED.
C                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
C
C     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
C                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
C
C     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLMAS=N.
C                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,IRTRN)
C                    REAL*8 X,Y(N)
C                    ....  
C                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
C                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
C                    THE FIRST GRID-POINT).
C                 "XOLD" IS THE PRECEEDING GRID-POINT.
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM.
C           
C          -----  CONTINUOUS OUTPUT: -----
C                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
C                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
C                 THE REAL*8 FUNCTION
C                        >>>   CONTS4(I,S)   <<<
C                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
C                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
C                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
C
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST
C                             N*(LJAC+LMAS+LE1+12)+7
C                 WHERE
C                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
C                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.)
C                 AND                  
C                    LMAS=0              IF IMAS=0
C                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
C                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.)
C                 AND
C                    LE1=N               IF MLJAC=N (FULL JACOBIAN)
C                    LE1=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.).
C                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
C                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
C                 STORAGE REQUIREMENT IS 
C                             LWORK = 2*N*N+12*N+7.
C
C     LWORK       DECLARED LENGHT OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGHT "LIWORK".
C                 "LIWORK" MUST BE AT LEAST 2*N+4.
C
C     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
C
C     LRCONT      DECLARED LENGTH OF COMMON BLOCK
C                  >>>  COMMON /CONT/ICONT(4),RCONT(LRCONT)  <<<
C                 WHICH MUST BE DECLARED IN THE CALLING PROGRAM.
C                 "LRCONT" MUST BE AT LEAST
C                             5*N+2 .
C                 THIS IS USED FOR STORING THE COEFFICIENTS OF THE
C                 CONTINUOUS SOLUTION AND MAKES THE CALLING LIST FOR THE
C                 FUNCTION "CONTS4" AS SIMPLE AS POSSIBLE.
C
C ----------------------------------------------------------------------
C 
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK 
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(7)
C              AS WELL AS IWORK(1),..,IWORK(4) DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
C              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
C              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
C              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N)
C              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1). IT IS
C              ALSO NOT GOOD FOR SPARSE JACOBIANS.
C
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
C
C    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE
C              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP.
C              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7.
C
C    IWORK(4)  SWITCH FOR THE COEFFICIENTS OF THE METHOD
C              IWORK(4)=1  COEFFICIENTS WITH GAMMA=0.25
C              IWORK(4)=2  COEFFICIENTS WITH GAMMA=4./15.
C              THE DEFAULT VALUE (FOR IWORK(4)=0) IS 2.
C               
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
C              DEFAULT 0.9D0.
C
C    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
C              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
C              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER 
C              (0.001D0, SAY).      
C              DEFAULT 0.001D0.
C
C    WORK(4)   STOPPING CRIERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
C              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER.
C              DEFAULT 0.03D0.
C
C    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE
C              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A
C              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR
C              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE
C              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS
C              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD.
C              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 .
C
C    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
C
C-----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS 
C     ----------------- 
C     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND).
C
C     Y(N)        NUMERICAL SOLUTION AT X
C 
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID=1  COMPUTATION SUCCESSFUL,
C                   IDID=-1 COMPUTATION UNSUCCESSFUL.
C
C-----------------------------------------------------------------------
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS 
C *** *** *** *** *** *** *** *** *** *** *** *** *** 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N),ATOL(1),RTOL(1),WORK(LWORK),IWORK(LIWORK)
      LOGICAL IMPLCT,JBAND,ARRET
      EXTERNAL FCN,JAC,MAS,SOLOUT
      COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
C --- COMMON STAT CAN BE INSPECTED FOR STATISTICAL PURPOSES:
C ---    NFCN      NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                  EVALUATION OF THE JACOBIAN ARE NOT COUNTED)  
C ---    NJAC      NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                  OR NUMERICALLY)
C ---    NSTEP     NUMBER OF COMPUTED STEPS
C ---    NACCPT    NUMBER OF ACCEPTED STEPS
C ---    NREJCT    NUMBER OF REJECTED STEPS (AFTER AT LEAST ONE STEP
C                  HAS BEEN ACCEPTED)
C ---    NDEC      NUMBER OF LU-DECOMPOSITIONS
C ---    NSOL      NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS 
C *** *** *** *** *** *** ***
       NFCN=0
       NJAC=0
       NSTEP=0
       NACCPT=0
       NREJCT=0
       NDEC=0
       NSOL=0
       ARRET=.FALSE.
C -------- SWITCH FOR TRANSFORMATION OF JACOBIAN TO HESSIAN FORM ---
      NHESS=IWORK(1)
      IF (N.LE.2) NHESS=0
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF(IWORK(2).EQ.0)THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(2)
         IF(NMAX.LE.0)THEN
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS
      IF(IWORK(3).EQ.0)THEN
         NIT=7
      ELSE
         NIT=IWORK(3)
         IF(NIT.LE.0)THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C -------- METH    SWITCH FOR THE COEFFICIENTS OF THE METHOD 
      IF(IWORK(4).EQ.0)THEN
         METH=2
      ELSE
         METH=IWORK(4)
         IF(METH.LE.0.OR.METH.GE.3)THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(4)=',IWORK(4)
            ARRET=.TRUE.
         END IF
      END IF
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  
      IF(WORK(1).EQ.0.D0)THEN
         UROUND=1.D-16
      ELSE
         UROUND=WORK(1)
         IF(UROUND.LE.1.D-19.OR.UROUND.GE.1.D0)THEN
            WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION
      IF(WORK(2).EQ.0.D0)THEN
         SAFE=0.9D0
      ELSE
         SAFE=WORK(2)
         IF(SAFE.LE..001D0.OR.SAFE.GE.1.D0)THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(2)=',WORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
      IF(WORK(3).EQ.0.D0)THEN
         THET=0.001D0
      ELSE
         THET=WORK(3)
      END IF
C --- FNEWT   STOPPING CRIERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
      IF(WORK(4).EQ.0.D0)THEN
         FNEWT=0.03D0
      ELSE
         FNEWT=WORK(4)
      END IF
C --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST.
      IF(WORK(5).EQ.0.D0)THEN
         QUOT1=1.D0
      ELSE
         QUOT1=WORK(5)
      END IF
      IF(WORK(6).EQ.0.D0)THEN
         QUOT2=1.2D0
      ELSE
         QUOT2=WORK(6)
      END IF
C -------- MAXIMAL STEP SIZE
      IF(WORK(7).EQ.0.D0)THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(7)
      END IF
C --------- CHECK IF TOLERANCES ARE O.K.
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
      ELSE
          DO 15 I=1,N
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
  15      CONTINUE
      END IF
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C         COMPUTATION OF ARRAY ENTRIES
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C ---- IMPLICIT, BANDED OR NOT ?
      IMPLCT=IMAS.NE.0
      JBAND=MLJAC.NE.N
      ARRET=.FALSE.
C -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS ---
C -- JACOBIAN 
      IF(JBAND)THEN
         LDJAC=MLJAC+MUJAC+1
      ELSE
         LDJAC=N
      END IF
C -- MATRIX E FOR LINEAR ALGEBRA
      IF(JBAND)THEN
         LDE=2*MLJAC+MUJAC+1
      ELSE
         LDE=N
      END IF
C -- MASS MATRIX
      IF (IMPLCT) THEN
          IF (MLMAS.NE.N) THEN
              LDMAS=MLMAS+MUMAS+1
          ELSE
              LDMAS=N
          END IF
C ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
              WRITE (6,*) 'BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF
     & "JAC"'
            ARRET=.TRUE.
          END IF
      ELSE
          LDMAS=0
      END IF
      LDMAS2=MAX(1,LDMAS)
C ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN
      IF ((IMPLCT.OR.JBAND).AND.NHESS.NE.0) THEN
         WRITE(6,*)' HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH 
     &FULL JACOBIAN'
         ARRET=.TRUE.
      END IF
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEYHAT=8
      IEZ=IEYHAT+N
      IEY0=IEZ+N
      IEZ1=IEY0+N
      IEZ2=IEZ1+N
      IEZ3=IEZ2+N
      IEZ4=IEZ3+N
      IEZ5=IEZ4+N
      IESCAL=IEZ5+N
      IEF1=IESCAL+N
      IEG1=IEF1+N
      IEH1=IEG1+N
      IEJAC=IEH1+N
      IEMAS=IEJAC+N*LDJAC
      IEE=IEMAS+N*LDMAS
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IEE+N*LDE-1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
      IEIP=5
      IEHES=IEIP+N
C --------- TOTAL REQUIREMENT ---------------
      ISTORE=IEHES+N-1
      IF(ISTORE.GT.LIWORK)THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C --------- CONTROL OF LENGTH OF COMMON BLOCK "CONT" -------
      IF(LRCONT.LT.(5*N+2))THEN
         WRITE(6,*)' INSUFF. STORAGE FOR RCONT, MIN. LRCONT=',5*N+2
         ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL SDICOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,METH,NHESS,
     &   IMPLCT,JBAND,LDJAC,LDE,LDMAS2,
     &   WORK(IEYHAT),WORK(IEZ),WORK(IEY0),WORK(IEZ1),WORK(IEZ2),
     &   WORK(IEZ3),WORK(IEZ4),WORK(IEZ5),WORK(IESCAL),WORK(IEF1),
     &   WORK(IEG1),WORK(IEH1),WORK(IEJAC),WORK(IEE),
     &   WORK(IEMAS),IWORK(IEIP),IWORK(IEHES))
C ----------- RETURN -----------
      RETURN
      END
C
C
C
C  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
C
      SUBROUTINE SDICOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,METH,NHESS,
     &   IMPLCT,BANDED,LDJAC,LE,LDMAS,
     &   YHAT,Z,Y0,Z1,Z2,Z3,Z4,Z5,SCAL,F1,G1,H1,FJAC,E,FMAS,IP,IPHES)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR SDIRK4
C     PARAMETERS SAME AS IN SDIRK4 WITH WORKSPACE ADDED 
C ---------------------------------------------------------- 
C         DECLARATIONS 
C ---------------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(N),YHAT(N),Z(N),Y0(N),Z1(N),Z2(N),Z3(N),Z4(N),Z5(N)
      REAL*8 SCAL(N),F1(N),G1(N),H1(N)
      REAL*8 FJAC(LDJAC,N),E(LE,N),FMAS(LDMAS,N)
      REAL*8 ATOL(1),RTOL(1)
      INTEGER IP(N),IPHES(N)
      LOGICAL REJECT,FIRST,IMPLCT,BANDED,CALJAC,NEWTRE
      COMMON /CONT/NN,NN2,NN3,NN4,XOLD,HSOL,CONT(1)
      COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
C --------- DUPLIFY N FOR COMMON BLOCK CONT -----
      NN=N
      NN2=2*N
      NN3=3*N
      NN4=4*N
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
      IF(IMPLCT)CALL MAS(N,FMAS,LDMAS)
C ---------- CONSTANTS ---------
      MBDIAG=MUMAS+1
      IF (BANDED) THEN
          MLE=MLJAC
          MUE=MUJAC
          MBJAC=MLJAC+MUJAC+1
          MBB=MLMAS+MUMAS+1
          MDIAG=MLE+MUE+1
          MDIFF=MLE+MUE-MUMAS
      END IF
      IF (METH.EQ.1) THEN
C -------- METHOD WITH GAMMA =1/4 ----------
          GAMMA=0.25D0
          C2=0.75D0
          C3=11.0D0/20.0D0
          C4=0.5D0
          ALPH21=2.0D0
          ALPH31=42.0D0/25.0D0
          ALPH32=-4.0D0/25.0D0
          ALPH41=89.0D0/68.0D0
          ALPH42=-25.0D0/136.0D0
          ALPH43=15.0D0/136.0D0
          ALPH51=-37.0D0/12.0D0
          ALPH52=-103.0D0/24.0D0
          ALPH53=275.0D0/8.0D0
          ALPH54=-85.0D0/3.0D0
          E1=-23.0D0/6.0D0
          E2=-17.0D0/12.0D0
          E3=125.0D0/4.0D0
          E4=-85.0D0/3.0D0
          D11=61.D0/27.D0
          D12=-185.D0/54.D0
          D13=2525.D0/18.D0
          D14=-3740.D0/27.D0
          D15=-44.D0/9.D0
          D21=2315.D0/81.D0
          D22=1049.D0/162.D0
          D23=-27725.D0/54.D0
          D24=40460.D0/81.D0
          D25=557.D0/27.D0
          D31=-6178.D0/81.D0
          D32=-1607.D0/81.D0
          D33=20075.D0/27.D0
          D34=-56440.D0/81.D0
          D35=-718.D0/27.D0
          D41=3680.D0/81.D0
          D42=1360.D0/81.D0
          D43=-10000.D0/27.D0
          D44=27200.D0/81.D0
          D45=320.D0/27.D0
          ETA1= 3.D0
          ANU1= 0.88D0
          ANU2= 0.44D0
          AMU1= 3.D0/17.D0
          AMU3= 155.D0/187.D0
      END IF
      IF (METH.EQ.2) THEN
C ---------- METHOD WITH GAMMA = 4/15 ---------------
          GAMMA=4.0D0/15.0D0
          C2=23.0D0/30.0D0
          C3=17.0D0/30.0D0
          C4=2881.0D0/28965.0D0+GAMMA
          ALPH21=15.0D0/8.0D0
          ALPH31=1577061.0D0/922880.0D0
          ALPH32=-23427.0D0/115360.0D0
          ALPH41=647163682356923881.0D0/2414496535205978880.0D0
          ALPH42=-593512117011179.0D0/3245291041943520.0D0
          ALPH43=559907973726451.0D0/1886325418129671.0D0
          ALPH51=724545451.0D0/796538880.0D0
          ALPH52=-830832077.0D0/267298560.0D0
          ALPH53=30957577.0D0/2509272.0D0
          ALPH54=-69863904375173.0D0/6212571137048.0D0
          E1=7752107607.0D0/11393456128.0D0
          E2=-17881415427.0D0/11470078208.0D0
          E3=2433277665.0D0/179459416.0D0
          E4=-96203066666797.0D0/6212571137048.0D0
          D11= 24.74416644927758D0
          D12= -4.325375951824688D0
          D13= 41.39683763286316D0
          D14= -61.04144619901784D0
          D15= -3.391332232917013D0
          D21= -51.98245719616925D0
          D22= 10.52501981094525D0
          D23= -154.2067922191855D0
          D24= 214.3082125319825D0
          D25= 14.71166018088679D0
          D31= 33.14347947522142D0
          D32= -19.72986789558523D0
          D33= 230.4878502285804D0
          D34= -287.6629744338197D0
          D35= -18.99932366302254D0
          D41= -5.905188728329743D0
          D42= 13.53022403646467D0
          D43= -117.6778956422581D0
          D44= 134.3962081008550D0
          D45= 8.678995715052762D0
         ETA1=23.D0/8.D0
         ANU1= 0.9838473040915402D0
         ANU2= 0.3969226768377252D0
         AMU1= 0.6563374010466914D0
         AMU3= 0.3372498196189311D0
      END IF
      POSNEG=SIGN(1.D0,XEND-X)
      HMAX1=MIN(ABS(HMAX),ABS(XEND-X))
      IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6
      H=MIN(ABS(H),HMAX1) 
      H=SIGN(H,POSNEG) 
      HOLD=H
      CFAC=SAFE*(1+2*NIT)
      NEWTRE=.FALSE.
      REJECT=.FALSE.
      FIRST=.TRUE.
      FACCO1=1.D0
      FACCO2=1.D0
      FACCO3=1.D0
      FACCO4=1.D0
      FACCO5=1.D0
      NSING=0
      XOLD=X
      IF (IOUT.NE.0) THEN
          IRTRN=1
          NRSOL=1
          XOSOL=XOLD
          XSOL=X
          DO 7 I=1,N
  7       CONT(I)=Y(I)
          NSOLU=N
          HSOL=HOLD
          CALL SOLOUT(NRSOL,XOSOL,XSOL,CONT,NSOLU,IRTRN)
          IF (IRTRN.LT.0) GOTO 79
      END IF
      IF (ITOL.EQ.0) THEN
          DO 8 I=1,N
   8      SCAL(I)=ATOL(1)+RTOL(1)*DABS(Y(I))
      ELSE
          DO 9 I=1,N
   9      SCAL(I)=ATOL(I)+RTOL(I)*DABS(Y(I))
      END IF
      IF (IJAC.EQ.0) CALL FCN(N,X,Y,Y0)
      IF (IJAC.EQ.0) NFCN=NFCN+1
C --- BASIC INTEGRATION STEP  
  10  CONTINUE
C *** *** *** *** *** *** ***
C  COMPUTATION OF THE JACOBIAN
C *** *** *** *** *** *** ***
      NJAC=NJAC+1
      IF (IJAC.EQ.0) THEN
C --- COMPUTE JACOBIAN MATRIX NUMERICALLY
          IF (BANDED) THEN
C --- JACOBIAN IS BANDED
              MUJACP=MUJAC+1
              MD=MIN(MBJAC,N)
              DO 16 K=1,MD
              J=K
 12           Z(J)=Y(J)
              G1(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J))))
              Y(J)=Y(J)+G1(J)
              J=J+MD
              IF (J.LE.N) GOTO 12 
              CALL FCN(N,X,Y,CONT)
              J=K
              LBEG=MAX(1,J-MUJAC)
 14           LEND=MIN(N,J+MLJAC)
              Y(J)=Z(J)
              MUJACJ=MUJACP-J
              DO 15 L=LBEG,LEND
 15           FJAC(L+MUJACJ,J)=(CONT(L)-Y0(L))/G1(J) 
              J=J+MD
              LBEG=LEND+1
              IF (J.LE.N) GOTO 14
 16           CONTINUE
          ELSE
C --- JACOBIAN IS FULL
              DO 18 I=1,N
              YSAFE=Y(I)
              DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
              Y(I)=YSAFE+DELT
              CALL FCN(N,X,Y,CONT)
              DO 17 J=1,N
  17          FJAC(J,I)=(CONT(J)-Y0(J))/DELT
  18          Y(I)=YSAFE
          END IF
      ELSE
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
          CALL JAC(N,X,Y,FJAC,LDJAC)
      END IF
      CALJAC=.TRUE.
      IF (NHESS.NE.0) CALL ELMHES (LDJAC,N,1,N,FJAC,IPHES) 
  20  CONTINUE
C *** *** *** *** *** *** ***
C  COMPUTE THE MATRIX E AND ITS DECOMPOSITION
C *** *** *** *** *** *** ***
      FAC1=1.D0/(H*GAMMA)
      IF (IMPLCT) THEN
          IF (BANDED) THEN
C --- THE MATRIX E (MAS IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX)
              DO 127 J=1,N
              I1J=MAX0(1,MUJAC+2-J)
              I2J=MIN(MBJAC,MUJAC+1-J+N)
              DO 125 I=I1J,I2J
 125          E(I+MLE,J)=-FJAC(I,J)
              I1B=MAX0(1,MUMAS+2-J)
              I2B=MIN0(MBB,MUMAS+1-J+N)
              DO 126 I=I1B,I2B
              IB=I+MDIFF
 126          E(IB,J)=E(IB,J)+FAC1*FMAS(I,J)
 127          CONTINUE
              CALL DECB(N,LE,E,MLE,MUE,IP,IER)
              IF (IER.NE.0) GOTO 78
          ELSE
              IF (MLMAS.NE.N) THEN
C --- THE MATRIX E (MAS IS A BANDED MATRIX, JACOBIAN A FULL MATRIX)
                  DO 225 J=1,N
                  DO 225 I=1,N
 225              E(I,J)=-FJAC(I,J)
                  DO 226 J=1,N
                  I1=MAX0(1,J-MUMAS)
                  I2=MIN0(N,J+MLMAS)
                  DO 226 I=I1,I2
 226              E(I,J)=E(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
                  CALL DEC(N,LE,E,IP,IER)
                  IF (IER.NE.0) GOTO 78
              ELSE
C --- THE MATRIX E (MAS IS A FULL MATRIX, JACOBIAN A FULL MATRIX)
                  IF (MLJAC.EQ.N) THEN
                      DO 324 J=1,N
                      DO 324 I=1,N
 324                  E(I,J)=FMAS(I,J)*FAC1-FJAC(I,J)
                      CALL DEC(N,LE,E,IP,IER)
                      IF (IER.NE.0) GOTO 78
                  ELSE
C --- THE MATRIX E (B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX)
                  END IF
              END IF
          END IF
      ELSE
          IF (BANDED) THEN
C --- THE MATRIX E (MAS=IDENTITY, JACOBIAN A BANDED MATRIX)
              DO 427 J=1,N
              I1J=MAX0(1,MUJAC+2-J)
              I2J=MIN(MBJAC,MUJAC+1-J+N)
              DO 425 I=I1J,I2J
 425          E(I+MLE,J)=-FJAC(I,J)
 427          E(MDIAG,J)=E(MDIAG,J)+FAC1
              CALL DECB(N,LE,E,MLE,MUE,IP,IER)
              IF (IER.NE.0) GOTO 78
          ELSE
C --- THE MATRIX E (MAS=IDENTITY, JACOBIAN A FULL MATRIX)
              IF (NHESS.EQ.0) THEN
                  DO 526 J=1,N
                  DO 525 I=1,N
 525              E(I,J)=-FJAC(I,J)
 526              E(J,J)=E(J,J)+FAC1
                  CALL DEC(N,LE,E,IP,IER)
                  IF (IER.NE.0) GOTO 79 
              ELSE
                  DO 624 J=1,N-1
                  J1=J+1
 624              E(J1,J)=-FJAC(J1,J)
                  DO 626 J=1,N
                  DO 625 I=1,J
 625              E(I,J)=-FJAC(I,J)
 626              E(J,J)=E(J,J)+FAC1
                  CALL DECH(N,LE,E,1,IP,IER)
                  IF (IER.NE.0) GOTO 79 
              END IF
          END IF
      END IF
      NDEC=NDEC+1
  30  CONTINUE
      IF (NSTEP.GT.NMAX.OR.X+.1D0*H.EQ.X.OR.ABS(H).LE.UROUND) GOTO 79
      XPH=X+H
C --- LOOP FOR THE 5 STAGES
      FACCO1=DMAX1(FACCO1,UROUND)**0.8D0
      FACCO2=DMAX1(FACCO2,UROUND)**0.8D0
      FACCO3=DMAX1(FACCO3,UROUND)**0.8D0
      FACCO4=DMAX1(FACCO4,UROUND)**0.8D0
      FACCO5=DMAX1(FACCO5,UROUND)**0.8D0
C *** *** *** *** *** *** ***
C  STARTING VALUES FOR NEWTON ITERATION
C *** *** *** *** *** *** ***
      DO 59 ISTAGE=1,5
      IF (ISTAGE.EQ.1) THEN
          XCH=X+GAMMA*H
          IF (FIRST.OR.NEWTRE) THEN
              DO 132 I=1,N
 132          Z(I)=0.D0
          ELSE
              S=1.D0+GAMMA*H/HOLD
              DO 232 I=1,N
 232          Z(I)=S*(CONT(I+NN)+S*(CONT(I+NN2)+S*(CONT(I+NN3)
     &             +S*CONT(I+NN4))))-YHAT(I)
          END IF
          DO 31 I=1,N
  31      G1(I)=0.D0 
          FACCON=FACCO1
      END IF
      IF (ISTAGE.EQ.2) THEN
          XCH=X+C2*H
          DO 131 I=1,N
          Z1I=Z1(I)
          Z(I)=ETA1*Z1I
 131      G1(I)=ALPH21*Z1I
          FACCON=FACCO2
      END IF
      IF (ISTAGE.EQ.3) THEN
          XCH=X+C3*H
          DO 231 I=1,N
          Z1I=Z1(I)
          Z2I=Z2(I)
          Z(I)=ANU1*Z1I+ANU2*Z2I
 231      G1(I)=ALPH31*Z1I+ALPH32*Z2I
          FACCON=FACCO3
      END IF
      IF (ISTAGE.EQ.4) THEN
          XCH=X+C4*H
          DO 331 I=1,N
          Z1I=Z1(I)
          Z3I=Z3(I)
          Z(I)=AMU1*Z1I+AMU3*Z3I
 331      G1(I)=ALPH41*Z1I+ALPH42*Z2(I)+ALPH43*Z3I
          FACCON=FACCO4
      END IF
      IF (ISTAGE.EQ.5) THEN
          XCH=XPH
          DO 431 I=1,N
          Z1I=Z1(I)
          Z2I=Z2(I)
          Z3I=Z3(I)
          Z4I=Z4(I)
          Z(I)=E1*Z1I+E2*Z2I+E3*Z3I+E4*Z4I
          YHAT(I)=Z(I)
 431      G1(I)=ALPH51*Z1I+ALPH52*Z2I+ALPH53*Z3I+ALPH54*Z4I
          FACCON=FACCO5
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
            NEWT=0
            THETA=ABS(THET)
            IF (REJECT) THETA=2*ABS(THET)
  40        CONTINUE
            IF (NEWT.GE.NIT) THEN
                H=H/2.D0
                REJECT=.TRUE.
                NEWTRE=.TRUE.
                IF (CALJAC) GOTO 20
                GOTO 10
            END IF
C ---     COMPUTE THE RIGHT-HAND SIDE
            DO 41 I=1,N
            H1(I)=G1(I)-Z(I)
  41        CONT(I)=Y(I)+Z(I)
            CALL FCN(N,XCH,CONT,F1)
            NFCN=NFCN+1
C ---     SOLVE THE LINEAR SYSTEMS
            IF (IMPLCT) THEN
                IF (MLMAS.NE.N) THEN
                    DO 146 I=1,N
                    S1=0.0D0
                    J1B=MAX0(1,I-MLMAS)
                    J2B=MIN0(N,I+MUMAS)
                    DO 145 J=J1B,J2B
 145                S1=S1+FMAS(I-J+MBDIAG,J)*H1(J)
 146                F1(I)=S1*FAC1+F1(I)
                    IF (BANDED) THEN
                        CALL SOLB(N,LE,E,MLE,MUE,F1,IP)
                    ELSE
                        CALL SOL(N,LE,E,F1,IP)
                    END IF
                ELSE
                    DO 246 I=1,N
                    S1=0.0D0
                    DO 245 J=1,N
 245                S1=S1+FMAS(I,J)*H1(J)
 246                F1(I)=S1*FAC1+F1(I)
                    CALL SOL(N,LE,E,F1,IP)
                END IF
            ELSE
                DO 345 I=1,N
 345            F1(I)=H1(I)*FAC1+F1(I)
                IF (BANDED) THEN
                    CALL SOLB(N,LE,E,MLE,MUE,F1,IP)
                ELSE
                    IF (NHESS.EQ.0) THEN
                        CALL SOL(N,LE,E,F1,IP)
                    ELSE
                        DO 140 MMM=N-2,1,-1
                        MP=N-MMM
                        MP1=MP-1
                        I=IPHES(MP)
                        IF (I.EQ.MP) GOTO 110
                        ZSAFE=F1(MP)
                        F1(MP)=F1(I)
                        F1(I)=ZSAFE
 110                    CONTINUE
                        DO 100 I=MP+1,N 
 100                    F1(I)=F1(I)-FJAC(I,MP1)*F1(MP)
 140                    CONTINUE
                            CALL SOLH(N,LE,E,1,F1,IP)
                        DO 240 MMM=1,N-2
                        MP=N-MMM
                        MP1=MP-1
                        DO 200 I=MP+1,N 
 200                    F1(I)=F1(I)+FJAC(I,MP1)*F1(MP)
                        I=IPHES(MP)
                        IF (I.EQ.MP) GOTO 240
                        ZSAFE=F1(MP)
                        F1(MP)=F1(I)
                        F1(I)=ZSAFE
 240                    CONTINUE
                    END IF
                END IF
            END IF
            NEWT=NEWT+1
            DYNO=0.D0
            DO 57 I=1,N
            DENOM=SCAL(I)
  57        DYNO=DYNO+(F1(I)/DENOM)**2
            DYNO=DSQRT(DYNO/N)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
            IF (NEWT.GE.2.AND.NEWT.LT.NIT) THEN
                THETA=DYNO/DYNOLD
                IF (THETA.LT.0.99D0) THEN
                    FACCON=THETA/(1.0D0-THETA)
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)
                    QNEWT=DMAX1(1.0D-4,DMIN1(16.0D0,DYTH/FNEWT))
                    IF (QNEWT.GE.1.0D0) THEN
                         H=.8D0*H*QNEWT**(-1.0D0/(NIT-NEWT))
                         REJECT=.TRUE.
                         NEWTRE=.TRUE.
                         IF (CALJAC) GOTO 20
                         GOTO 10
                    END IF
                ELSE
                    NEWTRE=.TRUE.
                    GOTO 78
                END IF
            END IF
            DYNOLD=DYNO
            DO 58 I=1,N
  58        Z(I)=Z(I)+F1(I)
            NSOL=NSOL+1
            IF (FACCON*DYNO.GT.FNEWT) GOTO 40
C --- END OF SIMPILFIED NEWTON
      IF (ISTAGE.EQ.1) THEN
          DO 158 I=1,N
 158      Z1(I)=Z(I)
          FACCO1=FACCON
      END IF
      IF (ISTAGE.EQ.2) THEN
          DO 258 I=1,N
 258      Z2(I)=Z(I)
          FACCO2=FACCON
      END IF
      IF (ISTAGE.EQ.3) THEN
          DO 358 I=1,N
 358      Z3(I)=Z(I)
          FACCO3=FACCON
      END IF
      IF (ISTAGE.EQ.4) THEN
          DO 458 I=1,N
 458      Z4(I)=Z(I)
          FACCO4=FACCON
      END IF
      IF (ISTAGE.EQ.5) THEN
          DO 558 I=1,N
 558      Z5(I)=Z(I)
          FACCO5=FACCON
      END IF
  59  CONTINUE
C *** *** *** *** *** *** ***
C  ERROR ESTIMATION  
C *** *** *** *** *** *** ***
      NSTEP=NSTEP+1
      IF (IMPLCT) THEN
          DO 359 I=1,N
 359      H1(I)=FAC1*(Z5(I)-YHAT(I))
          IF (MLMAS.EQ.N) THEN
              DO 361 I=1,N
              SUM=0.D0
              DO 360 J=1,N
 360          SUM=SUM+FMAS(I,J)*H1(J)
 361          CONT(I)=SUM
          ELSE
              DO 363 I=1,N
              SUM=0.D0
              J1B=MAX0(1,I-MLMAS)
              J2B=MIN0(N,I+MUMAS)
              DO 362 J=J1B,J2B
 362          SUM=SUM+FMAS(I-J+MBDIAG,J)*H1(J)
 363          CONT(I)=SUM
          END IF
      ELSE
          DO 461 I=1,N 
 461      CONT(I)=FAC1*(Z5(I)-YHAT(I))
      END IF
      IF (BANDED) THEN
          CALL SOLB(N,LE,E,MLE,MUE,CONT,IP)
      ELSE
               IF (NHESS.EQ.0) THEN
                   CALL SOL(N,LE,E,CONT,IP)
                ELSE
                   DO 540 MMM=N-2,1,-1
                   MP=N-MMM
                   MP1=MP-1
                   I=IPHES(MP)
                   IF (I.EQ.MP) GOTO 510
                   ZSAFE=CONT(MP)
                   CONT(MP)=CONT(I)
                   CONT(I)=ZSAFE
 510               CONTINUE
                   DO 500 I=MP+1,N 
 500               CONT(I)=CONT(I)-FJAC(I,MP1)*CONT(MP)
 540               CONTINUE
                       CALL SOLH(N,LE,E,1,CONT,IP)
                   DO 640 MMM=1,N-2
                   MP=N-MMM
                   MP1=MP-1
                   DO 600 I=MP+1,N 
 600               CONT(I)=CONT(I)+FJAC(I,MP1)*CONT(MP)
                   I=IPHES(MP)
                   IF (I.EQ.MP) GOTO 640
                   ZSAFE=CONT(MP)
                   CONT(MP)=CONT(I)
                   CONT(I)=ZSAFE
 640               CONTINUE
               END IF
      END IF
      ERR=0.D0
      DO 64 I=1,N
  64  ERR=ERR+(CONT(I)/SCAL(I))**2
      ERR=DMAX1(DSQRT(ERR/N),1.D-10)
C --- COMPUTATION OF HNEW
C --- WE REQUIRE .25<=HNEW/H<=4.
      FAC=DMIN1(SAFE,CFAC/(NEWT+2*NIT))
      QUOT=DMAX1(.25D0,DMIN1(4.D0,(ERR)**.25D0/FAC))
      HNEW=H/QUOT
C *** *** *** *** *** *** ***
C  IS THE ERROR SMALL ENOUGH ?
C *** *** *** *** *** *** ***
      IF (ERR.LT.1.D0) THEN
C --- STEP IS ACCEPTED  
         FIRST=.FALSE.
         NACCPT=NACCPT+1
         HOLD=H
         XOLD=X
C --- COEFFICIENTS FOR CONTINUOUS SOLUTION
         DO 74 I=1,N 
          Z1I=Z1(I)
          Z2I=Z2(I)
          Z3I=Z3(I)
          Z4I=Z4(I)
          Z5I=Z5(I)
         CONT(I)=Y(I)
         Y(I)=Y(I)+Z5I  
         CONT(I+NN)=D11*Z1I+D12*Z2I+D13*Z3I+D14*Z4I+D15*Z5I
         CONT(I+NN2)=D21*Z1I+D22*Z2I+D23*Z3I+D24*Z4I+D25*Z5I
         CONT(I+NN3)=D31*Z1I+D32*Z2I+D33*Z3I+D34*Z4I+D35*Z5I
         CONT(I+NN4)=D41*Z1I+D42*Z2I+D43*Z3I+D44*Z4I+D45*Z5I
         YHAT(I)=Z5I
         IF (ITOL.EQ.0) THEN
           SCAL(I)=ATOL(1)+RTOL(1)*DABS(Y(I))
         ELSE
           SCAL(I)=ATOL(I)+RTOL(I)*DABS(Y(I))
         END IF
  74     CONTINUE
         X=XPH 
         IF (IOUT.NE.0) THEN
             NRSOL=NACCPT+1
             XSOL=X
             XOSOL=XOLD
             DO 77 I=1,N
  77         H1(I)=Y(I)
             NSOLU=N
             HSOL=HOLD
             CALL SOLOUT(NRSOL,XOSOL,XSOL,H1,NSOLU,IRTRN)
             IF (IRTRN.LT.0) GOTO 79
         END IF
         CALJAC=.FALSE.
         IF ((X-XEND)*POSNEG+UROUND.GT.0.D0) THEN
            H=HOPT
            IDID=1
            RETURN
         END IF
         IF (IJAC.EQ.0) CALL FCN(N,X,Y,Y0)
         NFCN=NFCN+1
         HNEW=POSNEG*DMIN1(DABS(HNEW),HMAX1)
         HOPT=HNEW
         IF (REJECT) HNEW=POSNEG*DMIN1(DABS(HNEW),DABS(H)) 
         REJECT=.FALSE.
         NEWTRE=.FALSE.
         IF ((X+HNEW/QUOT1-XEND)*POSNEG.GT.0.D0) THEN
            H=XEND-X
         ELSE
            QT=HNEW/H
            IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) GOTO 30
            H=HNEW 
         END IF
         IF (THETA.LE.THET) GOTO 20
         GOTO 10
      ELSE
C --- STEP IS REJECTED  
         REJECT=.TRUE.
         IF (FIRST) THEN
             H=H/10.D0
         ELSE
             H=HNEW
         END IF
         IF (NACCPT.GE.1) NREJCT=NREJCT+1
         IF (CALJAC) GOTO 20
         GOTO 10
      END IF
C --- UNEXPECTED STEP-REJECTION
  78  CONTINUE
      IF (IER.NE.0) THEN
          WRITE (6,*) ' MATRIX IS SINGULAR, IER=',IER,' X=',X,' H=',H
          NSING=NSING+1
          IF (NSING.GE.6) GOTO 79
      END IF
      H=H*0.5D0
      REJECT=.TRUE.
      IF (CALJAC) GOTO 20
      GOTO 10
C --- FAIL EXIT
  79  WRITE (6,979) X,H,IER
 979  FORMAT(' EXIT OF SDIRK4 AT X=',D14.7,'   H=',D14.7,'   IER=',I4)
      IDID=-1
      RETURN
      END
C
      REAL*8 FUNCTION CONTS4(I,X) 
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT FOR THE
C     SUBROUTINE SDIRK4. IT PROVIDES AN APPROXIMATION
C     TO THE I-TH COMPONENT OF THE SOLUTION AT X.
C ----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CONT/NN,NN2,NN3,NN4,XOLD,HSOL,CONT(1)
      S=(X-XOLD)/HSOL
      CONTS4=CONT(I)+S*(CONT(I+NN)+S*(CONT(I+NN2)+S*(CONT(I+NN3)
     &     +S*CONT(I+NN4))))
      RETURN
      END
