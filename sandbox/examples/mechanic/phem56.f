C**********************************************************
      SUBROUTINE PHEM56(NQ,NV,NU,NL,FPROB, 
     &                T,Q,V,U,A,RLAM,TEND,H,
     &                RTOL,ATOL,ITOL,
     &                SOLOUT,IOUT,
     &                WK,LWK,IWK,LIWK,IDID)
C-----------------------------------------------------------------------
C    NUMERICAL SOLUTION OF A CONSTRAINED MECHANICAL SYSTEM                
C                                                                        
C                 q' = T(q,t)v                                         
C           M(t,q)v' = f(t,q,v,u) - L(q,v,u,t)*lamda                    
C                 0  = H(t,q)v + k(t,q) 
C                 u' = d(q,v,u,lambda,t)  
C 
C                      VERSION OF NOVEMBER 1995
C
C    THE LOCAL ERROR ESTIMATION AND THE STEP SIZE CONTROL IS BASED ON
C    EMBEDDED FORMULAS OF ORDERS 5 AND 4. THIS METHOD IS PROVIDED WITH
C    A DENSE OUTPUT FORMULA OF ORDER 4.  
C
C               VERSION WITH SPARSE LINEAR ALGEBRA OPTION (MA28)                            
C                                                                        
C
C     THIS CODE IS BASED ON THE CODE HEM5 OF VALERY BRASEY, WHICH
C     IS DESCRIBED IN 
C         V. BRASEY , HALF-EXPLICIT METHOD FOR SEMI-EXPLICIT 
C         DIFFERENTIAL-ALGEBRAIC EQUATIONS OF INDEX 2. THESIS, 1994 
C         UNIV. DE GENEVE, SECT. DE MATHEMATIQUES   
C
C
C     AUTHOR OF PHEM56:
C              ANDER MURUA
C              INFORMATIKA FAKULTATEA, KONPUTAZIO ETA A.A. SAILA
C              EUSKAL HERRIKO UNIBERTSITATEA
C              DONOSTIA/SAN SEBASTIAN, SPAIN
C              E-MAIL: ander@si.ehu.es 
C                                                        
C     THE PHEM56 CODE IMPLEMENTS A 6-STAGE PARTITIONED HALF-EXPLICIT 
C     RUNGE-KUTTA METHOD OF ORDER 5. THIS METHOD IS DESCRIBED IN
C       A. MURUA: PARTITIONED HALF-EXPLICIT RUNGE-KUTTA METHODS FOR
C       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX 2.  
C       SUBMITTED TO COMPUTING (1995). 
C
C     TO GET THE LATEST VERSION OF PHEM56, CONTACT WITH THE AUTHOR OR WITH
C              ERNST HAIRER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND 
C              E-MAIL:  hairer@divsun.unige.ch 
C     
C     THE USER INTERFACE OF PHEM56 IS VERY SIMILAR TO THE INTERFACE 
C     OF HEM5. AN USER'S GUIDE OF PHEM56 (OBTAINED MODIFYING THE GUIDE OF
C     HEM5) IS PROVIDED WITH THE CODE.
C     
C     
C     MANY THANKS TO E. HAIRER AND G. WANNER FOR THEIR SUPPORT AND 
C     ENCOURAGEMENT. THIS WORK HAS BEEN SUPPORTED BY THE SWISS NATIONAL
C     SCIENCE FOUNDATION.
C
C          
C     DIMENSIONS
C     ----------
C 
C    NQ (SIZE OF POSITION VECTOR)
C    NV (SIZE OF VELOCITY VECTOR, NQ>=NV)
C    NL (SIZE OF LAGRANGE MULTIPLIER VECTOR) 
C    NU (SIZE OF EXTERNAL DYNAMIC VECTOR) 
C    IWK(1) = LDG (SPARSE ALGEBRA: NON ZERO ENTRIES)
C    IWK(2) = LDF 
C
C 
C     OUTPUT
C     ------
C
C    THE PARAMETER 'SOLOUT' MUST CONTAIN THE NAME OF THE USER SUPLIED
C    ROUTINE TO (OPTIONALLY) APPLY AT THE END OF EACH STEP.
C    IF IOUT=2, THE USER SUPLIED SUBROUTINE SOLOUT IS CALLED AFTER
C    EACH SUCCESFUL STEP. DENSE OUTPUT CAN BE OBTAINED CALLING
C    TO THE FUCTION POL4(I,FIRST,NQ,NV,NU,LRDO,X,DOWK) FROM SOLOUT.
C    IF IOUT=3, ONLY THE CURRENT VALUES OF Q,V,U,A,RLAM CAN BE USED.
C
C    CALLING SEQUENCE OF SOLOUT:
C       SOLOUT (NSTEP,NQ,NV,NU,NL,LRDO,Q,V,U,A,RLAM,DOWK)
C
C
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK 
C              WELL. THEY MAY BE DEFINED BY SETTING WK(1),..,WK(8)
C              AS WELL AS IWK(11) .. IWK(16) DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWK(11)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWK(11)=0) IS 100000.
C
C    IWK(12)  SWITCH FOR A PROJECTION TO ENSURE CONSISTENT INITIAL VALUES (Q,V).
C             FOR IWK(12)=1 AN INITIAL PROJECTION FOR (Q,V) IS PERFORMED.
C              NO PROJECTION FOR (Q,V) IS DONE IF IWK(12)=0.
C              THE DEFAULT VALUE FOR IWK(12) IS 0. 
C
C    IWK(13)  FOR IWK(13).GT.0 IT IS THE NUMBER OF STEPS BETWEEN 
C              TWO PROJECTIONS ON THE MANIFOLD  DEFINED BY 0 = g(q,t).
C              FOR IWK(13).LE.0 NO PROECTION IS PERFORMED.
C              THE DEFAULT VALUE FOR IWK(13) IS 0.
C
C    IWK(14)    MODE (=0: FULL LINEAR ALGEBRA WITH DEC, FL=GQ^T 
C                     =1: IDEM WITH GENERAL FL,
C                     =2: FULL LINEAR ALGEBRA WITH DGETRF, FL=GQ^T 
C                     =3: IDEM FOR GENERAL FL
C                     =4: SPARSE FL=GQ^T, 
C                     =5: SPARSE WITH GENERAL FL)
C
C    IWK(15)  IACC ( =1: PERFORM A PROJECTION AT EACH STEP TO ACCURATELY CALCULATE 
C                       THE ACCELERATION A AND THE LAGRANGE MULTIPLIERS RLAM,
C                    =0: DO NOT PERFORM PROJECTION FOR A AND RLAM AT EACH STEP.
C                        (ONLY ONE INITIAL PROJECTION IS PERFORMEND TO TO OBTAIN 
C                         CONSISTENT INITIAL VALUES FOR A AND RLAM). 
C                        IN THAT CASE, APPROXIMATIONS TO A AND RLAM WITH O(H^4)
C                        GLOBAL ERROR ARE COMPUTED AT EACH STEP.
C
C 
C
C
C    IWK(21->29)  IPAR
C    IPAR(1) = IWK(21) = NMRC (SIZE OF A BLOCK OF AM)
C    IPAR(2) = IWK(22) = NBLK (NUMBER OF BLOCKS OF AM)
C    IPAR(3) = IWK(23) = NPGP (0 IF GP AS THE SAME PATTERN AS PREVIOUS CALL)
C    IPAR(4) = IWK(24) = NPFL (0 IF FL AS THE SAME PATTERN AS PREVIOUS CALL)
C    IPAR(5) = IWK(25) = IS (SIZE OF INTEGER WORK SPACE FOR MA28 (MIN 13*NM))
C    IPAR(6) = IWK(26) = IXS (SIZE OF REAL WORK SPACE FOR MA28 (MIN NM+4*NZA))
C    IPAR(7) = IWK(27) = PREVL
C    IPAR(8) = IWK(28) = IO
C    IPAR(9) = FLAG TO INDICATE IF UMDFAC HAS BEEN CALLED AT LEAST ONCE
C
C    IWK(31->38) ISTAT
C    ISTAT(1) = IWK(31) = NSTEP
C    ISTAT(2) = IWK(32) = NACCPT
C    ISTAT(3) = IWK(33) = NREJCT
C    ISTAT(4) = IWK(34) = NFCN
C    ISTAT(5) = IWK(35) = NGQCN
C    ISTAT(6) = IWK(36) = NAMAT
C    ISTAT(7) = IWK(37) = NDEC
C    ISTAT(8) = IWK(38) = NSOL
C ---    NFCN      NUMBER OF f - EVALUATIONS
C ---    NGCN      NUMBER OF g - EVALUATIONS
C ---    NSTEP     NUMBER OF ALL COMPUTED STEPS
C ---    NACCPT    NUMBER OF ACCEPTED STEPS
C ---    NREJCT    NUMBER OF REJECTED STEPS (AFTER AT LEAST ONE STEP
C                  HAS BEEN ACCEPTED)
C ---    NDEC      NUMBER OF LU-DECOMPOSITIONS
C ---    NSOL      NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
C    IWK(41->41+NM-1) = IP
C    IWK(29->29+NM+2*LDG-1) = INDG(1..2*LDG)
C    IWK(29+NM+2*LDG..) = INDG1(1..2*LDG)
C    IWK(29+NM+4*LDG..) = INDG2(1..2*LDG)
C    IWK(29+NM+6*LDG..) = INDGD(1..2*LDG)
C    IWK(29+NM+8*LDG..) = INDFL(1..2*LDF)
C    IWK(29+NM+4*LDG+2*LDF..) = INDAM(1..2*NZA)
C    IUMF(1) = IWK(29+NM+4*LDG+2*LDF+2*NZA..) = IKEEP(1..5*NM)
C    IUMF(IIK) = IWK(29+NM+4*LDG+2*LDF+4*NZA+5*NM..) = IW(1..8*NM)
C
C    WK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
C              DEFAULT 0.9D0.
C
C    WK(3), WK(4)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C                 WK(3) <= HNEW/HOLD <= WK(4).
C              DEFAULT VALUES: WK(3)=0.2D0, WK(4)=10.D0
C
C    WK(6)   MAXIMAL STEP SIZE, DEFAULT TEND-T.
C
C    WK(7) = BETA, DEFAULT 0.D0
C
C    WK(8) = ALPHA, DEFAULT 1/5
C
C----------------------------------------------------------------------
C    XLA(1..NZA) = AVALUE(1..4*NZA)
C    XUMF(IXW) = XLA(1+2*NZA..) = W(NM)
C    LXLA = 4*NZA +NM
C
C-----------------------------------------------------------------------
C     FPROB
C     -----
C   If IFCN =
C
C     0 -> UDOT
C     1 -> f,M,QDOT, and GQ (if FL=GQ^T) or  FL (else)
C     2 -> QDOT
C     3 -> gt
C     4 -> g
C     5 -> f and GQQ,  If GQQ is not avalaible, put IFCN=0 on output 
C                    (this indicates that GQQ must be calculated numerically) 
C     6 -> GQ,gt
C     7 -> f,M, and GQQ,  If GQQ is not avalaible, put IFCN=0 on output 
C                    (this indicates that GQQ must be calculated numerically) 
C     8 -> f,M,(FL)
C     9 -> M
C     10 -> gt,GQ,M,QDOT,(FL)
C     11 -> gt,GQ,M,(FL)
C       
C--------------------------------------------------------------------------
C     OUTPUT PARAMETERS 
C     ----------------- 
C     T           T-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN T=TEND).
C
C     Q(NQ)        NUMERICAL APPROXIMATION OF POSITION VECTOR AT T.
C 
C     V(NV)        NUMERICAL APPROXIMATION OF VELOCITY VECTOR AT T.
C
C     U(NV)        NUMERICAL APPROXIMATION OF EXTERNAL DYNAMIC VECTOR AT T.
C
C     A(NV)        NUMERICAL APPROXIMATION OF ACCELERATION VECTOR AT T.
C
C     RLAM(NL)     NUMERICAL APPROXIMATION OF LAGRANGE MULTIPLIERS AT T.
C 
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP.
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS SINGULAR.
C                   IDID=-5  INITIAL PROJECTION: NO CONVERGENCE  
C-----------------------------------------------------------------------
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(NQ),V(NV),U(NU),A(NV),RLAM(NL) 
      DIMENSION ATOL(1),RTOL(1),WK(LWK),IWK(LIWK)
      LOGICAL ARRET
      EXTERNAL FPROB,SOLOUT
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS 
C *** *** *** *** *** *** ***
      ARRET=.FALSE.
      MODE=IWK(14)
C -------- MODE, THE CHOICE OF LINEAR ALGEBRA
      IF ((IWK(14).LE.-1).OR.(IWK(14).GE.6)) THEN
        WRITE(6,*)' WRONG CHOICE OF IWK(14) (MODE):'
     &             ,IWK(14)
        ARRET=.TRUE.
      ELSE 
        MODE=IWK(14)
      END IF
      LDG = IWK(1)  
      LDF = IWK(2)
      NM=NV+NL   
      NMRC = IWK(21)
      NBLK = IWK(22)
      NZA = NBLK*NMRC**2+LDF+LDG
      IF (MODE.EQ.4) LDF=LDG 
      IF (MODE.GE.4) THEN
        IF (LDG.LE.0) THEN
          WRITE(6,*)' IWK(1) (LDG) MUST BE POSITIVE'
          ARRET=.TRUE.
        END IF
        IF (LDF.LE.0) THEN
          WRITE(6,*)' IWK(2) (LDF) MUST BE POSITIVE'
          ARRET=.TRUE.
        END IF
        IF ((NMRC.EQ.0).OR.(NBLK.EQ.0)) THEN
          WRITE(6,*)' IWK(21) (NMRC) AND IWK(22) (NBLK) MUST
     &        BE POSITIVE'
          ARRET=.TRUE.
        END IF
        IF (IWK(28).EQ.0) IWK(28)=6
        IF (IWK(25).LT.13*NM) THEN
           WRITE(6,*)' INTEGER WORK SPACE (IWK(25))
     &           FOR MA28 TOO SMALL'
           IWK(25)=13*NM
        END IF
        IF (IWK(26).LT.(4*NZA+NM)) THEN
           WRITE(6,*)' REAL WORK SPACE (IWK(26))
     &           FOR MA28 TOO SMALL'
           IWK(26)=4*NZA+NM
        END IF
      END IF
      IF (MODE.LE.3) THEN
        LDG=NL
        LDF=NV
        LDA=NM 
        NDIM=NV
        MDIM=NL
        NMDIM=NM
      ELSE
        LDA = NBLK*NMRC
        NDIM=1
        MDIM=1
        NMDIM=NMRC
      END IF
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF(IWK(11).EQ.0)THEN
        NMAX=100000
      ELSE
        NMAX=IWK(11)
      END IF
C -------- IPCIV    SWITCH FOR INITIAL PROECTION
      IF (IWK(12).EQ.1) THEN
        IPCIV=1
      ELSE
        IPCIV=0
      ENDIF
C -------- IGII AND IACC
      IF (IWK(15).EQ.1) THEN
        IACC=1
      ELSE
        IACC=0
      ENDIF
      IF (IWK(16).EQ.1) THEN
        IGII=1
      ELSE
        IGII=0
      ENDIF
C -------- IPRO    NUMBERS OF STEPS BETWEEN TWO PROECTIONS
      IPROJ = IWK(13)
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  
      IF(WK(1).EQ.0.D0)THEN
        UROUND=1.D-16
      ELSE
        UROUND=WK(1)
        IF(UROUND.LE.1.D-35.OR.UROUND.GE.1.D0)THEN
         WRITE(6,*)' WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:'
     &                            ,WK(1)
         ARRET=.TRUE.
        END IF
      END IF
C -------  SAFETY FACTOR -------------
      IF(WK(2).EQ.0.D0)THEN
        SAFE=0.9D0
      ELSE
        SAFE=WK(2)
        IF(SAFE.GE.1.D0.OR.SAFE.LE.1.D-4)THEN
         WRITE(6,*)' CURIOUS INPUT FOR SAFETY FACTOR WK(2)='
     &                            ,WK(2)
         ARRET=.TRUE.
        END IF
      END IF
C -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
      IF(WK(3).EQ.0.D0)THEN
        F1=0.2D0
      ELSE
        F1=WK(3)
      END IF
      IF(WK(4).EQ.0.D0)THEN
        F2=10.D0
      ELSE
        F2=WK(4)
      END IF
C -------- MAXIMAL STEP SIZE
      IF(WK(6).EQ.0.D0)THEN
        HMAX=TEND-T
      ELSE
        HMAX=WK(6)
      END IF 
C -------- GUSTAFFSON STRATEGIE
      BETA = WK(7)
      IF(WK(8).EQ.0.D0)THEN
        ALPHA = 1/5.d0
      ELSE
        ALPHA=WK(8)
      END IF 
C -- WK SPACE
      LXUMF = NM
      LIUMF = 13*NM+4*nza
      LIPAR = 9
      LISTA = 9
      LRDO = 2 + 5 * (NQ+NV+NU)
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WK -----
      IQ1 = 11
      IQ2 = IQ1+NQ  
      IQQ2 = IQ2+NQ
      IQDOT1 = IQQ2+NQ
      IQDOT2 = IQDOT1+NQ
      IQDOT3 = IQDOT2+NQ
      IQDOT4 = IQDOT3+NQ
      IQDOT5 = IQDOT4+NQ
      IQDOT6 = IQDOT5+NQ
      IQDOT = IQDOT6+NQ
      IV1 = IQDOT+NQ
      IV2 = IV1+NV  
      IU1 = IV2+NV
      IU2 = IU1+NU
      IVP1 = IU2+NU
      IVP2 = IVP1+NV  
      IVP3 = IVP2+NV
      IVP4 = IVP3+NV
      IVP5 = IVP4+NV
      IVP6 = IVP5+NV
      IXL = IVP6+NV
      IUDOT1 = IXL+NL
      IUDOT2 = IUDOT1+NU
      IUDOT3 = IUDOT2+NU
      IUDOT4 = IUDOT3+NU
      IUDOT5 = IUDOT4+NU
      IUDOT6 = IUDOT5+NU
      IUDOT = IUDOT6+NU
      IGQ = IUDOT+NU
      IGQ0 = IGQ+LDG*NDIM
      IGQ1 = IGQ0+LDG*NDIM 
      IB = IGQ1+LDG*NDIM
      IX0 = IB+LDA*NMDIM
      ITEMP = IX0+NM
      IAM = ITEMP +NV
      IGT = IAM+LDA*NMDIM 
      IG = IGT+NL
      IFL = IG+NL
      IAV = IFL+LDF*MDIM
      IXUMF = IAV+2*NZA
       IGD = IXUMF+LXUMF
      IGTD = IGD+LDG*NDIM
      IQD = IGTD+NL
      IVD = IQD+NQ
      IUD = IVD+NV 
      IDO = IUD+NU
C ------ TOTAL STORAGE REQUIREMENT -----------
      IS=IDO+LRDO
      IF(IS.GT.LWK)THEN
       WRITE(6,*)' INSUFFICIENT STORAGE FOR WK, MIN. LWK=',IS
       ARRET=.TRUE.
      END IF
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN IWK -----
      IPA = 21
      ISTAT = 31
      IIP = 41
      ING = 41+NM
      ING0 = ING+2*LDG
      ING1 = ING0+2*LDG
      INGD = ING1+2*LDG
      INFL = INGD+2*LDG
      INAM = INFL+2*LDF
      INUMF =INAM+2*NZA
C ------ TOTAL STORAGE REQUIREMENT -----------
      IS=INUMF+LIUMF	
      IF(IS.GT.LIWK)THEN
       WRITE(6,*)' INSUFFICIENT STORAGE FOR IWK, MIN. LIWK=',IS
       ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
       IDID=-1
       RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL PHCOR56(
     & NQ,NV,NU,NL,NM,NDIM,MDIM,NMDIM,NZA,LRDO,LIPAR,
     & LISTA,LIUMF,LXUMF,LDG,LDF,LDA,MODE,NMAX,IPCIV,IPROJ,IOUT,
     & IGII,IACC,ITOL,IWK(IPA),IWK(ISTAT),IWK(IIP),IWK(ING),
     & IWK(ING0),IWK(ING1),IWK(INGD),IWK(INFL),IWK(INAM),IWK(INUMF),
     & FPROB,SOLOUT,UROUND,T,TEND,H,HMAX,RTOL,ATOL,SAFE,
     & ALPHA,BETA,F1,F2,Q,V,U,A,RLAM,WK(IQ1),WK(IQ2),
     & WK(IQQ2),WK(IQDOT1),WK(IQDOT2),
     & WK(IQDOT3),WK(IQDOT4),WK(IQDOT5),WK(IQDOT6),
     & WK(IQDOT),WK(IV1),WK(IV2),
     & WK(IU1),WK(IU2),
     & WK(IVP1),WK(IVP2),
     & WK(IVP3),WK(IVP4),WK(IVP5),WK(IVP6),WK(IXL),
     & WK(IUDOT1),WK(IUDOT2),WK(IUDOT3),WK(IUDOT4),WK(IUDOT5),
     & WK(IUDOT6),WK(IUDOT),WK(IGQ),WK(IGQ0),
     & WK(IGQ1),WK(IB),WK(IX0),WK(ITEMP),WK(IAM),
     & WK(IGT),WK(IG),WK(IFL),WK(IAV),WK(IXUMF),
     & WK(IGD),WK(IGTD),
     & WK(IQD),WK(IVD),WK(IUD),WK(IDO))
C ----------- RETURN -----------
      RETURN
      END
C
C
C  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
C
      SUBROUTINE PHCOR56(
     & NQ,NV,NU,NL,NM,NDIM,MDIM,NMDIM,NZA,LRDO,
     & LIPAR,LISTA,LIUMF,LXUMF,LDG,LDF,LDA,MODE,NMAX,
     & IPCIV,IPROJ,IOUT,IGII,IACC,ITOL,IPAR,ISTAT,IP,INDG,
     & INDG0,INDG1,INDGD,INDFL,INDA,IUMF,FPROB,SOLOUT,
     & UROUND,T,TEND,H,HMAX,RTOL,ATOL,SAFE,ALPHA,BETA,
     & FAC1,FAC2,Q,V,U,A,RLAM,Q1,Q2,QQ2,
     & QDOT1,QDOT2,QDOT3,QDOT4,QDOT5,QDOT6,QDOT,V1,V2,
     & U1,U2,VP1,VP2,
     & VP3,VP4,VP5,VP6,XL,UDOT1,UDOT2,UDOT3,UDOT4,
     & UDOT5,UDOT6,UDOT,GQ,GQ0,GQ1,B,X0,TEMP,
     & AM,GT,G,FL,AVALUE,XUMF,
     & GD,GTD,QD,VD,UD,DOWK)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR PHEM56
C     PARAMETERS SAME AS IN PHEM56 WITH WORKSPACE ADDED 
C ---------------------------------------------------------- 
C         DECLARATIONS 
C ---------------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(NQ),Q1(NQ),Q2(NQ)
      DIMENSION QQ2(NQ),V(NV),V1(NV),V2(NV),VP1(NV),VP2(NV)
      DIMENSION VP3(NV),VP4(NV),VP5(NV),VP6(NV)
      DIMENSION GQ(LDG,NDIM),GQ0(LDG,NDIM),GQ1(LDG,NDIM)
      DIMENSION B(LDA,NMDIM),AM(LDA,NMDIM),INDG(2*LDG)
      DIMENSION X0(NM),IP(NM),TEMP(NV),AVALUE(2*NZA)
      DIMENSION GT(NL),G(NL)
      DIMENSION GTD(NL),QD(NQ),VD(NV),DOWK(LRDO) 
      DIMENSION ATOL(1),RTOL(1),XUMF(LXUMF),FL(LDF,MDIM)
      DIMENSION GD(LDG,NDIM),IPAR(LIPAR)
      DIMENSION INDG0(2*LDG),INDG1(2*LDG),INDGD(2*LDG),A(NV)
      DIMENSION INDFL(2*LDF),INDA(2*NZA),IUMF(LIUMF),ISTAT(LISTA)
      DIMENSION UDOT1(NU),UDOT2(NU),UDOT3(NU),UDOT4(NU)
      DIMENSION UDOT5(NU),UDOT6(NU),UDOT(NU),RLAM(NL)
      DIMENSION QDOT1(NQ),QDOT2(NQ),QDOT3(NQ),QDOT4(NQ)
      DIMENSION QDOT5(NQ),QDOT6(NQ),QDOT(NQ)
      DIMENSION U(NU),U1(NU),U2(NU),UD(NU),XL(NL)
      LOGICAL REJECT,ACCEPT,LAST
      EXTERNAL FPROB,SOLOUT
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** *** 
       CALL  PCOEF56(
     &    C2,C3,C4,C5,C6,C7,C10,A21,A31,A32,
     &    A41,A42,A43,A51,A52,A53,A54,A61,A62,A63,A64,A65,
     &    A71,A72,A73,A74,A75,A76,BB1,BB3,BB4,BB5,BB6,BB7,
     &    B1,B2,B3,B4,B5,B6,B7,BC1,BC3,BC4,BC5,BC6,BC7,
     &    CC1,CC2,CC3,CC4,CC5,CC6,CC7,AA11,AA21,AA22,
     &    AA31,AA32,AA33,AA41,AA42,AA43,AA44,AA51,AA52,AA53,
     &    AA54,AA55, AA61,AA62,AA63,AA64,AA65,AA66,AA71,AA72,
     &    AA73,AA74,AA75,AA76,AA77) 
      DO I=1,9
        ISTAT(I)=0
      END DO
      NSTEP=0
      NQ2=2*NQ
      NQ3=3*NQ
      NQ4=4*NQ
      NV2=2*NV
      NV3=3*NV
      NV4=4*NV
      NU2=2*NU
      NU3=3*NU
      NU4=4*NU
      NQVU=NQ+NV+NU
      NQVU2=2*NQVU
      NQVU3=3*NQVU
      NQVU4=4*NQVU
      NBLK = IPAR(2)
      NMRC = IPAR(1)
cc
cc
      IPAR(9)=1
C --  NBS : COUNTS THE STEPS BETWEEN THE PROJECTIONS
      NBS = 0  
      LAST=.FALSE.
      ACCEPT = .TRUE.
      REJECT=.FALSE.
      FACOLD=1.D0 
      FACC1=1.D0/FAC1
      FACC2=1.D0/FAC2
      POSNEG=SIGN(1.D0,TEND-T)
      IRTRN=1
      HMAX=ABS(HMAX)
      H=MIN(MAX(1.D-4,ABS(H)),HMAX)
      H=SIGN(H,POSNEG)
      TOLD=T
      CALL FPROB(10,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &      IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &      INDFL(LDF+1),T,Q,V,U,XL,
     &      G,GQ,X0(1),X0(NV+1),GT,FL,QDOT1,UDOT1,B)
      ISTAT(5)=1
      ISTAT(6)=1
C --- PROJECTION TO ENSURE CONSISTENT INITIAL VALUES
        CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &         NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG,
     &         INDG,INDFL,INDA,IP,IUMF,XUMF,
     &         B,GQ,GQ,FL,AVALUE,IER)
        IF (IER.NE.0) GOTO 176
        IF (IPCIV.EQ.1) THEN
          CALL APROJ(MODE,NQ,NV,NU,NL,NM,NDIM,MDIM,NMDIM,
     *         NBLK,NMRC,LDG,LDF,LDA,NZA,LIPAR,LIUMF,LXUMF,
     *         LISTA,ISTAT,IPAR,IUMF,INDA,INDG,INDFL,FL,XUMF,
     *         AVALUE,T,FBROB,Q,Q,Q2,QD,V,V,V2,U,XL,UDOT,G,GT,
     *         GQ,AM,B,X0,IP,ATOL,RTOL,ITOL,ACCEPT)
          ISTAT(6) = ISTAT(6)+2 
          ISTAT(5) = ISTAT(5)+1 
          IF (.NOT.(ACCEPT)) GOTO 179
        END IF
C   Consistent initial values for A and RLAM at initial values
          ifcn = 5
          CALL FPROB(ifcn,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &       IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &       INDFL(LDF+1),T,Q,V,U,XL,
     &       G,GQ,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
          IF (ifcn.EQ.0) CALL GIINUM(MODE,NQ,NV,NU,NL,NDIM,
     &       MDIM,NMDIM,NM,
     &       LDG,LDF,LDA,NBLK,NMRC,IPAR(3),IPAR(4),INDG,INDGD,
     &       INDFL,T,UROUND,FPROB,Q,QD,QDOT1,V,U,GQ,GD,GT,GTD,
     &       FL,AM,X0(NV+1),G)
          CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &       IP,NZA,AVALUE,XUMF,B,X0,A,RLAM)
          ISTAT(4)=1
          ISTAT(8)=1
C
       IF(IOUT.GE.1) THEN
        DOWK(1)=TOLD
        DOWK(2)=H
        CALL SOLOUT(ISTAT(2)+1,NQ,NV,NU,NL,LRDO,
     &              Q,V,U,A,RLAM,DOWK)
      END IF
C --  FIRST STAGE : Q,V,U,QDOT1,A,RLAM CONTAIN THE INITIAL VALUES
      DO  I=1,NV
        Q1(I) = Q(I)
        V1(I) = V(I)
        VP1(I) = A(I)
      END DO
      DO I=NV+1,NQ
        Q1(I) = Q(I)
      END DO   
      DO I=1,NU
        U1(I) = U(I)
      END DO   
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,Q1,V1,U1,RLAM,
     &     G,GQ,X0(1),X0(NV+1),GT,FL,QDOT1,UDOT1,B)
C --- BASIC INTEGRATION STEP  
 909  continue
      IF(NSTEP.GT.NMAX) GOTO 178
      IF (0.1D0*ABS(H).LE.ABS(T)*UROUND) GOTO 177
      IF((T+H-TEND)*POSNEG.GE.0.D0) THEN
        H=TEND-T
        LAST=.TRUE.
      END IF
      ISTAT(1)=ISTAT(1)+1 
      ACCEPT = .TRUE.
C --  2d STAGE
      ipar(3)=0
cc
      TCH = T+C2*H      
      DO  I=1,NU
        U2(I) = U1(I) + H*A21*UDOT1(I)
      END DO
      DO I=1,NQ
        Q2(I) = Q1(I) + H*A21*QDOT1(I)
      END DO
      DO  I=1,NV
        V2(I) = V1(I) + H*A21*VP1(I)
        TEMP(I) = V1(I) + H*AA21*VP1(I)
      END DO
      CALL FPROB(1,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,RLAM,
     &     G,GQ1,X0(1),X0(NV+1),GT,FL,QDOT2,UDOT1,B)
      DO  I=1,NQ
        QQ2(I) = Q1(I) + H*(AA21*QDOT1(I)+AA22*QDOT2(I))
      END DO
      TTCH = T+CC2*H
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TTCH,QQ2,V2,U2,RLAM,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT2,UDOT1,B)
      FAC = -1.D0/(H*AA22)
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG0,
     &     GQ0,FAC,TEMP,GT,X0(NV+1))	
      CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &     NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG1,INDG0,
     &     INDFL,INDA,IP,IUMF,XUMF,B,GQ1,GQ0,FL,AVALUE,IER)
      IF (IER.NE.0) GOTO 176
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &     IP,NZA,AVALUE,XUMF,B,X0,VP2,XL)
C - II.3D STAGE
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT2,UDOT2,B)
      DO  I=1,NU
        U2(I) = U1(I) + H*(A31*UDOT1(I)+A32*UDOT2(I))
      END DO
      DO  I=1,NQ
        Q2(I) = Q1(I) + H*(A31*QDOT1(I)+A32*QDOT2(I))
      END DO
      TCH = T+C3*H
      DO  I=1,NV
        V2(I) = V1(I) + H*(A31*VP1(I)+A32*VP2(I))
        TEMP(I) = V1(I) + H*(AA31*VP1(I)+AA32*VP2(I))
      END DO
      CALL FPROB(1,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT3,UDOT2,B)
      DO  I=1,NQ
        QQ2(I) = Q1(I) + H*(AA31*QDOT1(I)+AA32*QDOT2(I)+
     &                      AA33*QDOT3(I))
      END DO
      TTCH = T+CC3*H
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TTCH,QQ2,V2,U2,XL,
     &     G,GQ1,X0(1),X0(NV+1),GT,FL,QDOT3,UDOT1,B)
      FAC = -1.D0/(H*AA33)
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG1,
     &     GQ1,FAC,TEMP,GT,X0(NV+1))	
      CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &     NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG0,INDG1,
     &     INDFL,INDA,IP,IUMF,XUMF,B,GQ0,GQ1,FL,AVALUE,IER)
      IF (IER.NE.0) GOTO 176
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &     IP,NZA,AVALUE,XUMF,B,X0,VP3,XL)
C
C - II.4th STAGE
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT3,UDOT3,B)
      DO  I=1,NU
        U2(I) = U1(I) + H*(A41*UDOT1(I)+A42*UDOT2(I)+
     &                     A43*UDOT3(I))
      END DO
      DO  I=1,NQ
        Q2(I) = Q1(I) + H*(A41*QDOT1(I)+A42*QDOT2(I)+A43*QDOT3(I))
      END DO
      TCH = T+C4*H
      DO  I=1,NV
        V2(I) = V1(I) + H*(A41*VP1(I)+A42*VP2(I)+
     &                    A43*VP3(I))
        TEMP(I) = V1(I) + H*(AA41*VP1(I)+AA42*VP2(I)+
     &                       AA43*VP3(I))
      END DO
      CALL FPROB(1,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ1,X0(1),X0(NV+1),GT,FL,QDOT4,UDOT3,B)
      DO  I=1,NQ
        QQ2(I) = Q1(I) + H*(AA41*QDOT1(I)+AA42*QDOT2(I)+
     &                      AA43*QDOT3(I)+AA44*QDOT4(I))
      END DO
      TTCH = T+CC4*H
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TTCH,QQ2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT4,UDOT3,B)
      FAC = -1.D0/(H*AA44)
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG0,
     &     GQ0,FAC,TEMP,GT,X0(NV+1))	
      CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &     NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG1,INDG0,
     &     INDFL,INDA,IP,IUMF,XUMF,B,GQ1,GQ0,FL,AVALUE,IER)
      IF (IER.NE.0) GOTO 176
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &     IP,NZA,AVALUE,XUMF,B,X0,VP4,XL)
C - II.5th STAGE
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT4,UDOT4,B)
      DO  I=1,NU
        U2(I) = U1(I) +  H*(A51*UDOT1(I)+A52*UDOT2(I)+
     &                      A53*UDOT3(I)+A54*UDOT4(I))
      END DO
      DO  I=1,NQ
         Q2(I) = Q1(I) +  H*(A51*QDOT1(I)+A52*QDOT2(I)+
     &                      A53*QDOT3(I)+A54*QDOT4(I))
      END DO
      TCH = T+C5*H
      DO  I=1,NV
        V2(I) = V1(I) +  H*(A51*VP1(I)+A52*VP2(I)+
     &                      A53*VP3(I)+A54*VP4(I))
        TEMP(I) = V1(I) +  H*(AA51*VP1(I)+AA52*VP2(I)+
     &                      AA53*VP3(I)+AA54*VP4(I))
      END DO
      CALL FPROB(1,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT5,UDOT4,B)
      DO  I=1,NQ
        QQ2(I) = Q1(I) + H*(AA51*QDOT1(I)+AA52*QDOT2(I)+
     &           AA53*QDOT3(I) + AA54*QDOT4(I)+AA55*QDOT5(I) )
      END DO
      TTCH = T+CC5*H
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TTCH,QQ2,V2,U2,XL,
     &     G,GQ1,X0(1),X0(NV+1),GT,FL,QDOT5,UDOT4,B)
      FAC = -1.D0/(H*AA55)
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG1,
     &     GQ1,FAC,TEMP,GT,X0(NV+1))	
      CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &     NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG0,INDG1,
     &     INDFL,INDA,IP,IUMF,XUMF,B,GQ0,GQ1,FL,AVALUE,IER)
      IF (IER.NE.0) GOTO 176
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &     IP,NZA,AVALUE,XUMF,B,X0,VP5,XL)
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT5,UDOT5,B)
C - II.6th STAGE
      DO  I=1,NU
        U2(I) = U1(I) +  H*(A61*UDOT1(I)+A62*UDOT2(I)+
     &           A63*UDOT3(I) + A64*UDOT4(I)+A65*UDOT5(I) )
      END DO
      DO  I=1,NQ
        Q2(I) = Q1(I) +  H*(A61*QDOT1(I)+A62*QDOT2(I)+
     &            A63*QDOT3(I) + A64*QDOT4(I)+A65*QDOT5(I) )
      END DO
      TCH = T+C6*H
      DO  I=1,NV
        V2(I) = V1(I) +  H*(A61*VP1(I)+A62*VP2(I)+
     &           A63*VP3(I) + A64*VP4(I)+A65*VP5(I) )
        TEMP(I) = V1(I) +  H*(AA61*VP1(I)+AA62*VP2(I)+
     &           AA63*VP3(I) + AA64*VP4(I)+AA65*VP5(I) )
      END DO
      CALL FPROB(1,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ1,X0(1),X0(NV+1),GT,FL,QDOT6,UDOT5,B)
      DO  I=1,NQ
        QQ2(I) = Q1(I) + H*(AA61*QDOT1(I)+AA62*QDOT2(I)+
     &                      AA63*QDOT3(I)+AA64*QDOT4(I)+
     &                      AA65*QDOT5(I)+AA66*QDOT6(I))
      END DO
      TTCH = T+CC6*H
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TTCH,QQ2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT6,UDOT5,B)
      FAC = -1.D0/(H*AA66)
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG0,
     &     GQ0,FAC,TEMP,GT,X0(NV+1))	
      CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &     NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG1,INDG0,
     &     INDFL,INDA,IP,IUMF,XUMF,B,GQ1,GQ0,FL,AVALUE,IER)
      IF (IER.NE.0) GOTO 176
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &     IP,NZA,AVALUE,XUMF,B,X0,VP6,XL)
C - LAST STAGE
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TCH,Q2,V2,U2,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT6,UDOT6,B)
      DO  I=1,NU
        U(I) = U1(I) +  H*(B1*UDOT1(I)+
     &                      B3*UDOT3(I)+B4*UDOT4(I)+
     &                      B5*UDOT5(I)+B6*UDOT6(I))
      END DO
      DO  I=1,NQ
C         Q2(I) = Q1(I) +  H*(B1*QDOT1(I)+
C     &                      B3*QDOT3(I)+B4*QDOT4(I)+
C     &                      B5*QDOT5(I)+B6*QDOT6(I))
      Q(I) = QQ2(I)
      END DO
      TH = T+H
      DO  I=1,NV
C        V2(I) = V1(I) +   H*(B1*VP1(I)+
C     &                      B3*VP3(I)+B4*VP4(I)+
C     &                      B5*VP5(I)+B6*VP6(I))
         V(I) = TEMP(I) + H * AA66 * VP6(I)
      END DO
C   FSAL STAGE
C --------- COMPUTE THE ACCELERATION VP AND XL
        IF (IACC.EQ.1) THEN
         CALL FPROB(2,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &        IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &        INDFL(LDF+1),TH,Q,V,U,XL,
     &        G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
          ifcn = 7
          CALL FPROB(ifcn,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &         IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &         INDFL(LDF+1),TH,Q,V,U,XL,
     &         G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
          IF (ifcn.eq.0) THEN
           CALL GIINUM(MODE,NQ,NV,NU,NL,NDIM,MDIM,NMDIM,NM,
     &          LDG,LDF,LDA,NBLK,NMRC,IPAR(3),IPAR(4),INDG0,
     &          INDGD,INDFL,TH,UROUND,FPROB,Q,QD,QDOT,V,U,GQ0,
     &          GD,GT,GTD,FL,AM,X0(NV+1),G)
          END IF
          CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,NZA,
     &         LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG0,INDG0,
     &         INDFL,INDA,IP,IUMF,XUMF,B,GQ0,GQ0,FL,AVALUE,IER)
          IF (IER.NE.0) GOTO 176
          CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,IP,
     &         NZA,AVALUE,XUMF,B,X0,A,RLAM)
          ISTAT(4)=ISTAT(4)+1
          ISTAT(6)=ISTAT(6)+1
          ISTAT(8)=ISTAT(8)+1
      ELSE
C  --------- COMPUTE VP AND XL WITH AN ADDITIONAL STEP
      DO I=1,NV
        TEMP(I) = V1(I) +  H*(AA71*VP1(I)+AA72*VP2(I)+
     &                      AA73*VP3(I)+AA74*VP4(I)+
     &                      AA75*VP5(I)+AA76*VP6(I))
      END DO
      CALL FPROB(1,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TH,Q,V,U,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT,UDOT6,B)
      DO  I=1,NQ
        QQ2(I) = Q1(I) + H*(AA71*QDOT1(I)+AA72*QDOT2(I)+
     &                      AA73*QDOT3(I)+AA74*QDOT4(I)+
     &         AA75*QDOT5(I)+AA76*QDOT6(I)+AA77*QDOT(I))
      END DO
      TTCH = T+CC7*H
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TTCH,QQ2,V,U,XL,
     &     G,GQ1,X0(1),X0(NV+1),GT,FL,QDOT,UDOT6,B)
      FAC = -1.D0/(H*AA77)
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG1,
     &     GQ1,FAC,TEMP,GT,X0(NV+1))	
      CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &     NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG0,INDG1,
     &     INDFL,INDA,IP,IUMF,XUMF,B,GQ0,GQ1,FL,AVALUE,IER)
      IF (IER.NE.0) GOTO 176
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &     IP,NZA,AVALUE,XUMF,B,X0,A,RLAM)
          ISTAT(4)=ISTAT(4)+1
          ISTAT(5) = ISTAT(5)+1
          ISTAT(6)=ISTAT(6)+1
          ISTAT(8)=ISTAT(8)+1
      ENDIF
C -----------------
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TH,Q,V,U,RLAM,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
C -- Statistics
      ISTAT(4) = ISTAT(4)+5
      ISTAT(5) = ISTAT(5)+5
      ISTAT(6) = ISTAT(6)+5
      ISTAT(8) = ISTAT(8)+5
C --- Error Estimation  
      ERR = 0.D0
c
      DO  I=1,NQ
        IF (ITOL.EQ.0) THEN
          SK1 = ATOL(1)+RTOL(1)*MAX(DABS(Q(I)),ABS(Q1(I)))
        ELSE
          SK1 = ATOL(I)+RTOL(I)*MAX(ABS(Q(I)),ABS(Q1(I)))
        END IF
        Q2(I) = Q1(I) + H*(BB1*QDOT1(I)+BB3*QDOT3(I)+
     &    BB4*QDOT4(I)+BB5*QDOT5(I)+BB6*QDOT6(I)+BB7*QDOT(I))
        ERR = ERR+((Q(I)-Q2(I))/SK1)**2
      END DO
      DO  I=1,NV
        V2(I) = V1(I) + H*(BB1*VP1(I)+BB3*VP3(I)+
     &    BB4*VP4(I)+BB5*VP5(I)+BB6*VP6(I)+BB7*A(I))
      END DO
c
c     One Newton iteration of projection of V2 = $\tilde v_1$. 
c                                           VP2= $\tilde v_1 - \hat v_1$
      DO I=1,NV
        X0(I) = 0.D0
      END DO
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG1(1),INDG1(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TH,Q2,V2,U,RLAM,
     &     G,GQ1,X0(1),X0(NV+1),GT,FL,QDOT,UDOT6,B)
      FAC = -1.D0
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG1,
     &     GQ1,FAC,V2,GT,X0(NV+1))	
      IF (IER.NE.0) GOTO 176
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,
     &     IP,NZA,AVALUE,XUMF,B,X0,VP2,XL)
      ISTAT(5) = ISTAT(5)+1
      ISTAT(8) = ISTAT(8)+1
c
c     end of projection
c
      DO I=1,NL
        XL(I) = RLAM(I) + XL(I)/(H*BB7)
      END DO
      CALL FPROB(0,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG0(1),INDG0(LDG+1),INDFL(1),
     &     INDFL(LDF+1),TH,Q,V,U,XL,
     &     G,GQ0,X0(1),X0(NV+1),GT,FL,QDOT2,UDOT2,B)
      DO I=1,NV
        IF (ITOL.EQ.0) THEN
          SK2 = ATOL(1)+RTOL(1)*MAX(DABS(V(I)),ABS(V1(I)))
        ELSE
          SK2 = ATOL(NQ+I)+RTOL(NQ+I)*MAX(ABS(V(I)),ABS(V1(I)))
        END IF
        ERR = ERR+((V(I)-V2(I)-VP2(I))/SK2)**2
      END DO
      NQV=NQ+NV
      DO  I=1,NU
        IF (ITOL.EQ.0) THEN
          SK1 = ATOL(1)+RTOL(1)*MAX(DABS(U(I)),ABS(U1(I)))
        ELSE
          SK1 = ATOL(NQV+I)+RTOL(NQV+I)*MAX(ABS(U(I)),ABS(U1(I)))
        END IF
        U2(I) = U1(I) + H*(BB1*UDOT1(I)+BB3*UDOT3(I)+
     &    BB4*UDOT4(I)+BB5*UDOT5(I)+BB6*UDOT6(I)+BB7*UDOT2(I))
        ERR = ERR+((U(I)-U2(I))/SK1)**2
      END DO
      ERR = SQRT(ERR/(NQV+NU))  
c  
      EPS=1.D-13
      IF (ERR.LT.EPS) THEN
        HNEW=1.2D0*H
        GOTO 222
      END IF
C --- COMPUTATION OF HNEW WITH GUSTAFSSON STABILIZATION
      FAC11=ERR**ALPHA
      FAC=FAC11/FACOLD**BETA
C --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
      FAC=MAX(FACC2,MIN(FACC1,FAC/SAFE))
      HNEW=H/FAC 
      IF (ERR.GT.1.D0) ACCEPT=.FALSE.
c
 222  CONTINUE 
c
      IF (ACCEPT) THEN 
        NBS = NBS+1
        IF (NBS.EQ.IPROJ) THEN
          NBS = 0
C -- PROJECTION ON MANIFOLD DEFINED BY g(q,t)=0
          CALL APROJ(MODE,NQ,NV,NU,NL,NM,NDIM,MDIM,NMDIM,NBLK,
     *         NMRC,LDG,LDF,LDA,NZA,LIPAR,LIUMF,LXUMF,LISTA,
     *         ISTAT,IPAR,IUMF,INDA,INDG0,INDFL,FL,XUMF,
     *         AVALUE,T,FPROB,Q,Q1,Q2,QD,V,V1,V2,U,XL,UDOT,
     *         G,GT,GQ0,AM,B,X0,IP,ATOL,RTOL,ITOL,ACCEPT)
           ISTAT(6) = ISTAT(6)+2 
           ISTAT(5) = ISTAT(5)+1 
        END IF 
        IF (ACCEPT) THEN
C --- STEP IS ACCEPTED  
         FACOLD=MAX(ERR,1.D-4)
         ISTAT(2)=ISTAT(2)+1
C
c Dense output
        IF (IOUT.EQ.2) THEN
C -- COMPUTE QD ,VD AND UD AT theta = 19/30
c
      DO  I=1,NQ
        QD(I) = Q1(I) + H*(BC1*QDOT1(I)+BC3*QDOT3(I)+
     &    BC4*QDOT4(I)+BC5*QDOT5(I)+BC6*QDOT6(I)+BC7*QDOT(I))
      END DO
      DO  I=1,NV
        VD(I) = V1(I) + H*(BC1*VP1(I)+BC3*VP3(I)+
     &    BC4*VP4(I)+BC5*VP5(I)+BC6*VP6(I)+BC7*A(I))
      END DO
c
      DO  I=1,NU
        UD(I) = U1(I) + H*(BC1*UDOT1(I)+BC3*UDOT3(I)+
     &    BC4*UDOT4(I)+BC5*UDOT5(I)+BC6*UDOT6(I)+BC7*UDOT(I))
      END DO
C-- STATISTICS
          ISTAT(4)=ISTAT(4)+1
          ISTAT(5)=ISTAT(5)+1
          ISTAT(6)=ISTAT(6)+1
          ISTAT(8)=ISTAT(8)+1
C --- PREPARES DOWK
          DO I=1,NV
            IP2=I+2
            DOWK(IP2)=Q1(I)
            DOWK(IP2+NQVU)=QDOT1(I)
            DOWK(IP2+NQVU2)=QDOT(I)
            DOWK(IP2+NQVU3)=QD(I)-Q1(I)
            DOWK(IP2+NQVU4)=Q(I)-Q1(I)
            IP3=IP2+NQ
            DOWK(IP3)=V1(I)
            DOWK(IP3+NQVU)= VP1(I)
            DOWK(IP3+NQVU2)= A(I)
            DOWK(IP3+NQVU3)=VD(I)-V1(I)
            DOWK(IP3+NQVU4)=V(I)-V1(I) 
          END DO 
          DO I=NV+1,NQ
            IP2=I+2
            DOWK(IP2)=Q1(I)
            DOWK(IP2+NQVU)=QDOT1(I)
            DOWK(IP2+NQVU2)=QDOT(I)
            DOWK(IP2+NQVU3)=QD(I)-Q1(I)
            DOWK(IP2+NQVU4)=Q(I)-Q1(I)
          END DO
          NNN=NQ+NV+2
          DO I=1,NU
            IP4=NNN+I
            DOWK(IP4)=U1(I)
            DOWK(IP4+NQVU)=UDOT1(I)
            DOWK(IP4+NQVU2)=UDOT(I)
            DOWK(IP4+NQVU3)=UD(I)-U1(I)
            DOWK(IP4+NQVU4)=U(I)-U1(I)
          END DO
          DOWK(1)=T
          DOWK(2)=H
          CALL SOLOUT(ISTAT(2)+1,NQ,NV,NU,NL,LRDO,
     &              Q,V,U,A,RLAM,DOWK)
          END IF
c   End of dense output
        IF(IOUT.EQ.3) THEN
          DOWK(1)=T
          DOWK(2)=H
          CALL SOLOUT(ISTAT(2)+1,NQ,NV,NU,NL,LRDO,
     &              Q,V,U,A,RLAM,DOWK)
        END IF
        IF (IRTRN.LT.0) GOTO 176
        IF(ABS(HNEW).GT.HMAX)HNEW=POSNEG*HMAX  
        IF(REJECT) HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
        REJECT=.FALSE. 
        TOLD=T
        T=TH 
        IF (LAST) THEN
           IDID=1
           RETURN
        END IF
        DO I=1,NV
          Q1(I) = Q(I)
          QDOT1(I) = QDOT(I)
          V1(I) = V(I) 
          VP1(I) = A(I)
        END DO
        DO I=NV+1,NQ
          Q1(I) = Q(I)
          QDOT1(I) = QDOT(I)
        END DO
        DO I=1,NU
          U1(I) = U(I)
          UDOT1(I) = UDOT(I)
        END DO
      ELSE  
C -- STEP IS REJECTED AFTER A PROJECTION
        HNEW=0.7D0*H
        REJECT=.TRUE.
        IF (ISTAT(2).GE.1) ISTAT(3)=ISTAT(3)+1
       END IF
       ELSE
C --- STEP IS REJECTED WITHOUT PROJECTION
        HNEW=H/DMIN1(FACC1,FAC11/SAFE)
        REJECT=.TRUE.  
        IF(ISTAT(2).GE.1)ISTAT(3)=ISTAT(3)+1
      END IF
      H = HNEW 
      LAST=.FALSE.
      GOTO 909
C --- FAIL EXIT
 176  CONTINUE
      WRITE(6,979)T   
      WRITE(6,*) ' MATRIX IN ASET IS  SINGULAR, IER=',IER
      IDID=-4
      RETURN
 177  CONTINUE
      WRITE(6,979)T   
      WRITE(6,*) ' STEP SIZE T0O SMALL, H=',H
      IDID=-3
      RETURN
 178  CONTINUE
      WRITE(6,979)T   
      WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED' 
      IDID=-2
      RETURN
C --- EXIT CAUSED BY SOLOUT
 179  CONTINUE
      WRITE(6,979)T
      WRITE(6,*) 'INITIAL PROJECTION: NO CONVERGENCE' 
      IDID=-5
      RETURN
 979  FORMAT(' EXIT OF PHEM56 AT X=',E18.4) 
      END


C******************************************************************
C --- FUNCTION FOR DENSE OUTPUT -- IT IS BASED ON HERMITE
C                INTERPOLATION
      FUNCTION POL4(I,FIRST,NQ,NV,NU,LRDO,X,DOWK)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL FIRST
      SAVE CP1,CP2,CP3,CP4
      DIMENSION DOWK(LRDO)
      XOLD = DOWK(1) 
      H= DOWK(2)
      S=(X-XOLD)/H
      NQVU=NQ+NV+NU
      NQVU2=2*NQVU
      NQVU3=3*NQVU
      NQVU4=4*NQVU
      IF (FIRST)  CALL COPOL4(S,H,CP1,CP2,CP3,CP4) 
      FIRST = .FALSE.
      IP2=I+2
      POL4=DOWK(IP2)+CP1*DOWK(IP2+NQVU)+CP2*DOWK(IP2+NQVU2)+
     &      CP3*DOWK(IP2+NQVU3)+CP4*DOWK(IP2+NQVU4)
      RETURN
      END
C******************************************************************
C --- COEFFICIENTS FOR HERMITE POLYNOMIAL
      SUBROUTINE COPOL4(X,H,CP1,CP2,CP3,CP4)
      IMPLICIT REAL*8 (A-H,O-Z)
C     T = 19.D0/30.D0
      T= 0.633333333333333333D0
      XM1=(X-1.D0)
      TMX=(T-X)
      FAC=X*XM1 
      DEN1=T
C      DEN2 = 1.D0-T 
      DEN2 = 0.3666666666666666667D0
C      DEN3 = (1.D0-T)**2*T**2
      DEN3 = 0.0539271604938271605D0
C      DEN4 = (1.D0-T)**2
      DEN4 = 0.1344444444444444444D0
C      P1 = T*(3.D0*T-4.D0)
      P1 = -1.33D0
C      P2 = 4.D0-2.D0*T**2
      P2 = 3.19777777777777777778D0
C      P3 = 2.D0*T-3.D0
      P3 = -1.7333333333333333333D0
      CP1 =H*FAC*XM1*TMX/DEN1 
      CP2 =-H*FAC*X*TMX/DEN2 
      CP3 =FAC**2/DEN3 
      CP4 =X**2*(P1+X*(P2+X*P3))/DEN4
      RETURN
      END
C ****************************************************************
C --- COEFFICIENTS OF THE METHOD 
       SUBROUTINE PCOEF56(
     &    C2,C3,C4,C5,C6,C7,C10,A21,A31,A32,
     &    A41,A42,A43,A51,A52,A53,A54,A61,A62,A63,A64,A65,
     &    A71,A72,A73,A74,A75,A76,BB1,BB3,BB4,BB5,BB6,BB7,
     &    B1,B2,B3,B4,B5,B6,B7,BC1,BC3,BC4,BC5,BC6,BC7,
     &    CC1,CC2,CC3,CC4,CC5,CC6,CC7,AA11,AA21,AA22,
     &    AA31,AA32,AA33,AA41,AA42,AA43,AA44,AA51,AA52,AA53,
     &    AA54,AA55, AA61,AA62,AA63,AA64,AA65,AA66,AA71,AA72,
     &    AA73,AA74,AA75,AA76,AA77)
       IMPLICIT REAL*8 (A-H,O-Z)
      write(*,*) 'PHERK method based on DOPRI5'
      c1 = 0.D0
      c2 =  1.D0/5.D0
      c3 =  3.D0/10.D0
      c4 =  4.D0/5.D0    
      c5 =  8.D0/9.D0    
      c6 = 1.D0          
      c7 = 1.D0
      b1 =  35.D0/384.D0
      b2 = 0.D0
       b3 =  500.D0/1113.D0
      b4 =  125.D0/192.D0
      b5 =  -2187.D0/6784.D0
      b6 =  11.D0/84.D0
       a21 =  1.D0/5.D0
      a31 =  3.D0/40.D0
      a32 =  9.D0/40.D0
      a41 =  44.D0/45.D0                   
       a42 =  -56.D0/15.D0                 
      a43 =  32.D0/9.D0                    
      a51 =  19372.D0/6561.D0              
       a52 =  -25360.D0/2187.D0            
      a53 =  64448.D0/6561.D0              
      a54 =  -212.D0/729.D0                
       a61 =  9017.D0/3168.D0              
      a62 =  -355.D0/33.D0                 
      a63 =  46732.D0/5247.D0              
       a64 =  49.D0/176.D0                 
      a65 =  -5103.D0/18656.D0             
      bb1 =  5179.D0/57600.D0 
      bb2 =  0.D0
      bb3 =  7571.D0/16695.D0
      bb4 =  393.D0/640.D0
      bb5 =  -92097.D0/339200.D0
      bb6 =  187.D0/2100.D0
      bb7 =  1.D0/40.D0
      bc1 =  1812923.D0/20736000.D0 
      bc2 =  0.D0
      bc3 =  45847.D0/100170.D0 
      bc4 =  56677.D0/414720.D0
      bc5 =  -289161.D0/13568000.D0
      bc6 =  -67507.D0/4536000.D0
      bc7 =  -3971.D0/324000.D0 
      cc1 =  0.D0
      cc2 =  3.D0/10.D0   
      cc3 =  3.D0/8.D0    
      cc4 =  5.D0/7.D0    
      cc5 =  1.D0         
      cc6 =  1.D0
      aa21 = 3.D0/40.D0                     
      aa22 = 9.D0/40.D0                     
      aa31 = 21.D0/256.D0                   
      aa32 = 45.D0/256.D0                   
      aa33 = 15.D0/128.D0                   
      aa41 = 19235.D0/98784.D0              
      aa42 = -225.D0/392.D0                 
      aa43 = 6235.D0/6174.D0                
      aa44 = 2755.D0/32928.D0               
      aa51 = -293.D0/1152.D0                
      aa52 = 85.D0/48.D0                    
      aa53 = -2275.D0/1908.D0               
      aa54 = 35.D0/32.D0                    
      aa55 = -2835.D0/6784.D0               
      aa61 =  b1
      aa62 =  b2
      aa63 =  b3
      aa64 =  b4
      aa65 =  b5
      aa66 =  b6
      cc7 =   19.D0/30.D0                      
      aa71 =  10752038939.D0/66307507200.D0    
      aa72 =  -19259.D0/50160.D0               
      aa73 =  19581770468.D0/24023520675.D0    
      aa74 =  9941999.D0/650073600.D0          
      aa75 =  180060831.D0/4820710400.D0       
      aa76 =  -11331659.D0/439538400.D0        
      aa77 =  1.D0/76.D0    
      return
      end
C*****************************************************************
      SUBROUTINE GIINUM(MODE,NQ,NV,NU,NL,NDIM,MDIM,NMDIM,NM, 
     &  LDG,LDF,LDA,NBLK,NMRC,NPGP,NPFL,INDG,INDGD,INDFL,T,
     &  UROUND,FPROB,Q,QD,QDOT,V,U,GQ,GQD,GT,GTD,
     &  FL,B,GII,RES)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FPROB
      DIMENSION INDG(2*LDG),INDGD(2*LDG),Q(NQ),QD(NQ)
      DIMENSION QDOT(NQ),V(NV),U(NU),GQ(LDG,NDIM)
      DIMENSION GQD(LDG,NDIM),GT(NL),GTD(NL),FL(LDF,MDIM)
      DIMENSION B(LDA,NMDIM),GII(NL),RES(NL),INDFL(2*LDF)
      DA=DSQRT(UROUND)
      CALL FPROB(3,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     NPGP,NPFL,INDGD(1),INDGD(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,Q,V,U,RES,
     &     RES,GQD,V,RES,GT,FL,QDOT,U,B)
      DO I=1,NQ
        QD(I)=Q(I)+DA*QDOT(I)
      END DO
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     NPGP,NPFL,INDGD(1),INDGD(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,QD,V,U,RES,
     &     RES,GQD,V,RES,GTD,FL,QDOT,U,B)
      CALL GPM2(MODE,NV,NL,NDIM,MDIM,LDG,INDG,INDGD,
     &     UROUND,GQ,GQD,V,RES)
      DO I=1,NL
        GII(I)=RES(I)+(GTD(I)-GT(I))/DA
      END DO
      CALL FPROB(6,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     NPGP,NPFL,INDGD(1),INDGD(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T+DA,Q,V,U,
     &     RES,RES,GQD,V,RES,GTD,FL,QDOT,U,B)
      CALL GPM2(MODE,NV,NL,NDIM,MDIM,LDG,INDG,INDGD,
     &     UROUND,GQ,GQD,V,RES)
      DO I=1,NL
        GII(I)=-(GII(I)+RES(I)+(GTD(I)-GT(I))/DA)
      END DO
      RETURN
      END
C**********************************************************
      SUBROUTINE GPM2(MODE,N,M,NDIM,MDIM,LDG,INDG,INDGD,
     &                UROUND,GQ,GQD,V,RES)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION INDG(2*LDG),INDGD(2*LDG),GQ(LDG,NDIM)
      DIMENSION GQD(LDG,NDIM),V(N),RES(M)
      DA=DSQRT(UROUND)
      IF ((MODE.GE.0).AND.(MODE.LE.3)) THEN
        GOTO 10
      ELSE IF ((MODE.EQ.4).OR.(MODE.EQ.5)) THEN
        GOTO 20
      ELSE
        STOP ' GPM2: INVALID MODE'
      ENDIF
 10   CONTINUE
      DO I=1,M
        SUM=0.D0
        DO J=1,N
          SUM=SUM+(GQD(I,J)-GQ(I,J))*V(J)/DA
        END DO
        RES(I)=SUM
      END DO
      RETURN
 20   CONTINUE
      DO I=1,M
        RES(I)=0.D0
      END DO
      DO K=1,LDG
        I=INDG(K)
        J=INDG(LDG+K)
        RES(I)=RES(I)+(GQD(K,1)-GQ(K,1))*V(J)/DA
      END DO
      RETURN
      END
C*******************************************************************
C*   SUBROUTINE ASET  CONSTRUCTS AND DECOMPOSES THE MATRIX     
C*       WITH :
C*                       |  AM    G0^t |                                   
C*                   B = |             |                                   
C*                       |  G1      0  |   
C*     
C*******************************************************************
      SUBROUTINE ASET(MODE,N,M,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &  NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG0,INDG1,
     &  INDFL,INDA,IP,IUMF,XUMF,B,G0,G1,FL,AVALUE,IER)
      IMPLICIT LOGICAL (A-Z)
      INTEGER MODE,N,M,NM,NDIM,MDIM,NMDIM,LDG,LDF,LDA,IP(NM),
     &  LIPAR,LIUMF,LXUMF,IPAR(LIPAR),IUMF(LIUMF),INDG0(2*LDG),
     &  INDG1(2*LDG),INDFL(2*LDF),NZA,INDA(2*NZA),IIRN,IICN,IIK,
     &  IREFAC,I,J,LISTA,ISTAT(LISTA),IER,LICN,IIW,ERROR,IICNM1
      DOUBLE PRECISION B(LDA,NMDIM),G0(LDG,NDIM),G1(LDG,NDIM),
     &  AVALUE(2*NZA),XUMF(LXUMF),FL(LDF,MDIM)
C    For a description of the internal parameters of the package
C    see: Harwell Subroutine Library Specification
C                  
C
      INTEGER LP, MP, IRNCP, ICNCP, MINIRN, MINICN, IRANK, IDISP(2)
      LOGICAL QBLOCK, QGROW, QABRT1, QABRT2
      DOUBLE PRECISION EPSMA, RMIN, RESID,THRSH1,THRSH2
C
      COMMON /MA28ED/ LP, MP, QBLOCK, QGROW
      COMMON /MA28FD/ EPSMA, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     $     IRANK, QABRT1, QABRT2
      COMMON /MA28GD/ IDISP
C
C                       Description see Harwell Subroutine
C                       Library Specification
C
      SAVE /MA28ED/, /MA28FD/, /MA28GD/
C
      IER=0
      IF (MODE.EQ.0) THEN
        GOTO 10
      ELSE IF (MODE.EQ.1) THEN
        GOTO 15
      ELSE IF (MODE.EQ.2) THEN
        GOTO 17
      ELSE IF (MODE.EQ.3) THEN
        GOTO 18
      ELSE IF ((MODE.EQ.4).OR.(MODE.EQ.5)) THEN
        GOTO 20
      ELSE
        STOP ' ASET: INVALID MODE'
      END IF
 10   CONTINUE
      DO  I=1,N
        DO  J=1,M
          B(I,N+J) = G0(J,I)
        END DO
      END DO
      DO  I=1,M
        DO  J=1,N
          B(N+I,J) = G1(I,J)
        END DO
        DO  J=1,M
          B(N+I,N+J) = 0.D0
        END DO
      END DO
      CALL DEC(NM,NM,B,IP,IER)
      ISTAT(7)=ISTAT(7)+1
      RETURN
 15   CONTINUE
      DO  I=1,N
        DO  J=1,M
          B(I,N+J) = FL(I,J)
        END DO
      END DO
      DO  I=1,M
        DO  J=1,N
          B(N+I,J) = G1(I,J)
        END DO
        DO  J=1,M
          B(N+I,N+J) = 0.D0
        END DO
      END DO
      CALL DEC(NM,NM,B,IP,IER)
      ISTAT(7)=ISTAT(7)+1
      RETURN
 17   CONTINUE
      DO  I=1,N
        DO  J=1,M
          B(I,N+J) = G0(J,I)
        END DO
      END DO
      DO  I=1,M
        DO  J=1,N
          B(N+I,J) = G1(I,J)
        END DO
        DO  J=1,M
          B(N+I,N+J) = 0.D0
        END DO
      END DO
      CALL DGETRF(NM,NM,B,NM,IP,IER)
      ISTAT(7)=ISTAT(7)+1
      RETURN
 18   CONTINUE
      DO  I=1,N
        DO  J=1,M
          B(I,N+J) = FL(I,J)
        END DO
      END DO
      DO  I=1,M
        DO  J=1,N
          B(N+I,J) = G1(I,J)
        END DO
        DO  J=1,M
          B(N+I,N+J) = 0.D0
        END DO
      END DO
      CALL DGETRF(NM,NM,B,NM,IP,IER)
      ISTAT(7)=ISTAT(7)+1
      RETURN
 20   CONTINUE
C
C
C
C  RAPPEL INDG(1->LDG)=LIGNE, INDG(LDG+1->2*LDG)=COLONNE
C
C
      THRSH1 = 1.D-2         
      THRSH2 = 1.D-6  
      LICN=2*NZA 
      IIRN=1
      IICN=IIRN+LICN
      IIK=IICN+LICN
      IIW=IIK+5*NM
      IICNM1=IICN-1
C ??      IFDEC = 0
      EPSMA =DMIN1 (1.D0, 10.D0*THRSH2)
C
      IF (IPAR(9).EQ.0) THEN
        IREFAC=IPAR(3)+IPAR(4)
      ELSE
        IPAR(9)=0
        IREFAC=1
      ENDIF
 100  CALL ACOPY(MODE,IPAR(1),IPAR(2),LDG,LDF,LDA,NZA,N,NMDIM,
     &     INDG0,INDG1,INDFL,INDA,AVALUE,B,G0,G1,FL)
       IF (IREFAC.GT.0) THEN
         WRITE(6,*)'MA28AD'
         DO I=1,NZA
          IUMF(I)=INDA(I)
          IUMF(IICNM1+I)=INDA(NZA+I)
         END DO
         CALL MA28AD(NM,NZA,AVALUE,LICN,IUMF(IIRN),LICN,IUMF(IICN),
     &        THRSH1,IUMF(IIK),IUMF(IIW),XUMF,ERROR)
         ISTAT(7)=ISTAT(7)+1
      ELSE
         CALL MA28BD(NM,NZA,AVALUE,LICN,INDA(1),INDA(NZA+1),IUMF(IICN),
     &           IUMF(IIK),IUMF(IIW),XUMF,ERROR )
         ISTAT(7)=ISTAT(7)+1
         IF(ERROR.LT.0) THEN
           WRITE(6,*)'ERROR IN UMDREFAC'
           IREFAC=1
           GOTO 100
         END IF
         IF(RMIN.LT.THRSH2) THEN
           IREFAC=1
           WRITE(6,*)'RMIN TOO SMALL'
           GOTO 100
         END IF
      END IF
      RETURN
      END  
C
C*******************************************************************
C*   SUBROUTINE GPMULT COMPUTES FAC*(AG*VECT+GT) = RES
C******************************************************************
      SUBROUTINE GPMULT(MODE,NDIM,N,M,LDG,INDG,
     &  AG,FAC,VECT,GT,RES)
      IMPLICIT LOGICAL (A-Z)
      INTEGER MODE,NDIM,N,M,LDG,INDG(2*LDG),
     &  I,J,K
      DOUBLE PRECISION AG(LDG,NDIM),VECT(N),GT(M),
     &  RES(M),FAC,SUM
      IF ((MODE.GE.0).AND.(MODE.LE.3)) THEN
        GOTO 10
      ELSE IF ((MODE.EQ.4).OR.(MODE.EQ.5)) THEN
        GOTO 20
      ELSE
        STOP ' GPMULT: INVALID MODE'
      END IF
 10   CONTINUE
      DO I=1,M
        SUM=0.D0
        DO J=1,N
          SUM=SUM+AG(I,J)*VECT(J)
        END DO
        RES(I)=FAC*(SUM+GT(I))
      END DO
      RETURN
 20   CONTINUE
      DO I=1,M
        RES(I)=FAC*GT(I)
      END DO
      DO K=1,LDG
        I=INDG(K)
        J=INDG(LDG+K)
        RES(I)=RES(I)+FAC*AG(K,1)*VECT(J)
      END DO
      RETURN
      END
C*******************************************************************
C*   SUBROUTINE AMULT COMPUTES FAC*AM*VECT = RES
C******************************************************************
      SUBROUTINE AMULT(MODE,N,NMDIM,LDA,NBLK,NMRC,FAC,AM,VECT,RES)
      IMPLICIT LOGICAL (A-Z)
      INTEGER MODE,N,NMDIM,NBLK,NMRC,I,J,K,KN,LDA
      DOUBLE PRECISION AM(LDA,NMDIM),VECT(N),FAC,SUM,RES(N)
      IF ((MODE.GE.0).AND.(MODE.LE.3)) THEN
        GOTO 10
      ELSE IF ((MODE.EQ.4).OR.(MODE.EQ.5)) THEN
        GOTO 20
      ELSE
        STOP ' AMULT: INVALID MODE'
      END IF
 10   CONTINUE
      DO I=1,N
        SUM=0.D0
        DO J=1,N
          SUM=SUM+AM(I,J)*VECT(J)
        END DO
        RES(I)=FAC*SUM
      END DO
      RETURN
 20   CONTINUE
      DO K=1,NBLK
        KN=(K-1)*NMRC
        DO I=1,NMRC
          SUM=0.D0
          DO J=1,NMRC
            SUM=SUM+AM(KN+I,J)*VECT(KN+J)
          END DO
          RES(KN+I)=FAC*SUM
        END DO
      END DO
      RETURN
      END
C*******************************************************************
C*    SUBROUTINE ACOPY
C*******************************************************************
      SUBROUTINE ACOPY(MODE,NMRC,NBLK,LDG,LDF,LDA,NZA,N,NMDIM,
     & INDG1,INDG2,INDFL,INDA,AVALUE,AM,G1,G2,FL)
      IMPLICIT LOGICAL (A-Z)
      INTEGER MODE,NMRC,NBLK,LDG,LDF,LDA,NZA,N,NMDIM,NMNB,I,J,K,
     &  L,KN,INDG1(2*LDG),INDG2(2*LDG),INDA(2*NZA),INDFL(2*LDF)
      DOUBLE PRECISION AVALUE(2*NZA),G1(LDG),G2(LDG),
     &  AM(LDA,NMDIM),FL(LDF)
      IF (MODE.EQ.4) THEN
        GOTO 10
      ELSE IF (MODE.EQ.5) THEN
        GOTO 20
      ELSE
        STOP ' ACOPY: INVALID MODE'
      END IF
 10   CONTINUE
      NMNB=NMRC*NBLK
      L=0
      DO K=1,NBLK
        KN=(K-1)*NMRC
        DO i=1,NMRC
          DO j=1,NMRC
            L=L+1
            AVALUE(L)=AM(I+KN,J)
            INDA(L)=KN+I
            INDA(NZA+L)=KN+J
          END DO
        END DO
      END DO
      DO I=1,LDG
        L=L+1
        AVALUE(L)=G1(I)
        INDA(L)=INDG1(LDG+I)
        INDA(NZA+L)=NMNB+INDG1(I)	
        L=L+1
        AVALUE(L)=G2(I)
        INDA(L)=NMNB+INDG2(I)
        INDA(NZA+L)=INDG2(LDG+I)
      END DO
      RETURN
 20   CONTINUE
      NMNB=NMRC*NBLK
      L=0
      DO K=1,NBLK
        KN=(K-1)*NMRC
        DO i=1,NMRC
          DO j=1,NMRC
            L=L+1
            AVALUE(L)=AM(I+KN,J)
            INDA(L)=KN+I
            INDA(NZA+L)=KN+J
          END DO
        END DO
      END DO
      DO I=1,LDF
        L=L+1
        AVALUE(L)=FL(I)
        INDA(L)=INDFL(I)
        INDA(NZA+L)=NMNB+INDFL(LDF+I)
      END DO
      DO I=1,LDG	
        L=L+1
        AVALUE(L)=G2(I)
        INDA(L)=NMNB+INDG2(I)
        INDA(NZA+L)=INDG2(LDG+I)
      END DO
      RETURN
      END	
C*******************************************************************
C*   SUBROUTINE ASOL SOLVES THE SYSTEM B*X0=R
C*       WITH THE MATRIX B DECOMPOSED IN ASET
C*       AND MOVES THE RESULT X0(I=1,..,N)   -> RES(I)                      
C*******************************************************************
      SUBROUTINE ASOL(MODE,NM,N,M,LIPAR,LIUMF,LXUMF,IPAR,
     &    IUMF,IP,NZA,AVALUE,XUMF,B,X0,RES1,RES2)
      IMPLICIT LOGICAL (A-Z)
      INTEGER MODE,N,NM,NZA,LIPAR,LIUMF,LXUMF,IPAR(LIPAR),
     &  IUMF(LIUMF),IP(NM),IIW,I,M,LICN,IICN,IIK,MTYPE
      DOUBLE PRECISION B(NM,NM),X0(NM),RES1(N),RES2(M),XUMF(LXUMF),
     &  AVALUE(2*NZA)
C -- umfpack variables:
      INTEGER MTYPE
C
      IF ((MODE.EQ.0).OR.(MODE.EQ.1)) THEN
        GOTO 10
      ELSE IF ((MODE.EQ.2).OR.(MODE.EQ.3)) THEN
        GOTO 15
      ELSE IF ((MODE.EQ.4).OR.(MODE.EQ.5)) THEN
        GOTO 20
      ELSE 
        STOP ' ASOL: INVALID MODE'
      END IF
 10   CONTINUE
      CALL SOL(NM,NM,B,X0,IP)
      GOTO 100
 15   CONTINUE
      CALL DGETRS('N',NM,1,B,NM,IP,X0,NM,IER)
      GOTO 100
 20   CONTINUE
      LICN=2*NZA
      IICN=1+LICN
      IIK=IICN+LICN
      IIW=IIK+5*NM
      MTYPE=1
      CALL MA28CD(NM,AVALUE,LICN,IUMF(IICN),IUMF(IIK),X0,
     &   XUMF,MTYPE)
 100   CONTINUE
      DO  I=1,N
        RES1(I) = X0(I)
      END DO
      DO  I=1,M
        RES2(I) = X0(N+I)
      END DO
      RETURN
      END	
C***********************************************************
C*  SUBROUTINE APROJ  PROJECTION ON THE MANIFOLD DEFINED
C*                     BY {g(q) = 0} 
C*               
C*     SOLVES BY NEWTON ITERATIONS A NON LINEAR SYSTEM
C***********************************************************
      SUBROUTINE APROJ(MODE,NQ,NV,NU,NL,NM,NDIM,MDIM,NMDIM,
     *   NBLK,NMRC,LDG,LDF,LDA,NZA,LIPAR,LIUMF,LXUMF,LISTA,
     *   ISTAT,IPAR,IUMF,INDA,INDG,INDFL,FL,XUMF,AVALUE,T,
     *   FPROB,Q,Q1,Q2,QDOT,V,V1,V2,U,XL,UDOT, G,GT,GQ,AM,
     *   B,X0,IP,ATOL,RTOL,ITOL,ACCEPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(NQ),Q1(NQ),Q2(NQ),QDOT(NQ),V(NV),V1(NV)
      DIMENSION AM(LDA,NMDIM),G(NL),B(LDA,NMDIM),V2(NV),U(NU)
      DIMENSION X0(NM),GQ(LDG,NDIM),GT(NL),UDOT(NU),XL(NL)
      DIMENSION IP(NM),ATOL(1),RTOL(1),INDG(2*LDG),ISTAT(LISTA)
      DIMENSION INDFL(2*LDF),FL(LDF,MDIM),IPAR(LIPAR),IUMF(LIUMF)
      DIMENSION INDA(2*NZA),XUMF(LXUMF),AVALUE(2*NZA)
      EXTERNAL FPROB
      LOGICAL ACCEPT
      CALL FPROB(9,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,Q,V,U,XL,
     &     G,GQ,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,AM)
C --- INITIAL PREPARATION  ------------------------------------
      ERROLD=1.D20 
      FAC = -1.D0
      DO I=1,NV
        Q2(I)=Q(I)
        V2(I)=0.D0
      END DO
      DO I=NV+1,NQ
        Q2(I)=Q(I)
      END DO
      K=0
      ITMAX = 5
      EPS = 1.d-2
C --- COMPUTE DEFECT IN  G    --------------------------------------
      CALL FPROB(4,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,Q,V,U,XL,
     &     G,GQ,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
      DO  I=1,NV
        X0(I)=0.D0
      END DO
      DO  I=1,NL
        X0(NV+I)=-G(I)
      END DO 
 106  CALL  ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,IP,
     &      NZA,AVALUE,XUMF,B,X0,X0,XL)
      ISTAT(8)=ISTAT(8)+1
      DO I=1,NV
        V2(I)=V2(I)+X0(I)
      END DO
      CALL FPROB(2,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,Q,X0,U,XL,
     &     G,GQ,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
      DO I=1,NQ
        Q2(I)=Q2(I)+QDOT(I)
      END DO
      ERR=0.d0
      DO  I=1,NQ      
      IF (ITOL.EQ.0) THEN
        SK=ATOL(1)+RTOL(1)*MAX(ABS(Q(I)),ABS(Q1(I)))
      ELSE
        SK=ATOL(I)+RTOL(I)*MAX(ABS(Q(I)),ABS(Q1(I)))
      ENDIF
      ERR = ERR+(QDOT(I)/SK)**2
      END DO
      ERR=SQRT(ERR/NQ)
C
C ---    TEST FOR CONVERGENCE    --------------------------------------
C
      IF(ERR.LT.EPS) GOTO 200
      IF(ERR.GE.ERROLD*2.D0 .OR. K.GE.ITMAX) THEN
        WRITE(6,*)'NO CONVERGENCE IN PROJECTION'
C ---    NO CONVERGENCE STEP IS REJECTED   --------------------------------------------
        ACCEPT = .FALSE.
        RETURN
      ENDIF
C
C ---    NEXT ITERATION    --------------------------------------------
C
      ERROLD=ERR
      K=K+1 
      CALL FPROB(4,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,Q2,V,U,XL,
     &     G,GQ,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
      CALL AMULT(MODE,NV,NMDIM,LDA,NBLK,NMRC,FAC,AM,V2,X0)
      DO  I=1,NL
        X0(NV+I)=-G(I)
      END DO
      GOTO 106
 200  CONTINUE
C
c ---    PROJECTION OF VELOCITY     -----------------------------------
c
cc
cc
      DO  I=1,NQ
        Q(I)=Q2(I)
      END DO
      CALL FPROB(11,NQ,NV,NU,NL,LDG,LDF,LDA,NBLK,NMRC,
     &     IPAR(3),IPAR(4),INDG(1),INDG(LDG+1),INDFL(1),
     &     INDFL(LDF+1),T,Q,V,U,XL,
     &     G,GQ,X0(1),X0(NV+1),GT,FL,QDOT,UDOT,B)
      CALL ASET(MODE,NV,NL,NM,NMDIM,NDIM,MDIM,LDG,LDF,LDA,
     &     NZA,LIPAR,LIUMF,LXUMF,LISTA,ISTAT,IPAR,INDG,INDG,
     &     INDFL,INDA,IP,IUMF,XUMF,B,GQ,GQ,FL,AVALUE,IER)
      IF (IER.NE.0) THEN
        WRITE(6,*)'MATRIX SINGULAR IN APROJ'
        ACCEPT=.FALSE.
        RETURN
      END IF
      DO I=1,NV
        X0(I)=0.D0
      END DO
      FAC=-1.D0
      CALL GPMULT(MODE,NDIM,NV,NL,LDG,INDG,
     &     GQ,FAC,V,GT,X0(NV+1))
      CALL ASOL(MODE,NM,NV,NL,LIPAR,LIUMF,LXUMF,IPAR,IUMF,IP,
     &     NZA,AVALUE,XUMF,B,X0,X0,XL)
      ISTAT(8)=ISTAT(8)+1
      DO I=1,NV
         V(I)=V(I)+X0(I)
      END DO
C -- ERROR ESTIMATION 
      ERR=0.d0
      DO I=1,NV
        IF (ITOL.EQ.0) THEN
          SK=ATOL(1)+RTOL(1)*MAX(ABS(V(I)),ABS(V1(I)))
        ELSE
          SK=ATOL(I+NQ)+RTOL(I+NQ)*MAX(ABS(V(I)),ABS(V1(I)))
        ENDIF
        ERR = ERR+(X0(I)/SK)**2
      END DO
      ERR=SQRT(ERR/NV)
      IF(ERR.GT.EPS)ACCEPT=.FALSE.
      IF(ERR.GT.EPS) write(6,*) 'rejected after proj',err
      RETURN
      END
