C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON VDPOL PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=2,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2)
        EXTERNAL FVANDER,JVANDER,SOLOUT
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_rad5')
        REWIND 8
        RPAR=1.D-6
C --- LOOP FOR DIFFERENT TOLERANCES
        NTOLMN=2
        NTOLMX=10
        NTOLDF=4
        NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
        TOLST=0.1D0**NTOLMN
        TOLFC=0.1D0**(1.D0/NTOLDF)
        DO 30 NTOL=1,NRLOOP
C --- DIMENSION OF THE SYSTEM
        N=2
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=2.0D0
        Y(2)=0.0D0
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6 
C --- SET DEFAULT VALUES 
        DO I=1,20
           WORK(I)=0.D0
           IWORK(I)=0
           ISTAT(I)=0
        END DO
C --- ENDPOINT OF INTEGRATION
        XEND=1.0D0
        CALL DTIME(TARRAY)
        DO 20 I=1,11
C --- CALL OF THE SUBROUTINE RADAU5 (OR SDIRK4)
        CALL RADAU5(N,FVANDER,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JVANDER,IJAC,MLJAC,MUJAC,
     &                  FVANDER,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
        WRITE (8,*) Y(1)
        WRITE (8,*) Y(2)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
 20     XEND=XEND+1.D0
        CALL DTIME(TARRAY)
        WRITE(8,*)TARRAY(1)
        WRITE (8,*)(ISTAT(J),J=14,20)
        WRITE(6,*)' ***** TOL=',RTOL,'  ELAPSED TIME=',TARRAY(1),' ****'
        WRITE (6,91) (ISTAT(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
C -------- NEW TOLERANCE ---
        TOLST=TOLST*TOLFC
 30     CONTINUE
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
           WRITE (6,99) X,Y(1),Y(2),NR-1
 99     FORMAT(1X,'X =',F6.2,'    Y =',2F18.7,'    NSTEP =',I4)
        RETURN
        END

        SUBROUTINE FVANDER(N,X,Y,F,EPS,IPAR)
C --- RIGHT-HAND SIDE OF VANDERPOL EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=Y(2)
        PROD=1.D0-Y(1)*Y(1)
        F(2)=(PROD*Y(2)-Y(1))/EPS
        RETURN
        END 
C
        SUBROUTINE JVANDER(N,X,Y,DFY,LDFY,EPS,IPAR)
C --- JACOBIAN OF VANDERPOL EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        DFY(1,1)=0.D0
        DFY(1,2)=1.D0
        DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/EPS
        DFY(2,2)=(1.0D0-Y(1)**2)/EPS
        RETURN
        END

