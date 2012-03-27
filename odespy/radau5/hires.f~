C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 ON HIRES PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=8,LWORK=4*ND*ND+12*ND+20,LIWORK=3*ND+20)
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),ISTAT(20)
C -------- END PARAMETER LIST --------
        REAL*4 TARRAY(2)
        REAL*4 TIMERES
        EXTERNAL FHIRES,JHIRES,SOLOUT
c ------ FILE DE DONNEES ----------
        OPEN(8,FILE='res_rad5')
        REWIND 8
C --- LOOP FOR DIFFERENT TOLERANCES
        NTOLMN=2
        NTOLMX=10
        NTOLDF=4
        NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
        TOLST=0.1D0**NTOLMN
        TOLFC=0.1D0**(1.D0/NTOLDF)
        DO 30 NTOL=1,NRLOOP
C --- DIMENSION OF THE SYSTEM
        N=8
C --- COMPUTE THE JACOBIAN
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS NOT USED DURING INTEGRATION
        IOUT=0
C --- INITIAL VALUES
        X=0.0D0
      Y (1) = 1.D0
      Y (2) = 0.D0
      Y (3) = 0.D0
      Y (4) = 0.D0
      Y (5) = 0.D0
      Y (6) = 0.D0
      Y (7) = 0.D0
      Y (8) = 0.0057D0
C --- REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL*1.D-4
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
        XEND=321.8122D0
        CALL DTIME(TARRAY, TIMERES)
        DO 20 I=1,2
C --- CALL OF THE SUBROUTINE RADAU5
        CALL RADAU5(N,FHIRES,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JHIRES,IJAC,MLJAC,MUJAC,
     &                  FHIRES,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT SOLUTION
      DO 15 K=1,8
 15   WRITE (8,*)Y(K)
C --- PRINT STATISTICS
         DO J=14,20
            ISTAT(J)=ISTAT(J)+IWORK(J)
         END DO
 20     XEND=XEND+100.D0
        CALL DTIME(TARRAY, TIMERES)
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

C
      SUBROUTINE FHIRES (N, X, Y, DY, RPAR, IPAR)
      INTEGER N
      DOUBLE PRECISION X, Y, DY
      DIMENSION Y (8), DY (8)
      DY (1) = -1.71D0*Y(1) + 0.43D0*Y(2) + 8.32D0*Y(3) + 0.0007D0
      DY (2) = 1.71D0*Y(1) - 8.75D0*Y(2)
      DY (3) = -10.03D0*Y(3) + 0.43D0*Y(4) + 0.035D0*Y(5)
      DY (4) = 8.32D0*Y(2) + 1.71D0*Y(3) - 1.12D0*Y(4)
      DY (5) = -1.745D0*Y(5) + 0.43D0*Y(6) + 0.43D0*Y(7)
      DY (6) = -280.0D0*Y(6)*Y(8) + 0.69D0*Y(4) + 1.71D0*Y(5) -
     &         0.43D0*Y(6) + 0.69D0*Y(7)
      DY (7) = 280.0D0*Y(6)*Y(8) - 1.81D0*Y(7)
      DY (8) = -DY (7)
      RETURN 
      END


      SUBROUTINE JHIRES(N,X,Y,DFY,LDFY,RPAR,IPAR)
C -------- JACOBIAN FOR HIRES PROBLEM --------
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
C------ METTRE A ZERO -------
      DO 1 I=1,N
      DO 1 J=1,N
  1   DFY(I,J)=0.D0
C
      DFY(1,1)=  -1.71D0
      DFY(1,2)=   0.43D0
      DFY(1,3)=   + 8.32D0
C
      DFY(2,1)=   1.71D0
      DFY(2,2)=   - 8.75D0
C
      DFY(3,3)=   -10.03D0
      DFY(3,4)=   0.43D0
      DFY(3,5)=   + 0.035D0
C
      DFY(4,2)=    8.32D0
      DFY(4,3)=   + 1.71D0
      DFY(4,4)=   - 1.12D0
C
      DFY(5,5)=   -1.745D0
      DFY(5,6)=   + 0.43D0
      DFY(5,7)=   + 0.43D0
C
      DFY(6,4)=    + 0.69D0
      DFY(6,5)=    + 1.71D0
      DFY(6,6)=    - 0.43D0  -280.0D0*Y(8)
      DFY(6,7)=    + 0.69D0
      DFY(6,8)=     -280.0D0*Y(6)
C 
      DFY(7,6)=     280.0D0*Y(8)
      DFY(7,7)=      - 1.81D0
      DFY(7,8)=     280.0D0*Y(6)
C
      DFY(8,6)=    - 280.0D0*Y(8)
      DFY(8,7)=       1.81D0
      DFY(8,8)=    - 280.0D0*Y(6)
      RETURN
      END


