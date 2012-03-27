c-----------------------------------------------------------------------
c     additional linear algebra routines required by RADAU5
c-----------------------------------------------------------------------
C ******************************************
C     VERSION OF SEPTEMBER 18, 1995
C ******************************************
C
      SUBROUTINE DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          IP1(NM1),IPHES(N)
      LOGICAL CALHES
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO  I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
 45   MM=M1/M2
      DO J=1,M2
         DO I=1,NM1
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I,J)=E1(I,J)-SUM
         END DO
      END DO
      CALL DEC (NM1,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
  46  MM=M1/M2
      DO J=1,M2
         DO I=1,MBJAC
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I+MLE,J)=E1(I+MLE,J)-SUM
         END DO
      END DO
      CALL DECB (NM1,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,JM1)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      IF (CALHES) CALL ELMHES (LDJAC,N,1,N,FJAC,IPHES)
      CALHES=.FALSE.
      DO J=1,N-1
         J1=J+1
         E1(J1,J)=-FJAC(J1,J)
      END DO
      DO J=1,N
         DO I=1,J
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DECH(N,LDE1,E1,1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMR
C
C ***********************************************************
C
      SUBROUTINE DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),
     &          E2R(LDE1,NM1),E2I(LDE1,NM1),IP2(NM1)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECC (N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
  45  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,NM1
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            E2R(I,J)=E2R(I,J)-SUMR
            E2I(I,J)=E2I(I,J)-SUMI
         END DO
      END DO
      CALL DECC (NM1,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=BETAN
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=E2I(MDIAG,J)+BETAN
      END DO
  46  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,MBJAC
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            IMLE=I+MLE
            E2R(IMLE,J)=E2R(IMLE,J)-SUMR
            E2I(IMLE,J)=E2I(IMLE,J)-SUMI
         END DO
      END DO
      CALL DECBC (NM1,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO  J=1,N
         DO  I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
      END DO
      DO J=1,N
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            BB=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*BB
            E2I(I,J)=BETAN*BB
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            FFMA=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*FFMA
            E2I(I,J)=E2I(I,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         DO I=MAX(1,MUMAS+2-J),MIN(MBB,MUMAS+1-J+N)
            IB=I+MDIFF
            BB=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*BB
            E2I(IB,J)=BETAN*BB
         END DO
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            FFMA=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*FFMA
            E2I(IB,J)=E2I(IB,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            BB=FMAS(I,J)
            E2R(I,J)=BB*ALPHN-FJAC(I,J)
            E2I(I,J)=BB*BETAN
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=ALPHN*FMAS(I,J)-FJAC(I,JM1)
            E2I(I,J)=BETAN*FMAS(I,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO J=1,N-1
         J1=J+1
         E2R(J1,J)=-FJAC(J1,J)
         E2I(J1,J)=0.D0
      END DO
      DO J=1,N
         DO I=1,J
            E2I(I,J)=0.D0
            E2R(I,J)=-FJAC(I,J)
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECHC(N,LDE1,E2R,E2I,1,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMC
C
C ***********************************************************
C
      SUBROUTINE SLVRAR(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E1,LDE1,Z1,F1,IP1,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          IP1(NM1),IPHES(N),Z1(N),F1(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
 48   CONTINUE
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1)
 49   CONTINUE
      DO I=M1,1,-1
         Z1(I)=(Z1(I)+Z1(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
  45  CONTINUE
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=1,N
            S1=S1-FMAS(I,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         DO J=1,NM1
            S1=S1-FMAS(I,J)*F1(J+M1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N
             Z1(I)=Z1(I)-FJAC(I,MP1)*Z1(MP)
          END DO
       END DO
       CALL SOLH(N,LDE1,E1,1,Z1,IP1)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N
             Z1(I)=Z1(I)+FJAC(I,MP1)*Z1(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAR
C
C ***********************************************************
C
      SUBROUTINE SLVRAI(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,Z2,Z3,
     &          F2,F3,CONT,IP2,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),
     &          IP2(NM1),IPHES(N),Z2(N),Z3(N),F2(N),F3(N)
      DIMENSION E2R(LDE1,NM1),E2I(LDE1,NM1)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
 48   ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=1,NM1
               IM1=I+M1
               Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2)
 49   CONTINUE
      DO I=M1,1,-1
         MPI=M2+I
         Z2I=Z2(I)+Z2(MPI)
         Z3I=Z3(I)+Z3(MPI)
         Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO
         Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
  45  ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               IIMU=I+MUJAC+1-J
               Z2(IM1)=Z2(IM1)+FJAC(IIMU,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(IIMU,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            JM1=J+M1
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=1,N
            BB=FMAS(I,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S2=0.0D0
         S3=0.0D0
         DO J=1,NM1
            JM1=J+M1
            BB=FMAS(I,J)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z2(I)=Z2(I)-E1IMP*Z2(MP)
             Z3(I)=Z3(I)-E1IMP*Z3(MP)
          END DO
       END DO
       CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z2(I)=Z2(I)+E1IMP*Z2(MP)
             Z3(I)=Z3(I)+E1IMP*Z3(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAI
C
C ***********************************************************
C
      SUBROUTINE SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3,
     &          F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          E2R(LDE1,NM1),E2I(LDE1,NM1),IP1(NM1),IP2(NM1),
     &          IPHES(N),Z1(N),Z2(N),Z3(N),F1(N),F2(N),F3(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
 48   ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=1,NM1
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1
               Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1)
      CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2)
 49   CONTINUE
      DO I=M1,1,-1
         MPI=M2+I
         Z1(I)=(Z1(I)+Z1(MPI))/FAC1
         Z2I=Z2(I)+Z2(MPI)
         Z3I=Z3(I)+Z3(MPI)
         Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO
         Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
  45  ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               FFJA=FJAC(I+MUJAC+1-J,JKM)
               Z1(IM1)=Z1(IM1)+FFJA*SUM1
               Z2(IM1)=Z2(IM1)+FFJA*SUM2
               Z3(IM1)=Z3(IM1)+FFJA*SUM3
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1)
      CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         J1B=MAX(1,I-MLMAS)
         J2B=MIN(NM1,I+MUMAS)
         DO J=J1B,J2B
            JM1=J+M1
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,N
            BB=FMAS(I,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,NM1
            JM1=J+M1
            BB=FMAS(I,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)-E1IMP*Z1(MP)
             Z2(I)=Z2(I)-E1IMP*Z2(MP)
             Z3(I)=Z3(I)-E1IMP*Z3(MP)
          END DO
       END DO
       CALL SOLH(N,LDE1,E1,1,Z1,IP1)
       CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)+E1IMP*Z1(MP)
             Z2(I)=Z2(I)+E1IMP*Z2(MP)
             Z3(I)=Z3(I)+E1IMP*Z3(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAD
C
C ***********************************************************
C
      SUBROUTINE ESTRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD1,DD2,DD3,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,
     &          E1,LDE1,Z1,Z2,Z3,CONT,F1,F2,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),IP1(NM1),
     &     SCAL(N),IPHES(N),Z1(N),Z2(N),Z3(N),F1(N),F2(N),Y0(N),Y(N)
      DIMENSION CONT(N),RPAR(1),IPAR(1)
      LOGICAL FIRST,REJECT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      HEE1=DD1/H
      HEE2=DD2/H
      HEE3=DD3/H
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
   1  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO  I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  11  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
  48  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   2  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  12  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
  45  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   3  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  13  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   4  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  14  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 45
C
   5  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  15  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=1,NM1
            SUM=SUM+FMAS(I,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   6  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ------  THIS OPTION IS NOT PROVIDED
      RETURN
C
   7  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO MM=N-2,1,-1
         MP=N-MM
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 310
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 310     CONTINUE
         DO I=MP+1,N
            CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
         END DO
      END DO
      CALL SOLH(N,LDE1,E1,1,CONT,IP1)
      DO MM=1,N-2
         MP=N-MM
         DO I=MP+1,N
            CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 440
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 440     CONTINUE
      END DO
C
C --------------------------------------
C
  77  CONTINUE
      ERR=0.D0
      DO  I=1,N
         ERR=ERR+(CONT(I)/SCAL(I))**2
      END DO
      ERR=MAX(SQRT(ERR/N),1.D-10)
C
      IF (ERR.LT.1.D0) RETURN
      IF (FIRST.OR.REJECT) THEN
          DO I=1,N
             CONT(I)=Y(I)+CONT(I)
          END DO
          CALL FCN(N,X,CONT,F1,RPAR,IPAR)
          NFCN=NFCN+1
          DO I=1,N
             CONT(I)=F1(I)+F2(I)
          END DO
          GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB
C ------ FULL MATRIX OPTION
  31      CONTINUE
          CALL SOL(N,LDE1,E1,CONT,IP1)
          GOTO 88
C ------ FULL MATRIX OPTION, SECOND ORDER
 41      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=1,NM1
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOL(NM1,LDE1,E1,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
         GOTO 88
C ------ BANDED MATRIX OPTION
 32      CONTINUE
         CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
         GOTO 88
C ------ BANDED MATRIX OPTION, SECOND ORDER
 42      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ HESSENBERG MATRIX OPTION
  33      CONTINUE
          DO MM=N-2,1,-1
             MP=N-MM
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 510
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 510         CONTINUE
             DO I=MP+1,N
                CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
             END DO
          END DO
          CALL SOLH(N,LDE1,E1,1,CONT,IP1)
          DO MM=1,N-2
             MP=N-MM
             DO I=MP+1,N
                CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
             END DO
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 640
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 640         CONTINUE
          END DO
C -----------------------------------
   88     CONTINUE
          ERR=0.D0
          DO I=1,N
             ERR=ERR+(CONT(I)/SCAL(I))**2
          END DO
          ERR=MAX(SQRT(ERR/N),1.D-10)
       END IF
       RETURN
C -----------------------------------------------------------
  55   CONTINUE
       RETURN
       END
C
C     END OF SUBROUTINE ESTRAD
C
C ***********************************************************
C
      SUBROUTINE ESTRAV(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,NS,NNS,
     &          E1,LDE1,ZZ,CONT,FF,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),IP1(NM1),
     &     SCAL(N),IPHES(N),ZZ(NNS),FF(NNS),Y0(N),Y(N)
      DIMENSION DD(NS),CONT(N),RPAR(1),IPAR(1)
      LOGICAL FIRST,REJECT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
   1  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  11  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
  48  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   2  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  12  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
  45  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   3  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  13  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   4  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  14  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 45
C
   5  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  15  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=1,NM1
            SUM=SUM+FMAS(I,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   6  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ------  THIS OPTION IS NOT PROVIDED
      RETURN
C
   7  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO MM=N-2,1,-1
         MP=N-MM
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 310
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 310     CONTINUE
         DO I=MP+1,N
            CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
         END DO
      END DO
      CALL SOLH(N,LDE1,E1,1,CONT,IP1)
      DO MM=1,N-2
         MP=N-MM
         DO I=MP+1,N
            CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 440
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 440     CONTINUE
      END DO
C
C --------------------------------------
C
  77  CONTINUE
      ERR=0.D0
      DO  I=1,N
         ERR=ERR+(CONT(I)/SCAL(I))**2
      END DO
      ERR=MAX(SQRT(ERR/N),1.D-10)
C
      IF (ERR.LT.1.D0) RETURN
      IF (FIRST.OR.REJECT) THEN
          DO I=1,N
             CONT(I)=Y(I)+CONT(I)
          END DO
          CALL FCN(N,X,CONT,FF,RPAR,IPAR)
          NFCN=NFCN+1
          DO I=1,N
             CONT(I)=FF(I)+FF(I+N)
          END DO
          GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB
C ------ FULL MATRIX OPTION
 31      CONTINUE
         CALL SOL (N,LDE1,E1,CONT,IP1)
          GOTO 88
C ------ FULL MATRIX OPTION, SECOND ORDER
 41      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=1,NM1
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ BANDED MATRIX OPTION
 32      CONTINUE
         CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
          GOTO 88
C ------ BANDED MATRIX OPTION, SECOND ORDER
 42      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ HESSENBERG MATRIX OPTION
  33      CONTINUE
          DO MM=N-2,1,-1
             MP=N-MM
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 510
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 510         CONTINUE
             DO I=MP+1,N
                CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
             END DO
          END DO
          CALL SOLH(N,LDE1,E1,1,CONT,IP1)
          DO MM=1,N-2
             MP=N-MM
             DO I=MP+1,N
                CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
             END DO
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 640
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 640         CONTINUE
          END DO
C -----------------------------------
  88      CONTINUE
          ERR=0.D0
          DO I=1,N
             ERR=ERR+(CONT(I)/SCAL(I))**2
          END DO
          ERR=MAX(SQRT(ERR/N),1.D-10)
       END IF
       RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
       END
C
C     END OF SUBROUTINE ESTRAV
C
C ***********************************************************
C
      SUBROUTINE SLVROD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E,LDE,IP,DY,AK,FX,YNEW,HD,IJOB,STAGE1)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E(LDE,NM1),
     &          IP(NM1),DY(N),AK(N),FX(N),YNEW(N)
      LOGICAL STAGE1
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      IF (HD.EQ.0.D0) THEN
         DO  I=1,N
           AK(I)=DY(I)
         END DO
      ELSE
         DO I=1,N
            AK(I)=DY(I)+HD*FX(I)
         END DO
      END IF
C
      GOTO (1,2,3,4,5,6,55,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
 48   MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(AK(JKM)+SUM)/FAC1
            DO I=1,NM1
               IM1=I+M1
               AK(IM1)=AK(IM1)+FJAC(I,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE,E,AK(M1+1),IP)
      DO I=M1,1,-1
         AK(I)=(AK(I)+AK(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
  45  MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(AK(JKM)+SUM)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               AK(IM1)=AK(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE,E,MLE,MUE,AK(M1+1),IP)
      DO I=M1,1,-1
         AK(I)=(AK(I)+AK(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
      DO  I=1,N
         SUM=0.D0
         DO  J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,M1
            AK(I)=AK(I)+YNEW(I)
         END DO
         DO I=1,NM1
            SUM=0.D0
            DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
                SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J+M1)
            END DO
            IM1=I+M1
            AK(IM1)=AK(IM1)+SUM
         END DO
      END IF
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      IF (STAGE1) THEN
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,M1
            AK(I)=AK(I)+YNEW(I)
         END DO
         DO I=1,NM1
            SUM=0.D0
            DO J=1,NM1
               SUM=SUM+FMAS(I,J)*YNEW(J+M1)
            END DO
            IM1=I+M1
            AK(IM1)=AK(IM1)+SUM
         END DO
      END IF
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      IF (STAGE1) THEN
      DO 624 I=1,N
         SUM=0.D0
         DO 623 J=1,N
  623       SUM=SUM+FMAS(I,J)*YNEW(J)
  624    AK(I)=AK(I)+SUM
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      END IF
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVROD
C
C
C ***********************************************************
C
      SUBROUTINE SLVSEU(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E,LDE,IP,IPHES,DEL,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E(LDE,NM1),DEL(N)
      DIMENSION IP(NM1),IPHES(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,1,2,1,55,7,55,55,55,11,12,11,12,11), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      CALL SOL (N,LDE,E,DEL,IP)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(DEL(JKM)+SUM)/FAC1
            DO I=1,NM1
               IM1=I+M1
               DEL(IM1)=DEL(IM1)+FJAC(I,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE,E,DEL(M1+1),IP)
      DO I=M1,1,-1
         DEL(I)=(DEL(I)+DEL(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      CALL SOLB (N,LDE,E,MLE,MUE,DEL,IP)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(DEL(JKM)+SUM)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               DEL(IM1)=DEL(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE,E,MLE,MUE,DEL(M1+1),IP)
      DO I=M1,1,-1
         DEL(I)=(DEL(I)+DEL(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  HESSENBERG OPTION
      DO MMM=N-2,1,-1
         MP=N-MMM
         MP1=MP-1
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 110
         ZSAFE=DEL(MP)
         DEL(MP)=DEL(I)
         DEL(I)=ZSAFE
 110     CONTINUE
         DO I=MP+1,N
            DEL(I)=DEL(I)-FJAC(I,MP1)*DEL(MP)
         END DO
      END DO
      CALL SOLH(N,LDE,E,1,DEL,IP)
      DO MMM=1,N-2
         MP=N-MMM
         MP1=MP-1
         DO I=MP+1,N
            DEL(I)=DEL(I)+FJAC(I,MP1)*DEL(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 240
         ZSAFE=DEL(MP)
         DEL(MP)=DEL(I)
         DEL(I)=ZSAFE
 240     CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVSEU
C
      SUBROUTINE DEC (N, NDIM, A, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DEC -------------------------
      END
C
C
      SUBROUTINE SOL (N, NDIM, A, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
C----------------------- END OF SUBROUTINE SOL -------------------------
      END
c
c
      SUBROUTINE DECH (N, NDIM, A, LB, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J,LB,NA
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG
C  MATRIX WITH LOWER BANDWIDTH LB
C  INPUT..
C     N = ORDER OF MATRIX A.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1).
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A SLIGHT MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,NA
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,NA
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECH ------------------------
      END
C
C
      SUBROUTINE SOLH (N, NDIM, A, LB, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1,LB,NA
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX A.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH.
C    LB = LOWER BANDWIDTH OF A.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DECH HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
C----------------------- END OF SUBROUTINE SOLH ------------------------
      END
C
      SUBROUTINE DECC (N, NDIM, AR, AI, IP, IER)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
C  ------ MODIFICATION FOR COMPLEX MATRICES --------
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
C     (AR, AI) = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
C     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
C     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    REAL PART.
C     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    IMAGINARY PART.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(AR(I,K))+DABS(AI(I,K)) .GT.
     &          DABS(AR(M,K))+DABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(K,K)
        AI(M,K) = AI(K,K)
        AR(K,K) = TR
        AI(K,K) = TI
 20     CONTINUE
        IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = KP1,N
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        DO 50 J = KP1,N
          TR = AR(M,J)
          TI = AI(M,J)
          AR(M,J) = AR(K,J)
          AI(M,J) = AI(K,J)
          AR(K,J) = TR
          AI(K,J) = TI
          IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = KP1,N
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = KP1,N
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = KP1,N
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (DABS(AR(N,N))+DABS(AI(N,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECC ------------------------
      END
C
C
      SUBROUTINE SOLC (N, NDIM, AR, AI, BR, BI, IP)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
C    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    (BR,BI) = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    (BR,BI) = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        DO 10 I = KP1,N
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 10       CONTINUE
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K)
        PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K)
        PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        DO 30 I = 1,KM1
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 30       CONTINUE
 40     CONTINUE
 50     CONTINUE
        DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1)
        PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1)
        PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
      RETURN
C----------------------- END OF SUBROUTINE SOLC ------------------------
      END
C
C
      SUBROUTINE DECHC (N, NDIM, AR, AI, LB, IP, IER)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
C  ------ MODIFICATION FOR COMPLEX MATRICES --------
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
C     (AR, AI) = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
C     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
C     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    REAL PART.
C     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    IMAGINARY PART.
C     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (LB .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
          IF (DABS(AR(I,K))+DABS(AI(I,K)) .GT.
     &          DABS(AR(M,K))+DABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(K,K)
        AI(M,K) = AI(K,K)
        AR(K,K) = TR
        AI(K,K) = TI
 20     CONTINUE
        IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = KP1,NA
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        DO 50 J = KP1,N
          TR = AR(M,J)
          TI = AI(M,J)
          AR(M,J) = AR(K,J)
          AI(M,J) = AI(K,J)
          AR(K,J) = TR
          AI(K,J) = TI
          IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = KP1,NA
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = KP1,NA
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = KP1,NA
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (DABS(AR(N,N))+DABS(AI(N,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECHC -----------------------
      END
C
C
      SUBROUTINE SOLHC (N, NDIM, AR, AI, LB, BR, BI, IP)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
C    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    (BR,BI) = RIGHT HAND SIDE VECTOR.
C    LB = LOWER BANDWIDTH OF A.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    (BR,BI) = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      IF (LB .EQ. 0) GO TO 25
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        DO 10 I = KP1,MIN0(N,LB+K)
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 10       CONTINUE
 20     CONTINUE
 25     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K)
        PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K)
        PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        DO 30 I = 1,KM1
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 30       CONTINUE
 40     CONTINUE
 50     CONTINUE
        DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1)
        PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1)
        PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
      RETURN
C----------------------- END OF SUBROUTINE SOLHC -----------------------
      END
C
      SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
      REAL*8 A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  A.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J)
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECB ------------------------
      END
C
C
      SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
      REAL*8 A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    B      RIGHT HAND SIDE VECTOR.
C    IP     PIVOT VECTOR OBTAINED FROM DECB.
C  DO NOT USE IF DECB HAS SET IER .NE. 0.
C  OUTPUT..
C    B      SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
C----------------------- END OF SUBROUTINE SOLB ------------------------
      END
C
      SUBROUTINE DECBC (N, NDIM, AR, AI, ML, MU, IP, IER)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL
C                PART) AND AI (IMAGINARY PART)  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
      AR(I,J) = 0.D0
      AI(I,J) = 0.D0
  5   CONTINUE
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(AR(I,K))+DABS(AI(I,K)) .GT.
     &          DABS(AR(M,K))+DABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(MD,K)
        AI(M,K) = AI(MD,K)
        AR(MD,K) = TR
        AI(MD,K) = TI
 20     IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = MD1,MDL
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          TR = AR(M,J)
          TI = AI(M,J)
          IF (M .EQ. MM) GO TO 35
          AR(M,J) = AR(MM,J)
          AI(M,J) = AI(MM,J)
          AR(MM,J) = TR
          AI(MM,J) = TI
 35       CONTINUE
          IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48
          JK = J - K
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = MD1,MDL
            IJK = I - JK
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = MD1,MDL
            IJK = I - JK
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = MD1,MDL
            IJK = I - JK
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (DABS(AR(MD,N))+DABS(AI(MD,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECBC ------------------------
      END
C
C
      SUBROUTINE SOLBC (N, NDIM, AR, AI, ML, MU, BR, BI, IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B ,
C                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION.
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART).
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART).
C    IP     PIVOT VECTOR OBTAINED FROM DECBC.
C  DO NOT USE IF DECB HAS SET IER .NE. 0.
C  OUTPUT..
C    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART).
C-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(IMD) = BR(IMD) + PRODR
          BI(IMD) = BI(IMD) + PRODI
 10     CONTINUE
 20     CONTINUE
 25     CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        DEN=AR(MD,K)*AR(MD,K)+AI(MD,K)*AI(MD,K)
        PRODR=BR(K)*AR(MD,K)+BI(K)*AI(MD,K)
        PRODI=BI(K)*AR(MD,K)-BR(K)*AI(MD,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(IMD) = BR(IMD) + PRODR
          BI(IMD) = BI(IMD) + PRODI
 30       CONTINUE
 40     CONTINUE
        DEN=AR(MD,1)*AR(MD,1)+AI(MD,1)*AI(MD,1)
        PRODR=BR(1)*AR(MD,1)+BI(1)*AI(MD,1)
        PRODI=BI(1)*AR(MD,1)-BR(1)*AI(MD,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
 50   CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE SOLBC ------------------------
      END
c
C
      subroutine elmhes(nm,n,low,igh,a,int)
C
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real*8 a(nm,n)
      real*8 x,y
      real*8 dabs
      integer int(igh)
C
C     this subroutine is a translation of the algol procedure elmhes,
C     num. math. 12, 349-368(1968) by martin and wilkinson.
C     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
C
C     given a real general matrix, this subroutine
C     reduces a submatrix situated in rows and columns
C     low through igh to upper hessenberg form by
C     stabilized elementary similarity transformations.
C
C     on input:
C
C      nm must be set to the row dimension of two-dimensional
C        array parameters as declared in the calling program
C        dimension statement;
C
C      n is the order of the matrix;
C
C      low and igh are integers determined by the balancing
C        subroutine  balanc.      if  balanc  has not been used,
C        set low=1, igh=n;
C
C      a contains the input matrix.
C
C     on output:
C
C      a contains the hessenberg matrix.  the multipliers
C        which were used in the reduction are stored in the
C        remaining triangle under the hessenberg matrix;
C
C      int contains information on the rows and columns
C        interchanged in the reduction.
C        only elements low through igh are used.
C
C     questions and comments should be directed to b. s. garbow,
C     applied mathematics division, argonne national laboratory
C
C     ------------------------------------------------------------------
C
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
C
      do 180 m = kp1, la
       mm1 = m - 1
       x = 0.0d0
       i = m
C
       do 100 j = m, igh
          if (dabs(a(j,mm1)) .le. dabs(x)) go to 100
          x = a(j,mm1)
          i = j
  100   continue
C
       int(m) = i
       if (i .eq. m) go to 130
C    :::::::::: interchange rows and columns of a ::::::::::
       do 110 j = mm1, n
          y = a(i,j)
          a(i,j) = a(m,j)
          a(m,j) = y
  110   continue
C
       do 120 j = 1, igh
          y = a(j,i)
          a(j,i) = a(j,m)
          a(j,m) = y
  120   continue
C    :::::::::: end interchange ::::::::::
  130   if (x .eq. 0.0d0) go to 180
       mp1 = m + 1
C
       do 160 i = mp1, igh
          y = a(i,mm1)
          if (y .eq. 0.0d0) go to 160
          y = y / x
          a(i,mm1) = y
C
          do 140 j = m, n
  140      a(i,j) = a(i,j) - y * a(m,j)
C
          do 150 j = 1, igh
  150      a(j,m) = a(j,m) + y * a(j,i)
C
  160   continue
C
  180 continue
C
  200 return
C    :::::::::: last card of elmhes ::::::::::
      end
