C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR DMV10 FOR THE MOTION OF A RIGID BODY
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
        include 'dmv10.f'
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Q(4),AM(3),RPAR(20),IPAR(20)
        REAL TIME0,TIME1
        EXTERNAL POTENP
C --- INITIAL VALUES
C ANGULAR MOMENTUM
        AM(1)=1.8D0
        AM(2)=0.4D0
        AM(3)=-0.9D0
C QUATERNION
        Q(1)=1.0D0
        Q(2)=0.0D0
        Q(3)=0.0D0
        Q(4)=0.0D0
C --- MOMENTS OF INERTIA
        AI1=0.6D0
        AI2=0.8D0
        AI3=1.0D0
        RPAR(11)=AI1
        RPAR(12)=AI2
        RPAR(13)=AI3
C ---
        H=0.01
        XEND=10.0D0
        NSTEP=XEND/H
        H=XEND/NSTEP
        WRITE (6,*) 'XEND=',XEND,' H=',H,' NSTEP=',NSTEP
C ---
        DO I=1,10
          RPAR(I)=0.0D0
          IPAR(I)=0
        END DO
C
        WRITE (6,*) '--- INITIAL CONDITION'
        WRITE (6,*) '  AM    ',AM(1),AM(2),AM(3)
        WRITE (6,*) '  QQ    ',Q(1),Q(2),Q(3),Q(4)
        CALL HAMIL(Q,AM,HAM0,RPAR,IPAR)
        CALL CPU_TIME(TIME0)
C ---
        CALL DMV10 (AM,Q,POTENP,H,NSTEP,RPAR,IPAR)
C ---
        CALL CPU_TIME(TIME1)
        CALL HAMIL(Q,AM,HAM1,RPAR,IPAR)
        WRITE (6,*) '--- SOLUTION AT ENDPOINT'
        WRITE (6,*) '  AM    ',AM(1),AM(2),AM(3)
        WRITE (6,*) '  QQ    ',Q(1),Q(2),Q(3),Q(4)
        WRITE (6,*) 'ERR HAM=',(HAM1-HAM0)/HAM0,' TIME=',TIME1-TIME0
        STOP
        END
C
        SUBROUTINE HAMIL (Q,AM,HAM,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z) 
        DIMENSION Q(4),AM(3)
        DIMENSION IPAR(*),RPAR(*)
        CALL POTEN(Q,POT,RPAR,IPAR)
        HAM=AM(1)**2/RPAR(11)+AM(2)**2/RPAR(12)+AM(3)**2/RPAR(13)
        HAM=HAM/2.0D0+POT
        RETURN
        END 
c
        SUBROUTINE POTEN(Q,POT,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z) 
        DIMENSION Q(4)
        DIMENSION IPAR(*),RPAR(*)
        POT=Q(1)**2-Q(2)**2-Q(3)**2+Q(4)**2
        POT=0.0D0
        RETURN
        END 
c
        SUBROUTINE POTENP(Q,POTP,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z) 
        DIMENSION Q(4),POTP(3)
        DIMENSION IPAR(*),RPAR(*)
        POTP(1)=-2*(Q(1)*Q(2)+Q(3)*Q(4))
        POTP(2)=-2*(Q(1)*Q(3)-Q(2)*Q(4))
        POTP(1)=0.0d0
        POTP(2)=0.0d0
        POTP(3)=0.0d0
        RETURN
        END 
