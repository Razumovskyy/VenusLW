        FUNCTION PLANCK_oper(TT) 
!  Integrated Planck function for "separate" intervals *
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER NAME*17,CHANNEL_NUMBER*2,CHO*2   !  
        PARAMETER (NT=8001)
        DIMENSION T(NT),P(NT)
        COMMON/CHNN/CHANNEL_NUMBER
        DATA CHO/'00'/,IBG/0/
        IF (CHO/=CHANNEL_NUMBER)THEN 
        CHO=CHANNEL_NUMBER
		IF(IBG==1)CLOSE(10)
		IBG=1
        NAME='PLANCKseparate.'//CHANNEL_NUMBER  
        OPEN(10,FILE=NAME)
        READ(10,*)N_T
	IF(NT.NE.N_T)THEN
	WRITE(*,*)'N_T NT',N_T,NT
	STOP
	END IF
         DO I=1,N_T
         READ(10,*)T(I),P(I)
         END DO
          H=(T(N_T)-T(1))/(N_T-1)
          TS=T(1)
        END IF 
         I=(TT-TS)/H+1
          C=(TT-T(I))/H
         PLANCK_oper= (1.D0-C)*P(I)+C*P(I+1)
         END            
