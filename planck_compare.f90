        FUNCTION PLANCK_oper(TT)  
! * 3.01.2002 - Integrated Planck function for "separate" intervals *
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER NAME*17,CHANNEL_NUMBER0*3  
PARAMETER (NT=8001)
        DIMENSION TTT(NT),PPP(NT)
        COMMON/CHNN/CHANNEL_NUMBER0
        DATA IB/0/
        IF (IB.EQ.0)THEN 
	  IB=1
!*** 
        OPEN(10,FILE='PLANCKseparate.'//CHANNEL_NUMBER0)
        READ(10,*)N_T
	IF(NT.NE.N_T)THEN
	WRITE(*,*)'N_T NT',N_T,NT
	STOP
	END IF
         DO I=1,N_T
         READ(10,*)TTT(I),PPP(I)
         END DO
          H=(TTT(N_T)-TTT(1))/(N_T-1)
          TS=TTT(1)
        END IF 
         I=(TT-TS)/H+1
          C=(TT-TTT(I))/H
         PLANCK_oper= (1.D0-C)*PPP(I)+C*PPP(I+1)
         END            
