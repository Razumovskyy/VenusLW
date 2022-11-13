!* FILE='D_tau_SEP_Q.FOR' *** 08.07.2022 - April 24 - May 11,2001 - January 5,2002 *
!* N-channels (Total) variant to obtain effectiv DELTA(TAU)i from LBL-Separate DOWN Fluxes *
 IMPLICIT REAL*8 (A-H,O-Z)
 INTEGER*2 result
 CHARACTER FILNAME*260,MET*2
CHARACTER CHANNEL_NUMBER0*3,FI*20,BAN*3,CHANNEL_NUMBER*2,ATM_PR*3,CHN_N(20)*2 
 PARAMETER (JZMD=201,N_ATM=1) ! Number of the atmospheric models
       CHARACTER*80 ATM_PROF(N_ATM),LBL(N_ATM),K_FILE(N_ATM)
DIMENSION D_TAU(JZMD),C_K(JZMD),FDO_LBL(JZMD),FUP_LBL(JZMD),Q_LBL(JZMD),FDO_KDS(JZMD),FUP_KDS(JZMD),Q_KDS(JZMD)
REAL*4 TDO_LBL(JZMD),TUP_LBL(JZMD),TQ_LBL(JZMD),TDO_KDS(JZMD),TUP_KDS(JZMD),TQ_KDS(JZMD)
       COMMON/ZPT/Z(JZMD),P(JZMD),T(JZMD),RO(JZMD),JZM
       COMMON/PLAN/PL(JZMD)
	COMMON/CHNN/CHANNEL_NUMBER
      COMMON/PLANET/CP_CONST
DATA CHN_N/'1_','2_','3_','4_','5_','6_','7_','8_','9_','10' &
          ,'11','12','13','14','15','16','17','18','19','20'/  
!###  CP_CONST=8.442 ; WRITE(*,*)' *** ',CP_CONST,' ->  Earth ***'
!###   CP_CONST=3.963 ; WRITE(*,*)' *** ',CP_CONST,' ->  Mars ***'
CP_CONST=3.963 ; WRITE(*,*)' *** ',CP_CONST,' ->  Venus ***'


!***         Fit Setting        ***
!********** BAND SETTINGS **********!
!********** BAND SETTINGS **********!
OPEN(2001,FILE='band_settings.txt')
READ(2001,2002)CHANNEL_NUMBER0
READ(2001,2002)BAN
2002 FORMAT(A3)
READ(2001,*)VSTART,VFINISH
READ(2001,2022)MET
2022 FORMAT(A2)
READ(2001,2003)FI
2003 FORMAT(A20)
CLOSE(2001)
!****************************************!

! * -------------------------------------------------------- *
            DO N_A=1,N_ATM
ATM_PROF(N_A)='./Atmospheres/'//FI  ! ### FI(NA) <- This time for a SINGLE ATMOSPHERE ***
LBL(N_A)='Fl_N-'//FI
OPEN(2001,FILE=LBL(1))
READ(2001,*)NCH_TOT 
CLOSE(2001)
     ATM_PR=ATM_PROF(N_A)(15:17)
TDO_LBL=0. ; TUP_LBL=0. ; TQ_LBL=0. ; TDO_KDS=0. ; TUP_KDS=0. ; TQ_KDS=0.
  DO NTERM=1,NCH_TOT      ! *** Loop over TERMS ***
CHANNEL_NUMBER=CHN_N(NTERM) 
FDO_LBL=0. ; FUP_LBL=0. ; Q_LBL=0. ; D_TAU=0.
CALL TAU_(ATM_PROF(N_A),LBL(N_A),D_TAU,FDO_LBL,FUP_LBL,Q_LBL,NTERM)
FDO_KDS=0. ; FUP_KDS=0. ; Q_KDS=0.
CALL FROMTAUTOFLUX(NTERM,FDO_KDS,FUP_KDS,Q_KDS,D_TAU,Q_LBL)  ! "Fup,Fdo,Q" from "TAU" 
	GOTO 100           ! ### CORRECTIN NOT USED !!!
       WRITE(*,*)' With correction !!! '  ! ### CORRECTIN NOT USED !!!
! ********************** CORRECTION of D_TAU for Q improving ***********************************
      DO N_COR=1,100
       J_COR=1 ; JZM_C=99
      JZM_1=JZM-1
        DO J=JZM_C,J_COR,-1
	D_TAU(J)=D_TAU(J) -(Q_LBL(J+1)-Q_KDS(J+1))/CP_CONST*(P(J+1)-P(J))*2./(PL(J)+PL(J+1))
	  END DO
CALL FROMTAUTOFLUX(NTERM,FDO_KDS,FUP_KDS,Q_KDS,D_TAU,Q_LBL)  ! "Fup,Fdo,Q" from "TAU" 
	END DO
!**********************************************************************************************
 100   CONTINUE

!* K-coefficients *
      JZM_1=JZM-1
      OPEN(21,FILE='K(Z).'//CHANNEL_NUMBER) 
	C_K(JZM)=0.D0
	DO J=JZM_1,1,-1
	DZ=Z(J+1)-Z(J)
	C_K(J)=D_TAU(J)*2.D0/DZ-C_K(J+1)
	END DO
	DO J=1,JZM_1
      C_K(J)=C_K(J)/RO(J)
	WRITE(21,23)Z(J),C_K(J),P(J),T(J),RO(J)
	END DO
  23  FORMAT(F7.1,E14.6,F16.8,F8.2,E12.4)
	CLOSE(21)

!* ------ Graphics ------------- *
GOTO 297  ! ### < Graphics in each channel is OMITTED
!* ------ Fdown ---------------*
         FILNAME=ATM_PR//'_Fdown_'//CHANNEL_NUMBER
        OPEN(8,FILE=FILNAME)
	N=JZM
	M=2
	WRITE(8,*)M,N
             DO J=1,N !JZM
	L=J  ! JZM-J+1
      FLBL= FDO_LBL(L) !DLOG(FDO_LBL(L))
      FKDS= FDO_KDS(L) !DLOG(FDO_KDS(L))
      WRITE(8,*)Z(L),FKDS,FLBL
              END DO
        CLOSE(8)
!  	result=RUNQQ('ToDraw.exe',filname)
	if(result==-1_2)stop
!* ------ Q ---------------*
         FILNAME=ATM_PR//'_COOLING_'//Channel_number
        OPEN(8,FILE=FILNAME)
	WRITE(8,*)M,N
             DO J=1,N !JZM
	L=J  ! JZM-J+1
      FLBL= Q_LBL(L)
      FKDS= Q_KDS(L)
	IF(FKDS.GT.3.)FKDS=3.
      WRITE(8,*)Z(L),FKDS,FLBL
              END DO
        CLOSE(8)
!  	result=RUNQQ('ToDraw.exe',filname)
	if(result==-1_2)stop
!* ------ UP ---------------*
         FILNAME=ATM_PR//'_Fup_'//CHANNEL_NUMBER
        OPEN(8,FILE=FILNAME)
	WRITE(8,*)M,N
             DO J=1,N !JZM
	L=J  ! JZM-J+1
      FLBL=FUP_LBL(L)
      FKDS=FUP_KDS(L)
      WRITE(8,*)Z(L),FKDS,FLBL
              END DO
        CLOSE(8)
!  	result=RUNQQ('ToDraw.exe',filname)
	if(result==-1_2)stop

297 CONTINUE
!* -------------------------------------------------------- *
TDO_LBL=TDO_LBL+FDO_LBL 
TUP_LBL=TUP_LBL+FUP_LBL
TQ_LBL=TQ_LBL+Q_LBL
TDO_KDS=TDO_KDS+FDO_KDS 
TUP_KDS=TUP_KDS+FUP_KDS
TQ_KDS=TQ_KDS+Q_KDS
            END DO   ! *** DO NTERM=1,NCH_TOT

!* ------ TOTAL Graphics ------------- *
!* ------ Fdown ---------------*
         FILNAME=ATM_PR//'_Fdown_'//CHANNEL_NUMBER
        OPEN(8,FILE=FILNAME)
	N=JZM ! 100
	M=2
	WRITE(8,*)M,N
             DO J=1,N !JZM
	L=J  ! JZM-J+1
      FLBL= TDO_LBL(L) 
      FKDS= TDO_KDS(L) 
      WRITE(8,*)Z(L),FKDS,FLBL
              END DO
        CLOSE(8)
!  	result=RUNQQ('ToDraw.exe',filname)
	if(result==-1_2)stop
!* ------ Q ---------------*
         FILNAME=ATM_PR//'_COOLING_'//Channel_number
        OPEN(8,FILE=FILNAME)
N=JZM
	WRITE(8,*)M,N
             DO J=1,N !JZM
	L=J  ! JZM-J+1
      FLBL= TQ_LBL(L)
      FKDS= TQ_KDS(L)
!###	IF(FKDS.GT.3.)FKDS=3.
      WRITE(8,*)Z(L),FKDS,FLBL
              END DO
        CLOSE(8)
!  	result=RUNQQ('ToDraw.exe',filname)
	if(result==-1_2)stop
!* ------ UP ---------------*
         FILNAME=ATM_PR//'_Fup_'//CHANNEL_NUMBER
        OPEN(8,FILE=FILNAME)
	N=JZM ! 100
	WRITE(8,*)M,N
             DO J=1,N !JZM
	L=J  ! JZM-J+1
      FLBL=TUP_LBL(L)
      FKDS=TUP_KDS(L)
      WRITE(8,*)Z(L),FKDS,FLBL
              END DO
        CLOSE(8)
!  	result=RUNQQ('ToDraw.exe',filname)

                        END DO  !*** DO N_A=1,N_ATM
       
          END
!************************************************************    
        SUBROUTINE TAU_(ATM_PROF,LBL,D_TAU,FDO_LBL,FUP_LBL,Q_LBL,NTERM)
        PARAMETER (JZMD=201)
       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*80 ATM_PROF,LBL,GAS,OLD_NEW,TITLE
      COMMON /ZPT/Z(JZMD),P(JZMD),T(JZMD),RO(JZMD),JZM
      COMMON/PLANET/CP_CONST
      DIMENSION C(JZMD),R0(32)  &
     ,D_TAU(JZMD),FDO(JZMD),FDO_LBL(JZMD),FUP_LBL(JZMD),Q_LBL(JZMD),TAU(JZMD)
 ALLOCATABLE FLUX_N(:,:)
  SAVE FLUX_N
       DATA OLD_NEW/'OOOOOOOOOOOOOOOOOO'/
!* ------ Z,R,T,RO reading from atmospheric profile ------ * 
               IF(OLD_NEW.NE.ATM_PROF)THEN 
      OLD_NEW=ATM_PROF
         OPEN(55,FILE=ATM_PROF)  ! Atmospheric profile 
           READ(55,'(A)')TITLE
      READ(55,*)NGAS,JZM
      IF(JZM.GT.JZMD) WRITE(*,*)' Number of levels',JZM,'>',JZMD,' in DIMENSION'
      IF(JZM.GT.JZMD) STOP
         DO I=1,NGAS
         READ(55,'(A)')GAS
         END DO
      DO J=1,JZM
      READ(55,*)Z(J),P(J),T(J),( R0(I), I = 1, NGAS )
      RO(J)=R0(1)                                      ! Attention - H2O
      P(J)=P(J)*1013.25
      END DO
      CLOSE(55)
!* ------------------------------------- *

!* ----- Fdown from LBL-data ----------- *
      OPEN(55,FILE=LBL)  ! LBL benchmark calculation
	  READ(55,*)NCH_TOT
	  NCH2=NCH_TOT*2
ALLOCATE (FLUX_N(JZM,NCH2))
      DO J=1,JZM
      READ(55,*)ZH,FLUX_N(J,:)
      END DO
      CLOSE(55) 

             END IF
!*******  Iterations ***********************************
                      coefficient=1E-15 ! ###

JI=NTERM*2-1 
      FDO_LBL(1:JZM)=FLUX_N(:,JI)    
      FUP_LBL(1:JZM)=FLUX_N(:,JI+1)     
      FDO(1:JZM)=FDO_LBL(1:JZM)
	  !*  -----------  QOOLING RATES-------------- *
!*   'LAYER'      Derivatives                 *
         DO J=2,JZM
		 EFJU=FDO_LBL(J)-FUP_LBL(J)
		 EFJM=FDO_LBL(J-1)-FUP_LBL(J-1)
      Q_LBL(J)=(EFJU-EFJM)/(P(J)-P(J-1))*CP_CONST 
         END DO
        Q_LBL(1)=Q_LBL(2)
  
!
      TAU=0.0  
	D_TAU=0.0  
	JZMM=JZM-1
jumax=jzm-1
	DO J=JZMM,1,-1
	                   IF(J.GT.0)THEN
	DT1=TAU(J+1)
	DT2=1E6
        	if(j.lt.jumax)then   !	  IF(J<90)THEN
	DDTT=TAU(J+1)-TAU(J+2)
	DT1=TAU(J+1)+DDTT*0.8D0   ! ### ATTENTION !
    DT2=DT1+DDTT*3.D0         ! ### ATTENTION !
      END IF 
	F1=FLUX(NTERM,TAU,T,J,DT1,JZM)
	F2=FLUX(NTERM,TAU,T,J,DT2,JZM)
	      FLU=FDO(J)
        IF(FLU>F1.AND.FLU>F2.OR.FLU<F1.AND.FLU.LT.F2)THEN
	WRITE(*,*)' ATTENTION! F1 FLU F2 J',F1,FLU,F2,J
	WRITE(*,*)'DT1 DT2',DT1,DT2
!!!	pause
         END IF

1     DT = (DT1+DT2)*0.5D0
    	FF=FLUX(NTERM,TAU,T,J,DT,JZM)
        IF(F1>=FLU.AND.FLU>=FF.OR.FLU<=FF.AND.F1<=FLU)THEN
        F2=FF
	DT2=DT
	ELSE
	F1=FF
	DT1=DT
	END IF
	IF(ABS(DT2-DT1)/DT.GT.0.000001D0) GOTO 1
	GOTO 2
	                                 ELSE
	                                 END IF
      DT=coefficient*(RO(J)+RO(J+1))*0.5+TAU(J+1)  ! ### coefficient
 2      D_TAU(J)=DT-TAU(J+1)
	TAU(J)=DT
	END DO
! ******* end of the iterations, Writing  **********************
 GOTO 630  ! ###  Graph (for control) not used !
 OPEN(629,FILE='LBL-KD.D')
 WRITE(629,*)' 2 100 ' 
  	DO J=1,JZMM
		Fdow=FLUX(NTERM,TAU,T,J,TAU(J),JZM)
 WRITE(629,628)Z(J),fdo(j),Fdow
 628 FORMAT(F7.1,2E12.4)
	end do
	CLOSE(629)
!  	RESULT=RUNQQ('ToDraw.exe','LBL-KD.D')
630 CONTINUE
      END
! ****************************************************************
      FUNCTION FLUX(NTERM,TAUMA,T,JJ,DT,JZM)  ! <--- DOWNWARD Fluxes
! ***  DT=TAUMA(JJ)-TAUMA(JJ+1) *** 
       IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(TLIM=25.,DELT_CRT=0.001,NANGLE=3,JZMD=201)  ! NANGLE=1)
      DIMENSION WEIGHT(NANGLE),COSINE(NANGLE),TAUMA(JZMD),T(JZMD)
!* ---------- 3-Rays for angle integration ------------- *
	DATA COSINE/0.887298334620741,0.5D0,0.112701665379259/	&
	,WEIGHT/ 0.277777777777778,0.444444444444444,0.277777777777778/
!CC*  ---------- 5-Rays for angle integration -------------  *
!CC      DATA COSINE/0.4691008E-1,0.2307653,0.5,0.7692347,0.9530899/  
!CC     *,WEIGHT/0.1184634,0.2393143,0.2844444,0.2393143,0.1184634/,TIB/0./
!CC          DATA COSINE/0.6024096/,WEIGHT/0.83/  ! 1/1.66    1.66/2  
      COMMON/PLAN/PL(201)
	  DATA NTRO/0/
        IF(NTRO/=NTERM)THEN
        NTRO=NTERM
                DO  L=1,JZM
        PL(L)=PLANCK_oper(T(L)) ! Attention- SEPARATE technique
	     end do
        END IF
!* * * * * * * * * * * * * * * *  
!*  *** The loop by ZI-points ***
!*
      I=JJ 
!****** Downward fluxes
            JZM_1=JZM-1
	FLUX=0.D0
           DO  J=I,JZM_1
         BB2=PL(J)
         BB1=PL(J+1)
         TATA2=TAUMA(J)
	   IF(J==I)TATA2=DT           
         TATA1=TAUMA(J+1) 
          TAU_DO=dt
          DELT_TAU=TATA2-TATA1
           TTT1=0.
           IF(TAU_DO-TATA2.LE.TLIM) THEN
             DO IAN=1,nangle
             X_X=COSINE(IAN)
           IF(DELT_TAU.LT.DELT_CRT)THEN
          EMIS=X_X*0.5*(BB1+BB2)*DELT_TAU/X_X*EXP((tata2-TAU_DO)/X_X)
		ELSE
		ALPHA=(BB1*TATA2-BB2*TATA1)/DELT_TAU
           BETTA=(BB2-BB1)/DELT_TAU
           EMIS= X_X*TRACE_U(ALPHA,BETTA,TATA1,TATA2,TAU_DO,X_X)
		 END IF 
         TTT1=TTT1+EMIS*WEIGHT(IAN)
              END DO
           END IF
!*                   
       FLUX=FLUX+TTT1
	   END DO 
       END

      FUNCTION TRACE_U(ALPHA,BETTA,TAU1,TAU2,TAU_UP,X)
      IMPLICIT REAL*8 (A-H,O-Z)
         EMIS=(ALPHA-BETTA*X+BETTA*TAU2+(BETTA*X-BETTA*TAU1-ALPHA)*DEXP((TAU1-TAU2)/X))
         TRACE_U=EMIS*DEXP((TAU2-TAU_UP)/X)
          END

