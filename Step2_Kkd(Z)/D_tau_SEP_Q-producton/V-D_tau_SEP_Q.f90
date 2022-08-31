!* FILE='D_tau_SEP_Q.FOR' *** 23.06.2022 - April 24 - May 11,2001 - January 5,2002 *
!* To obtain effectiv DELTA(TAU)i from LBL-Separate DOWN Fluxes *
       IMPLICIT REAL*8 (A-H,O-Z)
          integer*2 result
		  CHARACTER*260 FILNAME,MET*2
	CHARACTER CHANNEL_NUMBER0*3,FI*20,BAN*3 
       PARAMETER (JDI=101,N_ATM=1) ! Number of the atmospheric models
       CHARACTER*80 ATM_PROF(N_ATM),LBL(N_ATM),K_FILE(N_ATM)
DIMENSION D_TAU(JDI),C_K(JDI),FDO_LBL(JDI),FUP_LBL(JDI),Q_LBL(JDI),FDO_KDS(JDI),FUP_KDS(JDI),Q_KDS(JDI)
       COMMON /JMX/JZM   ;  COMMON/ZPT/Z(JDI),P(JDI),T(JDI),RO(JDI)
       COMMON/PLAN/PL(JDI)
	COMMON/CHNN/CHANNEL_NUMBER0
!***         Fit Setting        ***
!********** BAND SETTINGS **********!
OPEN(2001,FILE='BAND_V1-V2.PT')
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

ATM_PROF(1)='./ATMOSPHERES/'//FI
  LBL(1)='Flux-Term'//CHANNEL_NUMBER0
  K_FILE(1)='K(Z).'//CHANNEL_NUMBER0

! * -------------------------------------------------------- *
             weight=1.0
! *          READ(*,*) WEIGHT  ! Main fitting parameter "WEIGHT"
            IF(WEIGHT.LT.0.)STOP

             DO N_A=1,N_ATM
      CALL TAU_(WEIGHT,ATM_PROF(N_A),LBL(N_A),D_TAU,FDO_LBL,FUP_LBL,Q_LBL )
  !*  -----------  QOOLING RATES for LBL -------------- *
!**********************************************************************
!*   'LAYER'      Derivatives                 *
           DO J=2,JZM
		   XEJ=FUP_LBL(J)-FDO_LBL(J)
		   XEJM=FUP_LBL(J-1)-FDO_LBL(J-1)
        Q_LBL(J)=-(XEJ-XEJM)/(P(J)-P(J-1))*3.963/1013.25
           END DO
          Q_LBL(1)=Q_LBL(2)

      CALL FROMTAUTOFLUX(WEIGHT,FDO_KDS,FUP_KDS,Q_KDS,D_TAU,fup_lbl)  ! "Fup,Fdo,Q" from "TAU" 
      	GOTO 100
             WRITE(*,*)' With correction !!! '
! ********************** CORRECTION of D_TAU for Q improving ***********************************
            DO N_COR=1,100
             J_COR=1 ; JZM_C=99
            JZM_1=JZM-1
              DO J=JZM_C,J_COR,-1
      	D_TAU(J)=D_TAU(J) -(Q_LBL(J+1)-Q_KDS(J+1))/3.963*(P(J+1)-P(J))*2./(PL(J)+PL(J+1))/1013.25
      	  END DO
              CALL FROMTAUTOFLUX(WEIGHT,FDO_KDS,FUP_KDS,Q_KDS,D_TAU,fup_lbl) ! "Fup,Fdo,Q" from "TAU"
      	END DO
!**********************************************************************************************
       100   CONTINUE

!* K-coefficients *
            JZM_1=JZM-1
            OPEN(21,FILE=K_FILE(N_A))
      	C_K(JZM)=0.D0
      	DO J=JZM_1,1,-1
		DZ=Z(J+1)-Z(J)     ! <**** ATTENTION!
      	C_K(J)=D_TAU(J)*2.D0/DZ-C_K(J+1)
      	END DO
      	DO J=1,JZM_1
            C_K(J)=C_K(J)/RO(J)
      	WRITE(21,23)Z(J),C_K(J),P(J),T(J),RO(J)
      	END DO
        23  FORMAT(F7.1,E14.6,F16.8,F8.2,E12.4)
      	CLOSE(21)

!* ------ Graphics ------------- *
!* ------ Fdown ---------------*
                        FILNAME='_Fdown_'
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
                 	result=RUNQQ('ToDraw.exe',filname)
      !VVV  	if(result==-1_2)stop
!* ------ Q ---------------*
                      FILNAME='_COOLING'
                      OPEN(8,FILE=FILNAME)
  NQ=JZM
              	WRITE(8,*)M,NQ
                           DO J=1,NQ !JZM
              	L=J  ! JZM-J+1
                    FLBL= Q_LBL(L)
                    FKDS= Q_KDS(L)
              	IF(FKDS.GT.3.)FKDS=3.
                    WRITE(8,*)Z(L),FKDS,FLBL
                            END DO
                      CLOSE(8)
                	result=RUNQQ('ToDraw.exe',filname)
              	if(result==-1_2)stop
!* ------ UP ---------------*
                       FILNAME='_Fup_'
                      OPEN(8,FILE=FILNAME)
              	WRITE(8,*)M,N
                           DO J=1,JZM
              	L=J  ! JZM-J+1
                    FLBL=FUP_LBL(L)
                    FKDS=FUP_KDS(L)
                    WRITE(8,*)Z(L),FKDS,FLBL
                            END DO
                      CLOSE(8)
                	result=RUNQQ('ToDraw.exe',filname)
                   	if(result==-1_2)stop
           	if(result==-1_2)stop

!	                        STOP
!* -------------------------------------------------------- *
 !VVV                   END DO
                               END DO
             
          END
!************************************************************    
        SUBROUTINE TAU_(WES,ATM_PROF,LBL,D_TAU,FDO_LBL,FUP_LBL,Q_LBL)
       IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (JDI=101)
       CHARACTER*80 ATM_PROF,LBL,GAS,OLD_NEW
      COMMON /JMX/JZM   ;  COMMON/ZPT/Z(JDI),P(JDI),T(JDI),RO(JDI)
      DIMENSION C(JDI),R0(32)  &
     ,D_TAU(JDI),FDO(JDI),FDO_LBL(JDI),FUP_LBL(JDI),Q_LBL(JDI),TAU(JDI)
       DATA OLD_NEW/'OOOOOOOOOOOOOOOOOO'/
!* ------ Z,R,T,RO reading from atmospheric profile ------ * 
               IF(OLD_NEW.NE.ATM_PROF)THEN 
      OLD_NEW=ATM_PROF
OPEN(55,FILE=ATM_PROF)  ! Atmospheric profile   
           READ(55,'(A)')TITLE
      READ(55,*)NGAS,JZM
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
      DO J=1,JZM
      READ(55,*)ZZZ,X1,X2,X3,X4
      FUP_LBL(J)=X4    
      FDO_LBL(J)=X3    
      FDO(J)=FDO_LBL(J)  
      END DO
      CLOSE(55) 
              END IF
!*******  Iterations ***********************************
                      coefficient=1E-15 ! ???????????????
      TAU(JZM)=0.0
	D_TAU(JZM)=0.0
	JZMM=JZM-1
jumax=jzm-1
	DO J=JZMM,1,-1
	                   if(j.gt.0)then ! $$$$
	DT1=TAU(J+1)
	DT2=1E6
    if(j.lt.jumax)then 
	ddtt=tau(j+1)-tau(j+2)
	dt1=tau(j+1)+ddtt*0.8
	dt2=dt1+ddtt*3.d0
         	end if
	F1=FLUX(WES,TAU,J,DT1)
	F2=FLUX(WES,TAU,J,DT2)
	      FLU=FDO(J)

       if(flu.gt.f1.and.flu.gt.f2.or.flu.lt.f1.and.flu.lt.f2)then
	write(*,*)' Attention! flu f1 f2 j', flu,f1,f2,j
	!### write(*,*)'dt1 dt2',dt1,dt2
       	end if

1     DT = (DT1+DT2)*0.5D0
    	FF=FLUX(WES,TAU,J,DT)
    	if(f1.ge.flu.and.flu.ge.ff.or.flu.le.ff.and.f1.le.flu)then
        f2=ff
	dt2=dt
	else
	f1=ff
	dt1=dt
	end if
	IF(abs(DT2-DT1)/DT.GT.0.000001D0)GOTO 1
	goto 2
	                                 else
	                                 end if
      dt=coefficient*(ro(j)+ro(j+1))*0.5+tau(j+1)
 2      D_TAU(J)=DT-TAU(J+1)
	TAU(J)=DT
	END DO
! ******* end of the iterations, Writing  **********************
      END
! ****************************************************************
      FUNCTION FLUX(WES,TAUMA,JJ,DT)  ! <--- DOWNWARD Fluxes
! ***  DT=TAUMA(JJ)-TAUMA(JJ+1) *** 
       IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(TLIM=25.,DELT_CRT=0.001,NANGLE=3,JDI=101)
      DIMENSION WEIGHT(NANGLE),COSINE(NANGLE),TAUMA(JDI)
       COMMON /JMX/JZM   ;  COMMON/ZPT/Z(JDI),P(JDI),T(JDI),RO(JDI)
      COMMON/PLAN/PL(JDI)
!CC*  ---------- 5-Rays for angle integration -------------  *
!CC      DATA COSINE/0.4691008E-1,0.2307653,0.5,0.7692347,0.9530899/  
!CC     *,WEIGHT/0.1184634,0.2393143,0.2844444,0.2393143,0.1184634/,TIB/0./
 !CC       DATA COSINE/0.6024096/,WEIGHT/0.83/  ! 1/1.66    1.66/2  
!* ---------- 3-Rays for angle integration ------------- *
	DATA COSINE/0.887298334620741,0.5D0,0.112701665379259/	&
	,WEIGHT/ 0.277777777777778,0.444444444444444,0.277777777777778/
   DATA TIB/0./
        IF(TIB.NE.T(1))THEN
        TIB=T(1)
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
	   if(j.eq.i)tata2=dt           
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
 
      FLUX=FLUX+TTT1*WES
	   END DO 
       END

      FUNCTION TRACE_U(ALPHA,BETTA,TAU1,TAU2,TAU_UP,X)
      IMPLICIT REAL*8 (A-H,O-Z)
         EMIS=(ALPHA-BETTA*X+BETTA*TAU2+(BETTA*X-BETTA*TAU1-ALPHA)*DEXP((TAU1-TAU2)/X))
         TRACE_U=EMIS*DEXP((TAU2-TAU_UP)/X)
          END

