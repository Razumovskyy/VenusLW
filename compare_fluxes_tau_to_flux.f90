       SUBROUTINE FROMTAUTOFLUX(WES,FLUXDO,FLUXUP,QT,D_TAU,fluxuplbl)
!*** LONG-WAVE calculations by "K-Distribution technique" ***
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER NAME*10, AAA*80
      PARAMETER (JZM=101,NANGLE=3,TLIM=25.,DELT_CRT=0.001) ! 5
      DIMENSION FLUXUP(JZM),FLUXDO(JZM),QT(JZM),XEF(JZM),fluxuplbl(jzm) &
     ,WEIGHT(NANGLE),COSINE(NANGLE),TAUMA(JZM),D_TAU(JZM)
      COMMON /ZPT/Z(JZM),P(JZM),T(JZM),RO(JZM)
       COMMON/PLAN/PL(JZM)
!*  ---------- 5-Rays for angle integration -------------  *
!###      DATA COSINE/0.4691008E-1,0.2307653,0.5,0.7692347,0.9530899/  &
!###     ,WEIGHT/0.1184634,0.2393143,0.2844444,0.2393143,0.1184634/
 !* ---------- 3-Rays for angle integration ------------- *
	DATA COSINE/0.887298334620741,0.5D0,0.112701665379259/	&
	,WEIGHT/ 0.277777777777778,0.444444444444444,0.277777777777778/
        DO J=1,JZM
      FLUXUP(J)=0.
      FLUXDO(J)=0.
        END DO
!* ------------------------------------- *
         JZM_1=JZM-1
         TAUMA(JZM)=0.
         DO J=JZM_1,1,-1
         TAUMA(J)=TAUMA(J+1)+D_TAU(J) !Attention
         END DO
!* * * * * * * * * * * * * * * *  
!*  *** The loop by ZI-points ***
!*
      DO 410 I=1,JZM 
!****** Downward fluxes
          IF(I.NE.JZM) THEN
           DO  J=I,JZM_1
         BB2=PL(J)
         BB1=PL(J+1)
         TATA1=TAUMA(J+1)           
         TATA2=TAUMA(J) 
          TAU_DO=TAUMA(I)
          DELT_TAU=TATA2-TATA1
           TTT1=0.
           IF(TAU_DO-TATA2.LE.TLIM) THEN
             DO IAN=1,NANGLE
             X_X=COSINE(IAN)
           IF(DELT_TAU.LT.DELT_CRT)THEN
          EMIS=X_X*0.5*(BB1+BB2)*DELT_TAU/X_X*DEXP((TATA2-TAU_DO)/X_X)
		ELSE
		ALPHA=(BB1*TATA2-BB2*TATA1)/DELT_TAU
           BETTA=(BB2-BB1)/DELT_TAU
           EMIS= X_X*TRACE_DN(ALPHA,BETTA,TATA1,TATA2,TAU_DO,X_X)
		 END IF 
         TTT1=TTT1+EMIS*WEIGHT(IAN)
              END DO
           END IF
!*                   
      FLUXDO(I)=FLUXDO(I)+TTT1*WES
	   END DO 
         	END IF

!****** Upward fluxes: The loop by ZV-mesh points ******
        IF(I.GT.1)THEN       
        DO  J=2,I
         BB2=PL(J-1)
         BB1=PL(J)
         TATA2=TAUMA(J-1)           
         TATA1=TAUMA(J) 
          TAU_UP=TAUMA(I)
          DELT_TAU=TATA2-TATA1
                TTT1=0.
           IF(TATA1-TAU_UP.LE.TLIM) THEN
            DO IAN=1,NANGLE
            X_X=COSINE(IAN)
           IF(DELT_TAU.LT.DELT_CRT)THEN
           EMIS=X_X*0.5*(BB1+BB2)*DELT_TAU/X_X*EXP((TAU_UP-TATA1)/X_X)
		 ELSE
           ALPHA=(BB1*TATA2-BB2*TATA1)/DELT_TAU
           BETTA=(BB2-BB1)/DELT_TAU
           EMIS= X_X*TRACE_UN(ALPHA,BETTA,TATA1,TATA2,TAU_UP,X_X)
               END IF                   
         TTT1=TTT1+EMIS*WEIGHT(IAN)
               END DO
           END IF

      FLUXUP(I)=FLUXUP(I)+TTT1*WES
	   END DO
       END IF  
!* *** Surface's emmitance ***
          TTT1=0.
	 TAUSUR=TAUMA(1)-TAUMA(I)
                IF(TAUSUR.LE.TLIM)THEN
                   DO IAN=1,NANGLE
            X_X=COSINE(IAN)
          TTT1=TTT1+X_X*EXP(-TAUSUR/X_X)*WEIGHT(IAN)
                  END DO
       END IF
        FLUXUP(I)=FLUXUP(I)+PL(1)*TTT1*WES
!* ------------------------------------------------------------------ *

  410   CONTINUE

!* * * * * * * * * * * * * * * *  
  101    CONTINUE
!*  -----------  QOOLING RATES, PRINT, etc. -------------- *
            DO J=1,JZM
!ccc            XEF(J)=fluxuplbl(j)- FLUXDO(J) ! FLUXUP(J)-FLUXDO(J)
       XEF(J)= FLUXUP(J)-FLUXDO(J)
            END DO
!**********************************************************************
!*   'LAYER'      Derivatives                 *
         DO J=2,jZM
      QT(J)=-(XEF(J)-XEF(J-1))/(P(J)-P(J-1))*8.442/1013.25
         END DO
        QT(1)=QT(2)
        END            

!* Subroutine to calculate the layer (TAU1,TAU2) emissioN (up/DOWN) at 
!* (TAU_UP/DOWN) point for trace X=cos(zenith angle)
      FUNCTION TRACE_DN(ALPHA,BETTA,TAU1,TAU2,TAU_UP,X)
      IMPLICIT REAL*8 (A-H,O-Z)
         EMIS=(ALPHA-BETTA*X+BETTA*TAU2 +(BETTA*X-BETTA*TAU1-ALPHA)*DEXP((TAU1-TAU2)/X))
         TRACE_DN=EMIS*DEXP((TAU2-TAU_UP)/X)
          END
       FUNCTION TRACE_UN(ALPHA,BETTA,TAU1,TAU2,TAU_DO,X)
              IMPLICIT REAL*8 (A-Z) 
         EMIS=(ALPHA+BETTA*X+BETTA*TAU1 -(BETTA*X+BETTA*TAU2+ALPHA)*DEXP((TAU1-TAU2)/X))
         TRACE_UN=EMIS*DEXP((TAU_DO-TAU1)/X)
          END
