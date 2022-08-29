!******************************************************************
	SUBROUTINE FLUX_PT  !   5 Aug.2022 ( date of 12 April,2002 ) 
!!******************************************************************
 USE INITIAL_PT
	IMPLICIT REAL*8 (a-h,o-z)
	PARAMETER(NSTEP=10,KMAX=NSTEP,DELT_CRT=0.0001, 	& 
	STEP=1.D0,H=1.D0/2048.D0,NQ=2049,PI=3.141593,TLIM=0.001,	&
! #	NANGLE=1,ND=110) ! ND - dimension of internal arrays
	NANGLE=3,ND=110) ! ND - dimension of internal arrays
!* ND	- maximal number (for DIMENSION) of the J-grid.
!* NQ	- ... of the wavenumber grid (per inverse cm).
!* NSTEP - may be considered as 'constant' = 10
!*.... The name of FILE with absorption coefficients at J-levels:
!* ---------------------------------------------------------------*
	DIMENSION  COSINE(NANGLE),WEIGHT(NANGLE),&
	 CL0(ND),CL1(ND),PL(ND),TZI(ND),EMIS_ON(NANGLE),TAUMA(ND) &
     ,Emis_Up(ND,NANGLE),Emis_Down(ND,NANGLE),E_TAU(ND,NANGLE),FLU(ND),FLD(ND)
	COMMON/R_A/VSTART
	save TZI, CL0, CL1, WEIGHT, COSINE, TAUMA
	save JZM, JZM_1
		DATA IBEG/0/	&
!* ---------- 1-Rays for angle integration ------------- *
!*        ,COSINE/0.6024096/,WEIGHT/0.83/  ! 1/1.66    1.66/2  
!* ---------- 3-Rays for angle integration ------------- *
	,COSINE/0.887298334620741,0.5D0,0.112701665379259/	&
	,WEIGHT/ 0.277777777777778,0.444444444444444,0.277777777777778/
!* ---------- 5-Rays for angle integration ------------- *
!#	,COSINE/0.4691008E-1,0.2307653,0.5,0.7692347,0.9530899/	&
!#	,WEIGHT/0.1184634,0.2393143,0.2844444,0.2393143,0.1184634/
!* --------------------------------------------------------------- *
	PLANK2(V1,TTT1)=3.7403E-8*V1**3/(DEXP(1.43868*V1/TTT1)-1.)*2.
!*	* /(2.*PI)						! for INTENSITIES
	VSREL=VSTART
		IF(IBEG == 0)THEN
                        JZM=JMAX
         IF(JZM>ND)THEN
          WRITE(*,*)'******* Increase ND-parameter !!! *******'
          STOP
         END IF

			JZM_1=JZM-1
					DO J=1,JZM
					TZI(J)=TEMPER(Z(J))
					END DO

		IBEG=1

	NTT	= STEP/H + 1.01
	NT=NQ
	IF(NT /= NTT) PAUSE 0
!*
				DELTAV=STEP
				NTMAX=NT
				NN=NTMAX-1
!* --------------------------------------------- *
!******************* CL0,CL1 for TAUMA - matrix ****************
	CL0(1) = 0.D0
	CL1(1) = 0.D0
!*
	DO 1 J = 2, JZM
!*
	Z1=Z(J-1)
	Z2=Z(J)
	PSIL=Z1
	PSIL1=Z2
	GA=Z2-Z1
	GB=0.5*(Z2+Z1)*GA
	CC=PSIL-PSIL1
	CL1(J)=(PSIL*GA-GB)/CC
	CL0(J)=(GB-PSIL1*GA)/CC
!*
1	CONTINUE

! COS(1) > COS(2) ... etc.
         DO WHILE (IPR/=NANGLE)
		 IPR=1
		  DO IA=2,NANGLE
          CO1=COSINE(IA-1)
		  CO2=COSINE(IA)
		    IF(CO1<CO2)THEN
            W1=WEIGHT(IA-1)
			WEIGHT(IA-1)=WEIGHT(IA)
			WEIGHT(IA)=W1
			COSINE(IA-1)=CO2
			COSINE(IA)=CO1
			ELSE
			IPR=IPR+1
			END IF
		  END DO
		 END DO

nnnn = 0
	END IF
	nnnn = nnnn+1
!C ***********************************
!C
!C *** STEP BY STEP
!  upon the wave number intervals where Planck=const ***
 J_STEP=0
 DO  KKK = 1, KMAX
	write(*,*) 'nnnn=', nnnn
	if (nnnn>1) write(*,*) 'kkk=', kkk
  J_STEP=J_STEP+1
  H_STP=H
	
 !*
    FREC=VSREL+STEP*(KKK-0.5)
!	if (nnnn>1) write(*,*) 'passed'
    IFT=(KKK-1)*(NQ-1)
!	if (nnnn>1) write(*,*) 'passed'
		DO JJ=1,JZM
		T=TZI(JJ)
		if (nnnn>1) write(*,*) jj, T, frec
		PL(JJ)=PLANK2(FREC,T)
		if (nnnn>1) write(*,*) jj, pl(jj)
		!if (nnnn>1)	pause 100
		END DO
		if (nnnn>1) then
		WRITE(*,*) TZI
		!write(*,*) "PL=", PL; pause 10
		end if
                FLU=0.D0
                FLD=0.D0

!* Loop over wavenumber points *
    DO II=1,NQ
     IF(II.EQ.1.OR.II.EQ.NQ) THEN ! Simpson's weight
         WN_W=1.D0*H/3.D0
     ELSE
               IF(II/2*2==II)THEN
                WN_W=4.D0*H/3.D0
               ELSE
                WN_W=2.D0*H/3.D0
                END IF
     END IF  
    TAUMA(JZM)=0.D0   ! Optical thicknesses
     DO IAN=1,NANGLE
     E_TAU(JZM,IAN)=1.D0
     END DO        
    DO J=JZM_1,1,-1 ! 
    TAUMA(J)=TAUMA(J+1) +RABMA(II+IFT,J)*CL0(J+1) &
                       +RABMA(II+IFT,J+1)*CL1(J+1)
     TAT=-(TAUMA(J)-TAUMA(J+1))
     DO IAN=1,NANGLE
     E_TAU(J,IAN)=DEXP(TAT/COSINE(IAN))
     END DO
    END DO

! ------------- Emissivities --------------------
!	if (nnnn>1) then
!		pause 2;
!	end if
!  Downward 
      DO J=1,JZM_1
		BB2=PL(J)
		BB1=PL(J+1)
		TATA2=TAUMA(J)
		TATA1=TAUMA(J+1)
		DELT_TAU=TATA2-TATA1
			DO IAN=1,NANGLE
			X_X=COSINE(IAN)
		IF(DELT_TAU < DELT_CRT)THEN
		EMIS=0.5D0*(BB1+BB2)*DELT_TAU      ! X_X/X_X
		ELSE
		ALPHA=(BB1*TATA2-BB2*TATA1)/DELT_TAU
		BETTA=(BB2-BB1)/DELT_TAU
		EMIS=(ALPHA-BETTA*X_X+BETTA*TATA2 +	&
	(BETTA*X_X-BETTA*TATA1-ALPHA)*DEXP((TATA1-TATA2)/X_X))*X_X
		END IF
        Emis_Down(J,IAN)=EMIS*WEIGHT(IAN) ! Angle weights		
	        	END DO
       END DO

! Upward
			DO IAN=1,NANGLE    ! Surface
            Emis_Up(1,IAN)=PL(1)*COSINE(IAN)*WEIGHT(IAN) ! Angle weights
                        END DO
	DO J=2,JZM
		BB2=PL(J-1)
		BB1=PL(J)
		TATA2=TAUMA(J-1)
		TATA1=TAUMA(J)
		DELT_TAU=TATA2-TATA1
			DO IAN=1,NANGLE
			X_X=COSINE(IAN)
		IF(DELT_TAU < DELT_CRT)THEN
		EMIS=0.5*(BB1+BB2)*DELT_TAU
		ELSE
		ALPHA=(BB1*TATA2-BB2*TATA1)/DELT_TAU
		BETTA=(BB2-BB1)/DELT_TAU
		EMIS=(ALPHA+BETTA*X_X+BETTA*TATA1 -	&
 	(BETTA*X_X+BETTA*TATA2+ALPHA)*DEXP((TATA1-TATA2)/X_X))*X_X
		END IF
                Emis_Up(J,IAN)=EMIS*WEIGHT(IAN)  ! Angle weights
            END DO
        END DO
! --------------------------------------------------------

!		if (nnnn>1) then
!			pause 3;
!		end if
!****** Downward fluxes
        EMIS_ON=0.D0
	DO  LL = JZM_1,1,-1  !   Loop over ZI-points ***
             F_=0.D0
                      DO IAN=1,NANGLE
                      EMIS_ON(IAN)=Emis_DOWN(LL,IAN)+E_TAU(LL,IAN)*EMIS_ON(IAN)
             F_=F_+EMIS_ON(IAN) 
                      END DO
              FLUXDO(LL)=FLUXDO(LL)+F_*WN_W
     END DO! End of the loop over Zi points *

!		if (nnnn>1) then
!			pause 3;
!		end if
	 !****** Upward fluxes ******
             F_=0.D0
              DO IAN=1,NANGLE
              EMIS_ON(IAN)=Emis_UP(1,IAN)  ! Here ALBEDO can be taken into account
              F_=F_+EMIS_ON(IAN)
              END DO 
 FLUXUP(1)=FLUXUP(1)+F_*WN_W
	DO  LL= 2,JZM  !   Loop over ZI-points ***
             F_=0.D0
              DO IAN=1,NANGLE
             EMIS_ON(IAN)=Emis_UP(LL,IAN)+E_TAU(LL-1,IAN)*EMIS_ON(IAN)
              F_=F_+EMIS_ON(IAN) 
              END DO 
FLUXUP(LL)=FLUXUP(LL)+F_*WN_W
     END DO! End of the loop over Zi points *
 END DO ! End of the loop over wavenumber points *
 END DO ! Finish
 !write(*,*) FLUXUP(1), FLUXUP(10), FLUXUP(101), FLUXDO(1), FLUXDO(10), FLUXDO(101); pause 1
 
END       
!***************************************************
	FUNCTION TEMPER(ZZZ)
 USE INITIAL_PT
	REAL*8 temper
!*
		DO J=2,JMAX
		IF(ZZZ >= Z(J-1).AND.ZZZ <= Z(J)) GOTO 1
		END DO
		WRITE(*,*)' TEMPER '
		STOP
 1		C1=(Z(J)-ZZZ)/(Z(J)-Z(J-1))
		C2=1.-C1
		TEMPER=C1*T1(J-1)+C2*T1(J)
		END
! ***************************************************** !
