!! 3 Nov.,2001. ! Attention: NO SELF part in the CKD 2.4
module INITIAL ! Information from atmospheric profile
    parameter (NTH = 20481, NCOMP = 49)  ! #### DLT8
    character LINE_PATH*12, ATM_PATH*14
	character*80 TITLE
	character*5 MOLECULE(NCOMP)
    integer*4 :: NGAS, JMAX 
    real*4, allocatable :: Z(:), RABMA(:,:), P1(:), RO1(:,:), T1(:)
	save Z, RABMA, P1, RO1, T1
  	contains
  	subroutine ATM_PROF_READING(NATT)
  		character NATT*20
        open(66, file='k_coef.chk')
		open(99, file='k_coef.in')
		read(99,'(A)')ATM_PATH
		WRITE(66,*)' The atmospheric conditions from file = ',ATM_PATH
		READ(99,'(A)')LINE_PATH
		WRITE(66,*)' The LBL data base from file = ',LINE_PATH
		CLOSE(99)
		OPEN(55, file=ATM_PATH//NATT)
		READ(55,'(A)')TITLE
		WRITE(*,*)TITLE
		WRITE(66,*)TITLE
		READ(55,*)NGAS,JMAX
		WRITE(*,*)'The number of gases (NGAS) : ',NGAS
		WRITE(*,*)'The number of levels (JMAX) : ',JMAX
		WRITE(66,*)'The number of gases (NGAS) : ',NGAS
		WRITE(66,*)'The number of levels (JMAX) : ',JMAX
		DO I=1,NGAS
			READ(55,'(A)')MOLECULE(I)
		END DO
		WRITE(66,*)'Account atmosphere gases		: ',	&
				(MOLECULE(I),I=1,NGAS)
		ALLOCATE (Z(JMAX),P1(JMAX),T1(JMAX))
        ALLOCATE (RO1(NCOMP,JMAX),RABMA(NTH,JMAX),STAT=IERR)
        IF(IERR/=0)THEN
          WRITE(*,*)' Allocaion is wrong !!!'
          STOP
        END IF

		DO J=1,JMAX
			READ(55,*)Z(J),P1(J),T1(J),( RO1(I,J), I = 1, NGAS )
			WRITE(66,*)Z(J),P1(J),T1(J),( RO1(I,J), I = 1, NGAS )
			!*	WRITE(*,*)Z(J),P1(J),T1(J),( RO1(I,J), I = 1, NGAS )
		END DO
		CLOSE(55)
        CLOSE(66)
    end subroutine ATM_PROF_READING
end module INITIAL 
	!**********************************************************************
	!* This programs intended for LONGWAVE calculations using HITRAN-2016 *
	!* (10 cm^-1 interval with the spectral resolutin of 1/2048 cm^-1)    *
	!* Consists of:													      *
	!*			LBL2016												      *
	!*			LEFTLBL+CENTLBL+RIGHTLBL (1/2048 cm^-1)			          *
	!*	CATAL2016+ SEARC2016 + BLOCK DATA INIT2016,ST_SUM20016DBASE_2016  *
	!*			VOIGT, VV_LOR, DOPLER								      
	!*
	!**********************************************************************
subroutine LBL2016(GASCOMP, IGAS, VS1)
	character*5 MOLECULE
	character*(*) GASCOMP
	real*8 V,VS1,VS,VF,VS_SEARCH,VF_SEARCH,VI
	parameter ( OBR=10., NLIN = 150000,OBR25 = 25.)
	parameter ( BOUNDL=10.,BOUNDD=0.01,BET=1.438786,NT0=10,	&
	NT1=NT0*2,NT2=NT1*2,NT3=NT2*2,NT4=NT3*2,NT5= NT4*2,	&
	NT6=NT5*2,NT= NT6*4+1)
	dimension WISO(125),ISOIND(382)
	common/ST_SUM16/QofT(74,501) ! HITRAN-2017  (74 for the first 12 gases)
	common/PARLIO/V(NLIN),S(NLIN),ALFA(NLIN),ALFAS(NLIN),E(NLIN),	&
	FACT(NLIN),NISO(NLIN),SHIFT(NLIN)
	common/PHIPAR/T1, P1, RO1
	common/SHAPE/SL,AL,ADD,ALAD,VI,MOTYPE,TT,CTF2,SLL,ALAL ,APF2,APS2,Y_LC
	common/MESH/DELTA,H,H0,H1,H2,H3,H4,H5,H6,H7,H8,H9,RK(NT),RK0(NT0),RK0L(NT0),RK0P(NT0),RK1(NT1),RK1L(NT1),RK1P(NT1) &
	,RK2(NT2),RK2L(NT2),RK2P(NT2),RK3(NT3),RK3L(NT3),RK3P(NT3),RK4(NT4),RK4L(NT4),RK4P(NT4),RK5(NT5),RK5L(NT5),RK5P(NT5) &
	,RK6(NT6),RK6L(NT6),RK6P(NT6)
	external VOIGT,VV_LOR,DOPLER

    save I,PSELF,PFOREI,EPS,TST,EK,NNNNNN,ISOOLD,TOLD,APF,APS,A0,A1,A2,A3,J_0_3,A0s,A1s,A2s,A3s,J_0_3s
 	DATA PI/3.141593/, TS/296./, PS/1.0/	&
	,T, P, RO/3*0./, VS,VF/2*1.0D-04/,TOLD/-5./,ISOOLD/-5/
! --------------------------------------------------------------------------- !
	DATA WISO/				& ! Molecular masses (isotops as in GAMACHE's 2017data for 12 gases)
	18.010565,20.014811,19.014780,19.016740,21.020985,20.020956, 20.0,22.0,21.0 & ! H2O (1) 9
     ,43.989830,44.993185,45.994076,44.994045,46.997431,45.997400, 	& ! CO2 (2)- 13
	 47.998322,46.998291,                                         & 
      46.0,                                                        & ! 727 ???
      49.001675,48.0,47.0,44.0  	                                 &
	,47.984745,49.988991,49.988991,48.988960,48.988960,		& ! O3 (3)- 18
      52.,52.,51.,51.,51.,50.,50.,54.,53.,53.,52.,52.,51.      &
	,44.001062,44.998096,44.998096,46.005308,45.005278		& ! N2O (4)- 5
	,27.994915,28.998270,29.999161,28.999130,31.002516,30.002485, 	& ! CO (5) - 9
      30.,32.,33. &
	,16.031300,17.034655,17.037475,18.040830                        & ! CH4 (6) -4
	,31.989830,33.994076,32.994045,36.,35.,34.	& ! O2  (7) - 6
	,29.997989,30.995023,32.002234	& ! NO (8)-3
	,63.961901,65.957695            & ! SO2 (9)-2
	,45.992904	                & ! NO2(10)-1
	,17.026549,18.023583	        & ! NH3(11) -2
! ------  below the HITRAN-12 ----------------- !
	,62.995644		        & ! HNO3(12)
	,17.002740,19.006986,18.008915	& ! OH(13)
	,20.006229	                & ! HF(14)
	,35.976678,37.973729            & ! HCl(15)
	,79.926160,81.924115            & ! HBr(16)
	,127.912297   	                & ! HI(17)
	,50.963768,52.960819            & ! ClO(18)
	,59.966986,61.962780,60.970341,60.966371,61.971231	        & ! OCS(19)
	,30.010565,31.013920,32.014811	& ! H2CO(20)
	,51.971593,53.968644		& ! HOCl(21)
	,28.006147			& ! N2 (22)
	,27.010899,28.014254,28.007933	& ! HCN (23)
	,49.992328,51.989379            & ! CH3Cl (24)
	,34.005480			& ! H2O2 (25)
	,26.015650,27.019005		& ! C2H2 (26)
	,30.046950                      & ! C2H6 (27)
	,33.997238                      & ! PH3 (28)
	,65.991722                      & ! COF2 (29)
	,145.962492                     & ! SF6 (30)
	,33.987721,35.983515,34.987105  & ! H2S (31)
	,46.005480              	& ! HCOOH (32)
	,32.997655                      & ! HO2 (33)
	,15.994915                      & ! O   (34)
	,96.956672,98.953723            & ! ClONO2 (35)
	,29.997989                      & ! NO+ (36)
	,95.921076,97.919027            & ! HOBr (37)
	,28.031300,29.034655            & ! C2H4 (38)
	,32.026215                      & ! CH3OH (39) (2004)
        ,93.941811,95.939764            & ! CH3Br (40) (2008) 
        ,41.026549                      & ! CH3CN (41) (2008)
        ,87.993616/                       ! CF4 (42)   (2008)	
   
	DATA MOLECULE/'VODKA'/
! --------------------------------------------------------------------------- !
!*...................................Physical state parameters:
! -------------------------------------------------------------------------- !
	IF(VS /= VS1) THEN   ! ! #### DLT8
		DLT8=10.0 ; VS=VS1 ; VF=VS+DLT8 ; VS_SEARCH =VS-OBR25 ;  VF_SEARCH = VF + OBR25
	ENDIF

	IF(MOLECULE /= GASCOMP)THEN
		MOLECULE=GASCOMP ;	MOTYPE=3
		IF(MOLECULE == 'CO2')MOTYPE=2
		IF(MOLECULE == 'H2O')MOTYPE=1
	ENDIF

	IF(P /= P1 .OR. T /= T1 .OR. RO /= RO1) THEN
 		IF(T /= T1)THEN
  			T=T1 ; TT=T ; TST = TS/T ; EK = BET*(T-TS)/(T*TS)
			NTAB_G=(T-20.0)/2+1 ; T_G1=NTAB_G*2.0+18. ; C_G2=(T-T_G1)/2. ; C_G1=1.-C_G2    
!###  296K -> N_TAB=139 ###TIPS=C_G1*QofT(N_MOLIS,NTAB_G)+C_G2*QofT(N_MOLIS,NTAB_G+1)
 		END IF
		P=P1
		RO=RO1
		PSELF=RO*10./2.6872E25*T/273.15
		PFOREI=P-PSELF
		APF = PFOREI/PS
		APS = PSELF/PS
		APF2 = PFOREI/P ; APS2=PSELF/P
	ENDIF
!*
!*...................................Line parameters:

!*....................................Line-by-line cicle:
	IF(IGAS == 1)I=1
 910	I = I + 1
	SLSS=S(I)
	IF(SLSS == 0.)GOTO 911
	           APALF = APF*ALFA(I)
               APALS = APS*ALFAS(I)
!* Lorentz half width *
	AL=TST**FACT(I)*(APALF*APF2+APALS*APS2)
		ALAL=AL*AL
!* Doppler parameter ((half width)/((ln2)**0.5) *
		VI = V(I)+SHIFT(I)*P  ! Shift
		ISO=NISO(I)
          N_MOLIS=NISO(I)/100  ! ##### Prover
 IF(ISO==222)ISO=221 ! Attention:  the SECOND N2 isotop is treated as the FIRST!!		IF(ISO /= ISOOLD .OR. T /= TOLD)THEN
		IF(ISO /= ISOOLD .OR. T /= TOLD)THEN
		ISOOLD=ISO
		TOOLD=T
		DOPCON=4.29E-7*SQRT(T/WISO(N_MOLIS))
!* Total internal partition function *
STS3=C_G1*QofT(N_MOLIS,NTAB_G)+C_G2*QofT(N_MOLIS,NTAB_G+1)
STR3=QofT(N_MOLIS,139) ! AT 296K.
QDQ=STR3/STS3*RO/PI
END IF
		ADD=DOPCON*VI
		ALAD=AL/ADD
		SL=SLSS*QDQ
		IF(STR3 >= 0.)SL=SL*EXP(EK*E(I))


!========================================================!
! LTE:	Pure radiation  (van Vleck-Weisskopf-Huber) factor!!! 26 Feb.,2009  !!!	*
! b=C2=1.438786 the second radiation constant 
!  Fi(V)=S*[1/FACTOR(VI)]*FACTOR(V)] =
!  Sref*(...)*[(1+exp(-Vi*b/T)/[Vi*(1-exp(-Vi*b/Tref))] x       ! see below  
! *** ATTENTON=> there are T in the NUMENATOR    and  Tref in the DENUMERATOR *** LTE-case !!!     
!       x [(1-exp(-V*b/T)/(1+exp(-V*b/T))]*V        ! see K_COEFF
! <<< in this non-LTE program only T is used (not referenced !!!)
        BETVITS=BET*VI/TS
        BETVIT=BET*VI/T
		IF(BETVITS > 1E-5)THEN
        EXPVV=EXP(-BETVIT)
        EXPVVS=EXP(-BETVITS)
		FVVHSL=(1.+EXPVV)/(1.-EXPVVS)/VI
		ELSE
        FVVHSL=(2.+BETVIT)/BETVITS/VI   ! attention -> bila pliuha!
		END IF
		SL=FVVHSL*SL
		SLL=SL*AL
!========================================================!
!* ------------------- LBL - calculation ---------------- *
		IF(ALAD > BOUNDL) THEN
			IF(VI < VS)THEN
			CALL LEFTLBL(VS,VI,VV_LOR,E_S)
			ELSE
								IF(VI >= VF)THEN
			CALL RIGHLBL(VS,VI,VV_LOR,E_S)
								ELSE
			CALL CENTLBL(VS,VI,VV_LOR,E_S)
								END IF
			END IF
		ELSE
		IF(ALAD > BOUNDD) THEN
			IF(VI < VS)THEN
			CALL LEFTLBL(VS,VI,VOIGT,E_S)
			ELSE
								IF(VI >= VF)THEN
			CALL RIGHLBL(VS,VI,VOIGT,E_S)
								ELSE
			CALL CENTLBL(VS,VI,VOIGT,E_S)
								END IF
			END IF
		ELSE
			IF(VI < VS)THEN
			CALL LEFTLBL(VS,VI,DOPLER,E_S)
			ELSE
								IF(VI >= VF)THEN
			CALL RIGHLBL(VS,VI,DOPLER,E_S)
								ELSE
			CALL CENTLBL(VS,VI,DOPLER,E_S)
								END IF
			END IF
		END IF
		END IF
!*....................................End of line-by-line cicle:
		GOTO 910
 911			I=I+1
				END
! ======================================================== !

!**** LEFTLBL+CENTLBL+RIGHTLBL subroutines *****
	SUBROUTINE LEFTLBL(FREQ,UL,FSHAPE,EPS)
!********************************************************
!* THIS PROGRAM INTENDS FOR THE VOLUME ABSORPTION COEF- *
!* FICIENT CALCULATIONS AT THE ATMOSPHERIC LEVELS ( Z ) *
!* "DELTA" = 10/1 [CM**(-1)] IF " Z " < 20/100 [KM] *
!********************************************************
	IMPLICIT INTEGER*4 (I-N)
	PARAMETER (NT0=10,	&
	NT1=NT0*2,NT2=NT1*2,NT3=NT2*2,NT4=NT3*2,NT5= NT4*2,	&
	NT6=NT5*2,NT7=NT6*2,NT8=NT7*2,NT9=NT8*2,NT= NT9*4+1)
	REAL*8 FREQ,UL
	common/MESH/DELTA,H,H0,H1,H2,H3,H4,H5,H6,H7,H8,H9,RK(NT)	&
	,RK0(NT0),RK0L(NT0),RK0P(NT0),RK1(NT1),RK1L(NT1),RK1P(NT1)	&
	,RK2(NT2),RK2L(NT2),RK2P(NT2),RK3(NT3),RK3L(NT3),RK3P(NT3)	&
	,RK4(NT4),RK4L(NT4),RK4P(NT4),RK5(NT5),RK5L(NT5),RK5P(NT5)	&
	,RK6(NT6),RK6L(NT6),RK6P(NT6),RK7(NT7),RK7L(NT7),RK7P(NT7)	&
	,RK8(NT8),RK8L(NT8),RK8P(NT8),RK9(NT9),RK9L(NT9),RK9P(NT9)
!* ---------------------------------------------
!*
	UU=UL-FREQ

					IF(UU >= 0.)GOTO 30
	FF=FSHAPE(UU)
		IF(FF < EPS)		RETURN
	RK(1)=FF+RK(1)
		IF(-UU < H0)		GOTO 20
	RK0P(1)=RK0P(1)+FF
	FF=FSHAPE(UU-H1)
	RK0(1)=RK0(1)+FF
	FF=FSHAPE(UU-H0)
	RK0L(1)=RK0L(1)+FF
								GOTO 10
 20	IF(-UU < H1)			GOTO 21
	RK1P(1)=RK1P(1)+FF
	FF=FSHAPE(UU-H2)
	RK1(1)=RK1(1)+FF
	FF=FSHAPE(UU-H1)
	RK1L(1)=RK1L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 121
!* 21 IF(-UU < H2.and.a < h2)goto 22
 21	IF(-UU < H2)			GOTO 22
	RK2P(1)=RK2P(1)+FF
	FF=FSHAPE(UU-H3)
	RK2(1)=RK2(1)+FF
	FF=FSHAPE(UU-H2)
	RK2L(1)=RK2L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 122
 22	IF(-UU < H3)			GOTO 23
	RK3P(1)=RK3P(1)+FF
	FF=FSHAPE(UU-H4)
	RK3(1)=RK3(1)+FF
	FF=FSHAPE(UU-H3)
	RK3L(1)=RK3L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 123
 23	IF(-UU < H4)			GOTO 24
	RK4P(1)=RK4P(1)+FF
	FF=FSHAPE(UU-H5)
	RK4(1)=RK4(1)+FF
	FF=FSHAPE(UU-H4)
	RK4L(1)=RK4L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 124
 24	IF(-UU < H5)			GOTO 25
	RK5P(1)=RK5P(1)+FF
	FF=FSHAPE(UU-H6)
	RK5(1)=RK5(1)+FF
	FF=FSHAPE(UU-H5)
	RK5L(1)=RK5L(1)+FF
	IF(FF < EPS)			RETURN
								GOTO 125
 25	IF(-UU < H6)			GOTO 26
	RK6P(1)=RK6P(1)+FF
	FF=FSHAPE(UU-H7)
	RK6(1)=RK6(1)+FF
	FF=FSHAPE(UU-H6)
	RK6L(1)=RK6L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 126
 26	IF(-UU < H7)			GOTO 27
	RK7P(1)=RK7P(1)+FF
	FF=FSHAPE(UU-H8)
	RK7(1)=RK7(1)+FF
	FF=FSHAPE(UU-H7)
	RK7L(1)=RK7L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 127
 27	IF(-UU < H8)			GOTO 28
	RK8P(1)=RK8P(1)+FF
	FF=FSHAPE(UU-H9)
	RK8(1)=RK8(1)+FF
	FF=FSHAPE(UU-H8)
	RK8L(1)=RK8L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 128
 28	IF(-UU < H9)			GOTO 29
	RK9P(1)=RK9P(1)+FF
	FF=FSHAPE(UU-H-H)
	RK9(1)=RK9(1)+FF
	FF=FSHAPE(UU-H9)
	RK9L(1)=RK9L(1)+FF
		IF(FF < EPS)			RETURN
								GOTO 129
 29	RK(2)=RK(2)+FSHAPE(UU-H)
	RK(3)=RK(3)+FSHAPE(UU-H-H)
	RK(4)=RK(4)+FSHAPE(UU+H-H9)
	FF=FSHAPE(UU-H9)
!* CASCADE
	RK(5)=RK(5)+FF
 129 RK9P(2)=RK9P(2)+FF
	FF=FSHAPE(UU-H9-H-H)
	RK9(2)=RK9(2)+FF
	FF=FSHAPE(UU-H8)
	RK9L(2)=RK9L(2)+FF
 128 RK8P(2)=RK8P(2)+FF
	FF=FSHAPE(UU-H8-H9)
	RK8(2)=RK8(2)+FF
	FF=FSHAPE(UU-H7)
	RK8L(2)=RK8L(2)+FF
	IF(FF < EPS)		RETURN
 127 RK7P(2)=RK7P(2)+FF
	FF=FSHAPE(UU-H7-H8)
	RK7(2)=RK7(2)+FF
	FF=FSHAPE(UU-H6)
	RK7L(2)=RK7L(2)+FF
	IF(FF < EPS)		RETURN
 126 RK6P(2)=RK6P(2)+FF
	FF=FSHAPE(UU-H6-H7)
	RK6(2)=RK6(2)+FF
	FF=FSHAPE(UU-H5)
	RK6L(2)=RK6L(2)+FF
	IF(FF < EPS)		RETURN
 125 RK5P(2)=RK5P(2)+FF
	FF=FSHAPE(UU-H5-H6)
	RK5(2)=RK5(2)+FF
	FF=FSHAPE(UU-H4)
	RK5L(2)=RK5L(2)+FF
	IF(FF < EPS)		RETURN
 124 RK4P(2)=RK4P(2)+FF
	FF=FSHAPE(UU-H4-H5)
	RK4(2)=RK4(2)+FF
	FF=FSHAPE(UU-H3)
	RK4L(2)=RK4L(2)+FF
	IF(FF < EPS)		RETURN
 123 RK3P(2)=RK3P(2)+FF
	FF=FSHAPE(UU-H3-H4)
	RK3(2)=RK3(2)+FF
	FF=FSHAPE(UU-H2)
	RK3L(2)=RK3L(2)+FF
	IF(FF < EPS)		RETURN
 122 RK2P(2)=RK2P(2)+FF
	FF=FSHAPE(UU-H2-H3)
	RK2(2)=RK2(2)+FF
	FF=FSHAPE(UU-H1)
	RK2L(2)=RK2L(2)+FF
	IF(FF < EPS)		RETURN
 121 RK1P(2)=RK1P(2)+FF
	FF=FSHAPE(UU-H1-H2)
	RK1(2)=RK1(2)+FF
	FF=FSHAPE(UU-H0)
	RK1L(2)=RK1L(2)+FF
!*					*** FINISH LEFT SIDE ***
 10			CONTINUE
			XXX=H0
			DO I = 2,NT0
			RK0P(I)=RK0P(I)+FF
			FF=FSHAPE(UU-XXX-H1)
			RK0(I)=RK0(I)+FF
			XXX=XXX+H0
			FF=FSHAPE(UU-XXX)
			RK0L(I)=RK0L(I)+FF
			IF(FF < EPS) GOTO 500
			END DO
 30		CONTINUE
 500		RETURN
			END
	SUBROUTINE CENTLBL(FREQ,UL,FSHAPE,EPS)
!*					***	CENTRAL PART	***
	IMPLICIT INTEGER*4 (I-N)
	PARAMETER (NT0=10,	&
	NT1=NT0*2,NT2=NT1*2,NT3=NT2*2,NT4=NT3*2,NT5= NT4*2,	&
	NT6=NT5*2,NT7=NT6*2,NT8=NT7*2,NT9=NT8*2,NT= NT9*4+1)
	REAL*8 FREQ,UL
	common/MESH/DELTA,H,H0,H1,H2,H3,H4,H5,H6,H7,H8,H9,RK(NT)	&
	,RK0(NT0),RK0L(NT0),RK0P(NT0),RK1(NT1),RK1L(NT1),RK1P(NT1)	&
	,RK2(NT2),RK2L(NT2),RK2P(NT2),RK3(NT3),RK3L(NT3),RK3P(NT3)	&
	,RK4(NT4),RK4L(NT4),RK4P(NT4),RK5(NT5),RK5L(NT5),RK5P(NT5)	&
	,RK6(NT6),RK6L(NT6),RK6P(NT6),RK7(NT7),RK7L(NT7),RK7P(NT7)	&
	,RK8(NT8),RK8L(NT8),RK8P(NT8),RK9(NT9),RK9L(NT9),RK9P(NT9)
	UU=UL-FREQ
		IF(UU >= DELTA)		GOTO 40
	FF=FSHAPE(0.)
		IF(FF < EPS)		RETURN
!*							* left-wright side *
		NPOINT=1
		CONSER=UU-H
		FA=FSHAPE(UU)
		EPS4=EPS*0.25
		IF(FA > EPS4) RK(1)=RK(1)+FA
		IF(UU < H)			GOTO 211
		I=0
		UUU=UU
		IF(UUU < H0+H0)	GOTO 201
				DO I =1,NT0
				UUU=UUU-H0
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 190
				RK0P(I)=RK0P(I)+FA
				RK0(I)=RK0(I)+FSHAPE(UUU+H1)
				RK0L(I)=RK0L(I)+FF
 190			FA=FF
				IF(UUU-H0 < H0) GOTO 201
				END DO
 201			I=I*2
				IF(UUU < H0)	GOTO 202
				IB=I+1
				DO I =IB,NT1
				UUU=UUU-H1
				FF=FSHAPE(UUU)
				IF(FF < EPS) GOTO 191
				RK1P(I)=RK1P(I)+FA
				RK1(I)=RK1(I)+FSHAPE(UUU+H2)
				RK1L(I)=RK1L(I)+FF
 191			FA=FF
				IF(UUU-H1 < H1)GOTO 202
				END DO
 202			I=I*2
				IF(UUU < H1)		GOTO 203
				IB=I+1
				DO I =IB,NT2
				UUU=UUU-H2
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 192
				RK2P(I)=RK2P(I)+FA
				RK2(I)=RK2(I)+FSHAPE(UUU+H3)
				RK2L(I)=RK2L(I)+FF
 192			FA=FF
				IF(UUU-H2 < H2) GOTO 203
				END DO
 203			I=I*2
				IF(UUU < H2)	GOTO 204
				IB=I+1
				DO I =IB,NT3
				UUU=UUU-H3
				FF=FSHAPE(UUU)
				IF(FF < EPS) GOTO 193
				RK3P(I)=RK3P(I)+FA
				RK3(I)=RK3(I)+FSHAPE(UUU+H4)
				RK3L(I)=RK3L(I)+FF
 193			FA=FF
				IF(UUU-H3 < H3)GOTO 204
				END DO
 204			I=I*2
				IF(UUU < H3) GOTO 205
				IB=I+1
				DO I =IB,NT4
				UUU=UUU-H4
				FF=FSHAPE(UUU)
				IF(FF < EPS) GOTO 194
				RK4P(I)=RK4P(I)+FA
				RK4(I)=RK4(I)+FSHAPE(UUU+H5)
				RK4L(I)=RK4L(I)+FF
 194			FA=FF
				IF(UUU-H4 < H4) GOTO 205
				END DO
 205			I=I*2
				IF(UUU < H4)	GOTO 206
				IB=I+1
				DO I =IB,NT5
				UUU=UUU-H5
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 195
				RK5P(I)=RK5P(I)+FA
				RK5(I)=RK5(I)+FSHAPE(UUU+H6)
				RK5L(I)=RK5L(I)+FF
 195			FA=FF
				IF(UUU-H5 < H5) GOTO 206
				END DO
 206			I=I*2
				IF(UUU < H5)	GOTO 207
				IB=I+1
				DO I =IB,NT6
				UUU=UUU-H6
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 196
				RK6P(I)=RK6P(I)+FA
				RK6(I)=RK6(I)+FSHAPE(UUU+H7)
				RK6L(I)=RK6L(I)+FF
 196			FA=FF
				IF(UUU-H6 < H6) GOTO 207
				END DO
 207			I=I*2
				IF(UUU < H6)	GOTO 208
				IB=I+1
				DO I =IB,NT7
				UUU=UUU-H7
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 197
				RK7P(I)=RK7P(I)+FA
				RK7(I)=RK7(I)+FSHAPE(UUU+H8)
				RK7L(I)=RK7L(I)+FF
 197			FA=FF
				IF(UUU-H7 < H7) GOTO 208
				END DO
 208			I=I*2
				IF(UUU < H7)	GOTO 209
				IB=I+1
				DO I =IB,NT8
				UUU=UUU-H8
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 198
				RK8P(I)=RK8P(I)+FA
				RK8(I)=RK8(I)+FSHAPE(UUU+H9)
				RK8L(I)=RK8L(I)+FF
 198			FA=FF
				IF(UUU-H8 < H8) GOTO 209
				END DO
 209			I=I*2
				IF(UUU < H8)	GOTO 210
				IB=I+1
				DO I =IB,NT9
				UUU=UUU-H9
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 199
				RK9P(I)=RK9P(I)+FA
				RK9(I)=RK9(I)+FSHAPE(UUU+H+H)
				RK9L(I)=RK9L(I)+FF
 199			FA=FF
				IF(UUU-H9 < H9) GOTO 210
				END DO
 210			I=I*4
				IB=I+2
				CONSER=UU-(IB-1)*H
				DO ICON=IB,NT
				RK(ICON)=RK(ICON)+FSHAPE(CONSER)
				CONSER=CONSER-H
				IF(CONSER < 0.) GOTO 490
				END DO
									GOTO 11
 490				NPOINT=ICON
!*							* wright-LEFT side *
 211	NPOINT=NPOINT+1
		UUU=DELTA-UU
		FA=FSHAPE(UUU)
				III=0
			IF(UUU < H0+H0)		GOTO 301
				DO I =NT0,1,-1
				III=III+1
				UUU=UUU-H0
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 290
				RK0L(I)=RK0L(I)+FA
				RK0(I)=RK0(I)+FSHAPE(UUU+H1)
				RK0P(I)=RK0P(I)+FF
 290			FA=FF
				IF(UUU-H0 < H0) GOTO 1301
				END DO
 301			IF(UUU < H0)	GOTO 302
 1301			III=III*2
				IB=NT1-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H1
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 291
				RK1L(I)=RK1L(I)+FA
				RK1(I)=RK1(I)+FSHAPE(UUU+H2)
				RK1P(I)=RK1P(I)+FF
 291			FA=FF
				IF(UUU-H1 < H1) GOTO 1302
				END DO
 302	IF(UUU < H1)			GOTO 303
 1302			III=III*2
				IB=NT2-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H2
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 292
				RK2L(I)=RK2L(I)+FA
				RK2(I)=RK2(I)+FSHAPE(UUU+H3)
				RK2P(I)=RK2P(I)+FF
 292			FA=FF
				IF(UUU-H2 < H2) GOTO 1303
				END DO
 303			IF(UUU < H2)	GOTO 304
 1303	III=III*2
				IB=NT3-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H3
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 293
				RK3L(I)=RK3L(I)+FA
				RK3(I)=RK3(I)+FSHAPE(UUU+H4)
				RK3P(I)=RK3P(I)+FF
 293			FA=FF
				IF(UUU-H3 < H3) GOTO 1304
				END DO
 304			IF(UUU < H3)	GOTO 305
 1304			III=III*2
				IB=NT4-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H4
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 294
				RK4L(I)=RK4L(I)+FA
				RK4(I)=RK4(I)+FSHAPE(UUU+H5)
				RK4P(I)=RK4P(I)+FF
 294			FA=FF
					IF(UUU-H4 < H4) GOTO 1305
				END DO
 305			IF(UUU < H4)	GOTO 306
 1305			III=III*2
				IB=NT5-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H5
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 295
				RK5L(I)=RK5L(I)+FA
				RK5(I)=RK5(I)+FSHAPE(UUU+H6)
				RK5P(I)=RK5P(I)+FF
 295			FA=FF
				IF(UUU-H5 < H5) GOTO 1306
				END DO
 306			IF(UUU < H5)	GOTO 307
 1306			III=III*2
				IB=NT6-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H6
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 296
				RK6L(I)=RK6L(I)+FA
				RK6(I)=RK6(I)+FSHAPE(UUU+H7)
				RK6P(I)=RK6P(I)+FF
 296			FA=FF
				IF(UUU-H6 < H6) GOTO 1307
				END DO
 307			IF(UUU < H6)	GOTO 308
 1307			III=III*2
				IB=NT7-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H7
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 297
				RK7L(I)=RK7L(I)+FA
				RK7(I)=RK7(I)+FSHAPE(UUU+H8)
				RK7P(I)=RK7P(I)+FF
 297			FA=FF
				IF(UUU-H7 < H7) GOTO 1308
				END DO
 308			IF(UUU < H7)	GOTO 309
 1308			III=III*2
				IB=NT8-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H8
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 298
				RK8L(I)=RK8L(I)+FA
				RK8(I)=RK8(I)+FSHAPE(UUU+H9)
				RK8P(I)=RK8P(I)+FF
 298			FA=FF
				IF(UUU-H8 < H8) GOTO 1309
				END DO
 309			IF(UUU < H8)	GOTO 310
 1309			III=III*2
				IB=NT9-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H9
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 299
				RK9L(I)=RK9L(I)+FA
				RK9(I)=RK9(I)+FSHAPE(UUU+H+H)
				RK9P(I)=RK9P(I)+FF
 299			FA=FF
				IF(UUU-H9 < H9) GOTO 310
				END DO
 310			III=III*4
				I=NT-III
				DO 9 II=NPOINT,I
				RK(II)=RK(II)+FSHAPE(CONSER)
 9			CONSER=CONSER-H
!*						*** FINISH CENTRAL PART ***
 40					CONTINUE
 11							RETURN
					END
			SUBROUTINE RIGHLBL(FREQ,UL,FSHAPE,EPS)
!*					*** RIGHT SIDE	***
	IMPLICIT INTEGER*4 (I-N)
	PARAMETER (NT0=10,	&
	NT1=NT0*2,NT2=NT1*2,NT3=NT2*2,NT4=NT3*2,NT5= NT4*2,	&
	NT6=NT5*2,NT7=NT6*2,NT8=NT7*2,NT9=NT8*2,NT= NT9*4+1)
	REAL*8 FREQ,UL
	common/MESH/DELTA,H,H0,H1,H2,H3,H4,H5,H6,H7,H8,H9,RK(NT)	&
	,RK0(NT0),RK0L(NT0),RK0P(NT0),RK1(NT1),RK1L(NT1),RK1P(NT1)	&
	,RK2(NT2),RK2L(NT2),RK2P(NT2),RK3(NT3),RK3L(NT3),RK3P(NT3)	&
	,RK4(NT4),RK4L(NT4),RK4P(NT4),RK5(NT5),RK5L(NT5),RK5P(NT5)	&
	,RK6(NT6),RK6L(NT6),RK6P(NT6),RK7(NT7),RK7L(NT7),RK7P(NT7)	&
	,RK8(NT8),RK8L(NT8),RK8P(NT8),RK9(NT9),RK9L(NT9),RK9P(NT9)
	UU=UL-FREQ-DELTA ! CHANGE DIAP
					IF(UU >= 25.) RETURN
	FF=FSHAPE(UU)
		IF(FF < EPS)			RETURN
		IF(UU < H0)			GOTO 60
	RK0L(NT0)=RK0L(NT0)+FF
	FF=FSHAPE(UU+H1)
	RK0(NT0)=RK0(NT0)+FF
	FF=FSHAPE(UU+H0)
	RK0P(NT0)=RK0P(NT0)+FF
								GOTO 12
 60	IF(UU < H1)			GOTO 61
	RK1L(NT1)=RK1L(NT1)+FF
	FF=FSHAPE(UU+H2)
	RK1(NT1)=RK1(NT1)+FF
	FF=FSHAPE(UU+H1)
	RK1P(NT1)=RK1P(NT1)+FF
		IF(FF < EPS)			RETURN
								GOTO 131
!* 61	IF(-UU < H2.and.a < h2)goto 22
 61	IF(UU < H2)			GOTO 62
	RK2L(NT2)=RK2L(NT2)+FF
	FF=FSHAPE(UU+H3)
	RK2(NT2)=RK2(NT2)+FF
	FF=FSHAPE(UU+H2)
	RK2P(NT2)=RK2P(NT2)+FF
		IF(FF < EPS)			RETURN
								GOTO 132
 62	IF(UU < H3)			GOTO 63
	RK3L(NT3)=RK3L(NT3)+FF
	FF=FSHAPE(UU+H4)
	RK3(NT3)=RK3(NT3)+FF
	FF=FSHAPE(UU+H3)
	RK3P(NT3)=RK3P(NT3)+FF
		IF(FF < EPS)			RETURN
								GOTO 133
 63	IF(UU < H4)			GOTO 64
	RK4L(NT4)=RK4L(NT4)+FF
	FF=FSHAPE(UU+H5)
	RK4(NT4)=RK4(NT4)+FF
	FF=FSHAPE(UU+H4)
	RK4P(NT4)=RK4P(NT4)+FF
		IF(FF < EPS)			RETURN
								GOTO 134
 64	IF(UU < H5)			GOTO 65
	RK5L(NT5)=RK5L(NT5)+FF
	FF=FSHAPE(UU+H6)
	RK5(NT5)=RK5(NT5)+FF
	FF=FSHAPE(UU+H5)
	RK5P(NT5)=RK5P(NT5)+FF
		IF(FF < EPS)			RETURN
								GOTO 135
 65	IF(UU < H6)			GOTO 66
	RK6L(NT6)=RK6L(NT6)+FF
	FF=FSHAPE(UU+H7)
	RK6(NT6)=RK6(NT6)+FF
	FF=FSHAPE(UU+H6)
	RK6P(NT6)=RK6P(NT6)+FF
		IF(FF < EPS)			RETURN
								GOTO 136
 66	IF(UU < H7)			GOTO 67
	RK7L(NT7)=RK7L(NT7)+FF
	FF=FSHAPE(UU+H8)
	RK7(NT7)=RK7(NT7)+FF
	FF=FSHAPE(UU+H7)
	RK7P(NT7)=RK7P(NT7)+FF
		IF(FF < EPS)			RETURN
								GOTO 137
 67	IF(UU < H8)			GOTO 68
	RK8L(NT8)=RK8L(NT8)+FF
	FF=FSHAPE(UU+H9)
	RK8(NT8)=RK8(NT8)+FF
	FF=FSHAPE(UU+H8)
	RK8P(NT8)=RK8P(NT8)+FF
		IF(FF < EPS)			RETURN
								GOTO 138
 68	IF(UU < H9)			GOTO 69
	RK9L(NT9)=RK9L(NT9)+FF
	FF=FSHAPE(UU+H+H)
	RK9(NT9)=RK9(NT9)+FF
	FF=FSHAPE(UU+H9)
	RK9P(NT9)=RK9P(NT9)+FF
		IF(FF < EPS)			RETURN
								GOTO 139
 69	RK(NT)=RK(NT)+FSHAPE(UU)
	RK(NT-1)=RK(NT-1)+FSHAPE(UU+H)
	RK(NT-2)=RK(NT-2)+FSHAPE(UU+H+H)
	RK(NT-3)=RK(NT-3)+FSHAPE(UU+H9-H)
!* CASCADE
 139	N2=NT9-1
	FF=FSHAPE(UU+H9)
	RK9L(N2)=RK9L(N2)+FF
	FF=FSHAPE(UU+H9+H+H)
	RK9(N2)=RK9(N2)+FF
	FF=FSHAPE(UU+H8)
	RK9P(N2)=RK9P(N2)+FF
 138 N=NT8-1
	RK8L(N)=RK8L(N)+FF
	FF=FSHAPE(UU+H8+H9)
	RK8(N)=RK8(N)+FF
	FF=FSHAPE(UU+H7)
	RK8P(N)=RK8P(N)+FF
		IF(FF < EPS)			RETURN
 137 N=NT7-1
	RK7L(N)=RK7L(N)+FF
	FF=FSHAPE(UU+H7+H8)
	RK7(N)=RK7(N)+FF
	FF=FSHAPE(UU+H6)
	RK7P(N)=RK7P(N)+FF
		IF(FF < EPS)			RETURN
 136 N=NT6-1
	RK6L(N)=RK6L(N)+FF
	FF=FSHAPE(UU+H6+H7)
	RK6(N)=RK6(N)+FF
	FF=FSHAPE(UU+H5)
	RK6P(N)=RK6P(N)+FF
		IF(FF < EPS)			RETURN
 135 N=NT5-1
	RK5L(N)=RK5L(N)+FF
	FF=FSHAPE(UU+H5+H6)
	RK5(N)=RK5(N)+FF
	FF=FSHAPE(UU+H4)
	RK5P(N)=RK5P(N)+FF
		IF(FF < EPS)			RETURN
 134 N=NT4-1
	RK4L(N)=RK4L(N)+FF
	FF=FSHAPE(UU+H4+H5)
	RK4(N)=RK4(N)+FF
	FF=FSHAPE(UU+H3)
	RK4P(N)=RK4P(N)+FF
		IF(FF < EPS)			RETURN
 133 N=NT3-1
	RK3L(N)=RK3L(N)+FF
	FF=FSHAPE(UU+H3+H4)
	RK3(N)=RK3(N)+FF
	FF=FSHAPE(UU+H2)
	RK3P(N)=RK3P(N)+FF
		IF(FF < EPS)			RETURN
 132 N=NT2-1
	RK2L(N)=RK2L(N)+FF
	FF=FSHAPE(UU+H2+H3)
	RK2(N)=RK2(N)+FF
	FF=FSHAPE(UU+H1)
	RK2P(N)=RK2P(N)+FF
		IF(FF < EPS)			RETURN
 131 N=NT1-1
	RK1L(N)=RK1L(N)+FF
	FF=FSHAPE(UU+H1+H2)
	RK1(N)=RK1(N)+FF
	FF=FSHAPE(UU+H0)
	RK1P(N)=RK1P(N)+FF
!*
!*						*** FINISH RIGHT PART ***
 12		CONTINUE
		XXX=H0
		DO I = NT0-1,1,-1
		RK0L(I)=RK0L(I)+FF
		FF=FSHAPE(UU+XXX+H1)
		RK0(I)=RK0(I)+FF
		XXX=XXX+H0
		FF=FSHAPE(UU+XXX)
		RK0P(I)=RK0P(I)+FF
		IF(FF < EPS)GOTO 1000
		END DO
			RK(1)=RK(1)+FF
 1000							RETURN
		END





! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

!*** CATAL2016+ SEARC2016+ BLOCK DATA INIT2016,ST_SUM2016,DBASE_2016 ***
	SUBROUTINE CATAL2016(VS,MOLECULE,NGAS,OBR,OBR25)
	REAL*8 V,VS,VS_SEARCH,VF_SEARCH,VS1,VS2,VV3,VV2
	CHARACTER*5 MOLECULE(*)
	PARAMETER (NLIN=150000) !
	common/PARLIO/V(NLIN),S(NLIN),ALFA(NLIN),ALFAS(NLIN)	&
	,E(NLIN),FACT(NLIN),NISO(NLIN),SHIFT(NLIN)
!* --------------------------------------------------------------
	VS1=VS-OBR25
	VS2=VS+OBR25+OBR
				I=0
			DO M=1,NGAS
				I=I+1
			S(I)=-M
			V(I)=1D15
			CALL SEARC2016(MOLECULE(M),VS1,VS2,I)
			I=I+1
			S(I)=0.
			V(I)=1D15
			END DO
			I=I+1
			S(I)=-1E10
			V(I)=1D15
			END
	SUBROUTINE SEARC2016(MOLECULE, VSTART, VFINISH,I )
!*	<<< ODINARY PRECISION FOR ' S(I) '>>>
	REAL*8 V, VV, VVV, VSTART, VFINISH
	INTEGER*4 NZ,NS,N1,N2,NISO
!*
	CHARACTER*(*) MOLECULE
	CHARACTER*5	COMPONENT
	CHARACTER*17 FNAME
!*.....For UNIX
	PARAMETER ( NMOL = 49, NLIN = 150000, IDENT = 126) ! Line 33
!*
	common/PARLIO/ V(NLIN), S(NLIN), ALFA(NLIN),ALFAS(NLIN),E(NLIN)	&
					,FACT(NLIN),NISO(NLIN),SHIFT(NLIN)
	common/GASES/ COMPONENT(NMOL)
	common/FILES/FNAME(NMOL), NZ(NMOL)
		

		DO 100 K = 1, NMOL
!*
	IF(MOLECULE == COMPONENT(K))THEN
		OPEN(IDENT, ACCESS='DIRECT', FORM='UNFORMATTED',	&
			RECL=36, FILE="/srv/HITRAN16/"//FNAME(K))! 36 for OTHER FORTRANs
	NZCOMP	= K
	GO TO 10
								ENDIF
!*
100		CONTINUE
!*
	WRITE(*,*)'There are no ',MOLECULE,' lines in HITRAN'
!*...............................Searching for origin line record :
10	N1 = 1
	N2 = NZ(NZCOMP)
!*
1	NS = ( N1 + N2 )/2
	READ(IDENT, REC = NS ) VV
		IF(NS==N1) GOTO 3
		IF(VV < VSTART)GOTO 2
	N2 = NS
			GO TO 1
2	N1 = NS
			GO TO 1
3	NN = NS+1
!*...............................Searching for ending line record :
	N1 = NN
	N2 = NZ(NZCOMP)
!*
4	NS = ( N1 + N2 )/2
!*
	READ(IDENT, REC = NS ) VV
			IF(NS == N1)GOTO 6
		IF(VV < VFINISH)GOTO 5
	N2 = NS
			GO TO 4
5	N1 = NS
			GO TO 4
6	NF = NS
!*
!************** MAIN PART : READING *********************
		READ(IDENT, REC = NN ) VV
		READ(IDENT, REC = NF ) VVV
	IF(VV < VSTART .OR. VVV > VFINISH)THEN ! There are no lines in interval :
	CLOSE(IDENT)
	RETURN
	END IF
!*
	DO 61 J = NN, NF
!*
	I = I + 1
		IF(I >= NLIN)THEN
	WRITE(*,*)'<<< TOO MUCH LINES: ',I,NLIN,'>>>'
	WRITE(*,*)	&
	'(par. NLIN in SEARC2000,CATAL2016,LBL2016 should be increased)'
		PAUSE
		END IF
READ(IDENT,REC=J)V(I),S(I),ALFA(I),ALFAS(I),E(I),FACT(I),NISO(I),SHIFT(I)
		IF(V(I) < VSTART .OR. V(I) > VFINISH)THEN
		write(*,*)' attention',i,NN,NF,VV,VVV,vstart,v(i),vfinish
		pause
		END IF
	IF(ALFAS(I) < 0.0001)THEN
		MH20=NISO(I)/10
		IF(MH20 == 1)THEN
			ALFAS(I)=5.*ALFA(I)
					ELSE
			ALFAS(I)=ALFA(I)
		END IF
		END IF
!*	write(*,*)V(I),S(I),ALFA(I),ALFAS(I),E(I),FACT(I),NISO(I),shift
!*		pause
!*
61	CONTINUE
!*
	CLOSE(IDENT)
!*..............................End of the program
	RETURN
	END
!*****************************
	BLOCK DATA INIT2016 ! for HITRAN-2016
! *** ATTENTION: at present only 12 gases from HITRAN-2016.!!! *** ! 
	PARAMETER ( NMOL = 49) !49 molecules in HITRAN considered
	CHARACTER*5 COMPONENT  !names of molecules
	CHARACTER*17 FNAME     !names of files in HITRAN16 folder
!*
	common/GASES/ COMPONENT(NMOL)
	common/FILES/FNAME(NMOL), NZ(NMOL)

	DATA COMPONENT	&
		/'H2O', 'CO2', 'O3' , 'N2O' , 'CO' , 'CH4' , 'O2',	&
			'NO ', 'SO2', 'NO2' , 'NH3' , 'HNO3', 'OH' , 'HF',	&
			'HCL', 'HBR', 'HI' , 'CLO' , 'OCS' , 'H2CO', 'HOCL',	&
			'N2' , 'HCN', 'CH3CL', 'H2O2', 'C2H2', 'C2H6', 'PH3',	&

	'COF2','SF6','H2S','HCOOH','HO2','O','ClONO2','NO+','HOBr','C2H4','CH3OH', &
	'CH3Br','CH3CN','CF4', &
	 'C4H2','HC3N','H2','CS','SO3','C2N2','COCL2'/
!*.........................Corresponding line parameters file names:
       DATA  FNAME   &
  /'H16.01','H16.02','H16.03','H16.04','H16.05','H16.06','H16.07','H16.08', &
   'H16.09','H16.10','H16.11','H16.12', &
   'H12.13','H12.14','H12.15','H12.16',  &
   'H12.17','H12.18','H12.19','H12.20','H12.21','H12.22','H12.23','H12.24',   &
   'H12.25','H12.26','H12.27','H12.28','H12.29','H12.30','H12.31','H12.32',    &
   'H12.33','H12.34','H12.35','H12.36','H12.37','H12.38','H12.39', &
   'H12.40','H12.41','H12.42',                                                 &
   'H16.43','H16.44','H16.45','H16.46','H16.47','H16.48','H16.49' / !Fictiv

!*.........................Corresponding numbers of lines in files :
      DATA  NZ &
!*...............       H2O( 1):    CO2( 2):     O3( 3):    N2O( 4): ! HITRAN-2016
! 2008        /        69201,       314919,       409686,     47843,  & ! All the CO2 Lines
! 2012        /       224515,       471847,       422116,     47843,  & ! All the CO2 Lines 
              /       304225,       559874,       449570,    160287, &

!*...............         CO( 5):    CH4( 6):     O2( 7):     NO( 8):
! 2008                     4477,      290091,       6428,      105079,  &
! 2012                     4606,      468013,      13975,      105079,  &
                           5381,      450332,      14085,      105079,  &
!*...............       SO2( 9):    NO2(10):    NH3(11):   HNO3(12):
! 2008                     58250,     104223,      29084,     487254,  &
! 2012                     95121,     104223,      46392,     961962,  &
                           95121,     104223,      67148,    1008972,  &
!*...............         OH(13):     HF(14):    HCL(15):    HBR(16):
! 2008                     31976,        107,        613,       1293,  &
                           31979,      34376,      83691,       8980,  &

!*...............          HJ(17):    CLO(18):    OCS(19):   H2CO(20):
! 2008                      806,       11501,       29242,       37050,  &
                           4751,       11501,       29361,       44601,  &

!*...............         HOCL(21):     N2(22):    HCN(23):  CH3CL(24):
! 2008                    16276,        120,        4253,       196171,   &
                          16276,       1268,        4253,       212496,   &

!*...............           H2O2(25):   C2H2(26):   C2H6(27):    PH3(28):
! 2008                     126983,       11340,       22402,      20099,   &
                           126983,       20410,       49629,      22189,   &

!*...............           COF2(29):   SF6(30):    H2S(31)   HCOOH(32):
! 2008                    70601,      2889065,      20788,       62684,   &
                         184104,      2889065,      54235,       62684,   &

!*...............        HO2(33):     O(34):    ClONO2(35):  NO+(36):
! 2008                    38804,        2,        32199,     1206,     &
                          38804,        2,        32199,     1206,     &

!*...............        HOBr(37):    C2H4(38):  CH3OH(39)
! 2008                     4358,      18378,      19899,               &
                           4358,      18378,      19899,               &

!* 2008                                  CH3Br      CH3CN        CF4
!                       					36911,       3572,      60033/
                                         36911,       3572,      60033, &
  7*0/  ! <=== at present NOT ready!

      END
 


BLOCK DATA ST_SUM2017 ! Parameterization by R.R. Gamache
!*-------------------------------------------------------*
	common/ST_SUM16/QofT(74,501) ! HITRAN-2017  (74 for the first 12 gases)
! -------------------------------------------------------------------- !
 ! *** for Venus_19.07.2022 (up to 1000K). TIPS - Gamache 2017 (program by Fomin, 28 June, 2018). ***
 ! ***  N_MOLIS - total number of each isotop in my database ***
  !            1           1
 DATA(QofT(           1 ,J),J=1,501)/   3.34890007972717       &
, 0.38767E+01, 0.44281E+01, 0.49993E+01, 0.55878E+01, 0.61918E+01, 0.68103E+01, 0.74424E+01, 0.80877E+01, 0.87459E+01, 0.94170E+01 &
, 0.10101E+02, 0.10797E+02, 0.11506E+02, 0.12227E+02, 0.12962E+02, 0.13708E+02, 0.14467E+02, 0.15239E+02, 0.16023E+02, 0.16819E+02 &
, 0.17627E+02, 0.18447E+02, 0.19280E+02, 0.20124E+02, 0.20980E+02, 0.21847E+02, 0.22726E+02, 0.23617E+02, 0.24519E+02, 0.25432E+02 &
, 0.26356E+02, 0.27291E+02, 0.28237E+02, 0.29194E+02, 0.30161E+02, 0.31139E+02, 0.32127E+02, 0.33125E+02, 0.34134E+02, 0.35153E+02 &
, 0.36182E+02, 0.37221E+02, 0.38270E+02, 0.39328E+02, 0.40396E+02, 0.41474E+02, 0.42561E+02, 0.43658E+02, 0.44764E+02, 0.45880E+02 &
, 0.47004E+02, 0.48138E+02, 0.49281E+02, 0.50432E+02, 0.51593E+02, 0.52763E+02, 0.53941E+02, 0.55128E+02, 0.56324E+02, 0.57528E+02 &
, 0.58741E+02, 0.59963E+02, 0.61193E+02, 0.62431E+02, 0.63678E+02, 0.64932E+02, 0.66196E+02, 0.67467E+02, 0.68746E+02, 0.70034E+02 &
, 0.71330E+02, 0.72633E+02, 0.73945E+02, 0.75264E+02, 0.76591E+02, 0.77926E+02, 0.79269E+02, 0.80620E+02, 0.81978E+02, 0.83344E+02 &
, 0.84718E+02, 0.86099E+02, 0.87487E+02, 0.88883E+02, 0.90287E+02, 0.91698E+02, 0.93116E+02, 0.94542E+02, 0.95975E+02, 0.97415E+02 &
, 0.98863E+02, 0.10032E+03, 0.10178E+03, 0.10325E+03, 0.10472E+03, 0.10621E+03, 0.10770E+03, 0.10920E+03, 0.11070E+03, 0.11221E+03 &
, 0.11373E+03, 0.11525E+03, 0.11679E+03, 0.11833E+03, 0.11987E+03, 0.12142E+03, 0.12298E+03, 0.12455E+03, 0.12612E+03, 0.12770E+03 &
, 0.12929E+03, 0.13088E+03, 0.13248E+03, 0.13409E+03, 0.13570E+03, 0.13732E+03, 0.13895E+03, 0.14058E+03, 0.14222E+03, 0.14386E+03 &
, 0.14552E+03, 0.14717E+03, 0.14884E+03, 0.15051E+03, 0.15219E+03, 0.15387E+03, 0.15556E+03, 0.15726E+03, 0.15896E+03, 0.16067E+03 &
, 0.16239E+03, 0.16411E+03, 0.16584E+03, 0.16758E+03, 0.16932E+03, 0.17107E+03, 0.17282E+03, 0.17458E+03, 0.17635E+03, 0.17812E+03 &
, 0.17990E+03, 0.18168E+03, 0.18348E+03, 0.18527E+03, 0.18708E+03, 0.18889E+03, 0.19070E+03, 0.19253E+03, 0.19435E+03, 0.19619E+03 &
, 0.19803E+03, 0.19988E+03, 0.20173E+03, 0.20359E+03, 0.20545E+03, 0.20733E+03, 0.20920E+03, 0.21109E+03, 0.21298E+03, 0.21487E+03 &
, 0.21677E+03, 0.21868E+03, 0.22060E+03, 0.22252E+03, 0.22444E+03, 0.22637E+03, 0.22831E+03, 0.23026E+03, 0.23221E+03, 0.23416E+03 &
, 0.23613E+03, 0.23810E+03, 0.24007E+03, 0.24205E+03, 0.24404E+03, 0.24603E+03, 0.24803E+03, 0.25003E+03, 0.25205E+03, 0.25406E+03 &
, 0.25609E+03, 0.25812E+03, 0.26015E+03, 0.26219E+03, 0.26424E+03, 0.26629E+03, 0.26835E+03, 0.27042E+03, 0.27249E+03, 0.27457E+03 &
, 0.27665E+03, 0.27874E+03, 0.28084E+03, 0.28294E+03, 0.28505E+03, 0.28716E+03, 0.28929E+03, 0.29141E+03, 0.29354E+03, 0.29568E+03 &
, 0.29783E+03, 0.29998E+03, 0.30214E+03, 0.30430E+03, 0.30647E+03, 0.30865E+03, 0.31083E+03, 0.31302E+03, 0.31521E+03, 0.31741E+03 &
, 0.31962E+03, 0.32183E+03, 0.32405E+03, 0.32627E+03, 0.32851E+03, 0.33074E+03, 0.33299E+03, 0.33524E+03, 0.33749E+03, 0.33976E+03 &
, 0.34203E+03, 0.34430E+03, 0.34658E+03, 0.34887E+03, 0.35116E+03, 0.35346E+03, 0.35577E+03, 0.35808E+03, 0.36040E+03, 0.36273E+03 &
, 0.36506E+03, 0.36740E+03, 0.36974E+03, 0.37209E+03, 0.37445E+03, 0.37681E+03, 0.37918E+03, 0.38156E+03, 0.38394E+03, 0.38633E+03 &
, 0.38872E+03, 0.39113E+03, 0.39353E+03, 0.39595E+03, 0.39837E+03, 0.40080E+03, 0.40323E+03, 0.40567E+03, 0.40812E+03, 0.41057E+03 &
, 0.41303E+03, 0.41550E+03, 0.41797E+03, 0.42045E+03, 0.42294E+03, 0.42543E+03, 0.42793E+03, 0.43043E+03, 0.43295E+03, 0.43546E+03 &
, 0.43799E+03, 0.44052E+03, 0.44306E+03, 0.44561E+03, 0.44816E+03, 0.45072E+03, 0.45328E+03, 0.45585E+03, 0.45843E+03, 0.46102E+03 &
, 0.46361E+03, 0.46621E+03, 0.46881E+03, 0.47143E+03, 0.47405E+03, 0.47667E+03, 0.47930E+03, 0.48194E+03, 0.48459E+03, 0.48724E+03 &
, 0.48990E+03, 0.49257E+03, 0.49524E+03, 0.49792E+03, 0.50061E+03, 0.50331E+03, 0.50601E+03, 0.50871E+03, 0.51143E+03, 0.51415E+03 &
, 0.51688E+03, 0.51962E+03, 0.52236E+03, 0.52511E+03, 0.52787E+03, 0.53063E+03, 0.53340E+03, 0.53618E+03, 0.53896E+03, 0.54175E+03 &
, 0.54455E+03, 0.54736E+03, 0.55017E+03, 0.55299E+03, 0.55582E+03, 0.55866E+03, 0.56150E+03, 0.56435E+03, 0.56720E+03, 0.57007E+03 &
, 0.57294E+03, 0.57581E+03, 0.57870E+03, 0.58159E+03, 0.58449E+03, 0.58740E+03, 0.59031E+03, 0.59323E+03, 0.59616E+03, 0.59910E+03 &
, 0.60204E+03, 0.60499E+03, 0.60795E+03, 0.61091E+03, 0.61389E+03, 0.61687E+03, 0.61986E+03, 0.62285E+03, 0.62585E+03, 0.62886E+03 &
, 0.63188E+03, 0.63491E+03, 0.63794E+03, 0.64098E+03, 0.64403E+03, 0.64708E+03, 0.65014E+03, 0.65322E+03, 0.65629E+03, 0.65938E+03 &
, 0.66247E+03, 0.66557E+03, 0.66868E+03, 0.67180E+03, 0.67492E+03, 0.67805E+03, 0.68119E+03, 0.68434E+03, 0.68749E+03, 0.69066E+03 &
, 0.69383E+03, 0.69701E+03, 0.70019E+03, 0.70339E+03, 0.70659E+03, 0.70980E+03, 0.71301E+03, 0.71624E+03, 0.71947E+03, 0.72271E+03 &
, 0.72596E+03, 0.72922E+03, 0.73249E+03, 0.73576E+03, 0.73904E+03, 0.74233E+03, 0.74563E+03, 0.74893E+03, 0.75224E+03, 0.75557E+03 &
, 0.75890E+03, 0.76223E+03, 0.76558E+03, 0.76893E+03, 0.77229E+03, 0.77566E+03, 0.77904E+03, 0.78243E+03, 0.78582E+03, 0.78923E+03 &
, 0.79264E+03, 0.79606E+03, 0.79949E+03, 0.80292E+03, 0.80637E+03, 0.80982E+03, 0.81328E+03, 0.81675E+03, 0.82023E+03, 0.82372E+03 &
, 0.82721E+03, 0.83071E+03, 0.83423E+03, 0.83775E+03, 0.84128E+03, 0.84481E+03, 0.84836E+03, 0.85191E+03, 0.85548E+03, 0.85905E+03 &
, 0.86263E+03, 0.86622E+03, 0.86982E+03, 0.87342E+03, 0.87704E+03, 0.88066E+03, 0.88429E+03, 0.88793E+03, 0.89158E+03, 0.89524E+03 &
, 0.89891E+03, 0.90259E+03, 0.90627E+03, 0.90997E+03, 0.91367E+03, 0.91738E+03, 0.92110E+03, 0.92483E+03, 0.92857E+03, 0.93232E+03 &
, 0.93607E+03, 0.93984E+03, 0.94361E+03, 0.94740E+03, 0.95119E+03, 0.95499E+03, 0.95880E+03, 0.96262E+03, 0.96645E+03, 0.97029E+03 &
, 0.97414E+03, 0.97799E+03, 0.98186E+03, 0.98573E+03, 0.98962E+03, 0.99351E+03, 0.99742E+03, 0.10013E+04, 0.10052E+04, 0.10092E+04 &
, 0.10131E+04, 0.10171E+04, 0.10210E+04, 0.10250E+04, 0.10290E+04, 0.10330E+04, 0.10370E+04, 0.10410E+04, 0.10450E+04, 0.10490E+04 &
, 0.10530E+04, 0.10571E+04, 0.10611E+04, 0.10652E+04, 0.10693E+04, 0.10734E+04, 0.10774E+04, 0.10816E+04, 0.10857E+04, 0.10898E+04 &
, 0.10939E+04, 0.10981E+04, 0.11022E+04, 0.11064E+04, 0.11105E+04, 0.11147E+04, 0.11189E+04, 0.11231E+04, 0.11273E+04, 0.11315E+04 &
, 0.11358E+04, 0.11400E+04, 0.11443E+04, 0.11485E+04, 0.11528E+04, 0.11571E+04, 0.11614E+04, 0.11657E+04, 0.11700E+04, 0.11743E+04 &
, 0.11786E+04, 0.11830E+04, 0.11873E+04, 0.11917E+04, 0.11961E+04, 0.12004E+04, 0.12048E+04, 0.12092E+04, 0.12136E+04, 0.12181E+04 &
, 0.12225E+04, 0.12269E+04, 0.12314E+04, 0.12359E+04, 0.12403E+04, 0.12448E+04, 0.12493E+04, 0.12538E+04, 0.12584E+04, 0.12629E+04 /
 !            1           2
 DATA(QofT(           2 ,J),J=1,501)/   3.37280011177063       &
, 0.39054E+01, 0.44617E+01, 0.50378E+01, 0.56313E+01, 0.62404E+01, 0.68640E+01, 0.75014E+01, 0.81520E+01, 0.88157E+01, 0.94923E+01 &
, 0.10182E+02, 0.10884E+02, 0.11599E+02, 0.12326E+02, 0.13066E+02, 0.13819E+02, 0.14585E+02, 0.15362E+02, 0.16153E+02, 0.16956E+02 &
, 0.17771E+02, 0.18598E+02, 0.19437E+02, 0.20288E+02, 0.21151E+02, 0.22026E+02, 0.22913E+02, 0.23811E+02, 0.24720E+02, 0.25641E+02 &
, 0.26573E+02, 0.27516E+02, 0.28470E+02, 0.29434E+02, 0.30410E+02, 0.31396E+02, 0.32392E+02, 0.33399E+02, 0.34417E+02, 0.35444E+02 &
, 0.36482E+02, 0.37529E+02, 0.38587E+02, 0.39654E+02, 0.40731E+02, 0.41818E+02, 0.42915E+02, 0.44021E+02, 0.45136E+02, 0.46261E+02 &
, 0.47395E+02, 0.48538E+02, 0.49691E+02, 0.50852E+02, 0.52023E+02, 0.53202E+02, 0.54390E+02, 0.55587E+02, 0.56793E+02, 0.58008E+02 &
, 0.59231E+02, 0.60463E+02, 0.61703E+02, 0.62952E+02, 0.64209E+02, 0.65475E+02, 0.66748E+02, 0.68030E+02, 0.69321E+02, 0.70619E+02 &
, 0.71926E+02, 0.73240E+02, 0.74563E+02, 0.75893E+02, 0.77232E+02, 0.78578E+02, 0.79932E+02, 0.81294E+02, 0.82664E+02, 0.84042E+02 &
, 0.85427E+02, 0.86819E+02, 0.88220E+02, 0.89628E+02, 0.91043E+02, 0.92466E+02, 0.93896E+02, 0.95334E+02, 0.96779E+02, 0.98232E+02 &
, 0.99691E+02, 0.10116E+03, 0.10263E+03, 0.10411E+03, 0.10560E+03, 0.10710E+03, 0.10860E+03, 0.11011E+03, 0.11163E+03, 0.11315E+03 &
, 0.11468E+03, 0.11622E+03, 0.11777E+03, 0.11932E+03, 0.12088E+03, 0.12244E+03, 0.12402E+03, 0.12560E+03, 0.12718E+03, 0.12877E+03 &
, 0.13037E+03, 0.13198E+03, 0.13359E+03, 0.13521E+03, 0.13684E+03, 0.13847E+03, 0.14011E+03, 0.14176E+03, 0.14341E+03, 0.14507E+03 &
, 0.14674E+03, 0.14841E+03, 0.15009E+03, 0.15178E+03, 0.15347E+03, 0.15517E+03, 0.15687E+03, 0.15858E+03, 0.16030E+03, 0.16203E+03 &
, 0.16376E+03, 0.16550E+03, 0.16724E+03, 0.16899E+03, 0.17075E+03, 0.17251E+03, 0.17428E+03, 0.17605E+03, 0.17783E+03, 0.17962E+03 &
, 0.18142E+03, 0.18322E+03, 0.18502E+03, 0.18684E+03, 0.18866E+03, 0.19048E+03, 0.19231E+03, 0.19415E+03, 0.19599E+03, 0.19784E+03 &
, 0.19970E+03, 0.20156E+03, 0.20343E+03, 0.20531E+03, 0.20719E+03, 0.20908E+03, 0.21097E+03, 0.21287E+03, 0.21478E+03, 0.21669E+03 &
, 0.21861E+03, 0.22053E+03, 0.22246E+03, 0.22440E+03, 0.22634E+03, 0.22829E+03, 0.23024E+03, 0.23221E+03, 0.23417E+03, 0.23615E+03 &
, 0.23813E+03, 0.24011E+03, 0.24210E+03, 0.24410E+03, 0.24611E+03, 0.24812E+03, 0.25013E+03, 0.25215E+03, 0.25418E+03, 0.25622E+03 &
, 0.25826E+03, 0.26031E+03, 0.26236E+03, 0.26442E+03, 0.26648E+03, 0.26856E+03, 0.27063E+03, 0.27272E+03, 0.27481E+03, 0.27690E+03 &
, 0.27900E+03, 0.28111E+03, 0.28323E+03, 0.28535E+03, 0.28748E+03, 0.28961E+03, 0.29175E+03, 0.29389E+03, 0.29604E+03, 0.29820E+03 &
, 0.30037E+03, 0.30254E+03, 0.30471E+03, 0.30690E+03, 0.30908E+03, 0.31128E+03, 0.31348E+03, 0.31569E+03, 0.31790E+03, 0.32012E+03 &
, 0.32235E+03, 0.32458E+03, 0.32682E+03, 0.32906E+03, 0.33131E+03, 0.33357E+03, 0.33584E+03, 0.33811E+03, 0.34038E+03, 0.34266E+03 &
, 0.34495E+03, 0.34725E+03, 0.34955E+03, 0.35186E+03, 0.35417E+03, 0.35649E+03, 0.35882E+03, 0.36115E+03, 0.36349E+03, 0.36584E+03 &
, 0.36819E+03, 0.37055E+03, 0.37292E+03, 0.37529E+03, 0.37767E+03, 0.38005E+03, 0.38244E+03, 0.38484E+03, 0.38724E+03, 0.38965E+03 &
, 0.39207E+03, 0.39449E+03, 0.39692E+03, 0.39936E+03, 0.40180E+03, 0.40425E+03, 0.40671E+03, 0.40917E+03, 0.41164E+03, 0.41412E+03 &
, 0.41660E+03, 0.41909E+03, 0.42158E+03, 0.42409E+03, 0.42659E+03, 0.42911E+03, 0.43163E+03, 0.43416E+03, 0.43669E+03, 0.43924E+03 &
, 0.44179E+03, 0.44434E+03, 0.44690E+03, 0.44947E+03, 0.45205E+03, 0.45463E+03, 0.45722E+03, 0.45981E+03, 0.46242E+03, 0.46502E+03 &
, 0.46764E+03, 0.47026E+03, 0.47289E+03, 0.47553E+03, 0.47817E+03, 0.48082E+03, 0.48348E+03, 0.48614E+03, 0.48881E+03, 0.49149E+03 &
, 0.49418E+03, 0.49687E+03, 0.49957E+03, 0.50227E+03, 0.50498E+03, 0.50770E+03, 0.51043E+03, 0.51316E+03, 0.51590E+03, 0.51865E+03 &
, 0.52141E+03, 0.52417E+03, 0.52693E+03, 0.52971E+03, 0.53249E+03, 0.53528E+03, 0.53808E+03, 0.54088E+03, 0.54369E+03, 0.54651E+03 &
, 0.54934E+03, 0.55217E+03, 0.55501E+03, 0.55786E+03, 0.56071E+03, 0.56357E+03, 0.56644E+03, 0.56932E+03, 0.57220E+03, 0.57509E+03 &
, 0.57799E+03, 0.58089E+03, 0.58380E+03, 0.58672E+03, 0.58965E+03, 0.59259E+03, 0.59553E+03, 0.59848E+03, 0.60143E+03, 0.60440E+03 &
, 0.60737E+03, 0.61035E+03, 0.61333E+03, 0.61633E+03, 0.61933E+03, 0.62234E+03, 0.62535E+03, 0.62838E+03, 0.63141E+03, 0.63445E+03 &
, 0.63749E+03, 0.64055E+03, 0.64361E+03, 0.64668E+03, 0.64975E+03, 0.65284E+03, 0.65593E+03, 0.65903E+03, 0.66214E+03, 0.66525E+03 &
, 0.66838E+03, 0.67151E+03, 0.67465E+03, 0.67779E+03, 0.68095E+03, 0.68411E+03, 0.68728E+03, 0.69045E+03, 0.69364E+03, 0.69683E+03 &
, 0.70003E+03, 0.70324E+03, 0.70646E+03, 0.70968E+03, 0.71292E+03, 0.71616E+03, 0.71940E+03, 0.72266E+03, 0.72593E+03, 0.72920E+03 &
, 0.73248E+03, 0.73577E+03, 0.73906E+03, 0.74237E+03, 0.74568E+03, 0.74900E+03, 0.75233E+03, 0.75567E+03, 0.75901E+03, 0.76237E+03 &
, 0.76573E+03, 0.76910E+03, 0.77248E+03, 0.77586E+03, 0.77926E+03, 0.78266E+03, 0.78607E+03, 0.78949E+03, 0.79292E+03, 0.79636E+03 &
, 0.79980E+03, 0.80325E+03, 0.80671E+03, 0.81018E+03, 0.81366E+03, 0.81715E+03, 0.82064E+03, 0.82415E+03, 0.82766E+03, 0.83118E+03 &
, 0.83471E+03, 0.83824E+03, 0.84179E+03, 0.84535E+03, 0.84891E+03, 0.85248E+03, 0.85606E+03, 0.85965E+03, 0.86325E+03, 0.86685E+03 &
, 0.87047E+03, 0.87409E+03, 0.87773E+03, 0.88137E+03, 0.88502E+03, 0.88868E+03, 0.89234E+03, 0.89602E+03, 0.89970E+03, 0.90340E+03 &
, 0.90710E+03, 0.91081E+03, 0.91453E+03, 0.91826E+03, 0.92200E+03, 0.92575E+03, 0.92951E+03, 0.93327E+03, 0.93705E+03, 0.94083E+03 &
, 0.94462E+03, 0.94843E+03, 0.95224E+03, 0.95606E+03, 0.95989E+03, 0.96372E+03, 0.96757E+03, 0.97143E+03, 0.97529E+03, 0.97917E+03 &
, 0.98305E+03, 0.98695E+03, 0.99085E+03, 0.99476E+03, 0.99868E+03, 0.10026E+04, 0.10066E+04, 0.10105E+04, 0.10145E+04, 0.10184E+04 &
, 0.10224E+04, 0.10264E+04, 0.10304E+04, 0.10344E+04, 0.10384E+04, 0.10424E+04, 0.10465E+04, 0.10505E+04, 0.10546E+04, 0.10586E+04 &
, 0.10627E+04, 0.10668E+04, 0.10709E+04, 0.10750E+04, 0.10791E+04, 0.10832E+04, 0.10874E+04, 0.10915E+04, 0.10956E+04, 0.10998E+04 &
, 0.11040E+04, 0.11082E+04, 0.11123E+04, 0.11165E+04, 0.11208E+04, 0.11250E+04, 0.11292E+04, 0.11334E+04, 0.11377E+04, 0.11420E+04 &
, 0.11462E+04, 0.11505E+04, 0.11548E+04, 0.11591E+04, 0.11634E+04, 0.11677E+04, 0.11721E+04, 0.11764E+04, 0.11807E+04, 0.11851E+04 &
, 0.11895E+04, 0.11939E+04, 0.11982E+04, 0.12026E+04, 0.12071E+04, 0.12115E+04, 0.12159E+04, 0.12204E+04, 0.12248E+04, 0.12293E+04 &
, 0.12337E+04, 0.12382E+04, 0.12427E+04, 0.12472E+04, 0.12518E+04, 0.12563E+04, 0.12608E+04, 0.12654E+04, 0.12699E+04, 0.12745E+04 /
 !            1           3
 DATA(QofT(           3 ,J),J=1,501)/   20.1690006256104       &
, 0.23351E+02, 0.26675E+02, 0.30118E+02, 0.33665E+02, 0.37305E+02, 0.41032E+02, 0.44841E+02, 0.48730E+02, 0.52697E+02, 0.56740E+02 &
, 0.60860E+02, 0.65057E+02, 0.69329E+02, 0.73677E+02, 0.78101E+02, 0.82600E+02, 0.87175E+02, 0.91824E+02, 0.96548E+02, 0.10135E+03 &
, 0.10622E+03, 0.11116E+03, 0.11618E+03, 0.12126E+03, 0.12642E+03, 0.13165E+03, 0.13695E+03, 0.14232E+03, 0.14775E+03, 0.15325E+03 &
, 0.15882E+03, 0.16446E+03, 0.17016E+03, 0.17592E+03, 0.18175E+03, 0.18765E+03, 0.19360E+03, 0.19962E+03, 0.20570E+03, 0.21184E+03 &
, 0.21804E+03, 0.22430E+03, 0.23062E+03, 0.23700E+03, 0.24344E+03, 0.24993E+03, 0.25649E+03, 0.26310E+03, 0.26976E+03, 0.27648E+03 &
, 0.28326E+03, 0.29009E+03, 0.29698E+03, 0.30392E+03, 0.31092E+03, 0.31797E+03, 0.32507E+03, 0.33222E+03, 0.33943E+03, 0.34669E+03 &
, 0.35400E+03, 0.36136E+03, 0.36877E+03, 0.37623E+03, 0.38375E+03, 0.39131E+03, 0.39892E+03, 0.40659E+03, 0.41430E+03, 0.42206E+03 &
, 0.42986E+03, 0.43772E+03, 0.44563E+03, 0.45358E+03, 0.46158E+03, 0.46962E+03, 0.47772E+03, 0.48586E+03, 0.49404E+03, 0.50227E+03 &
, 0.51055E+03, 0.51887E+03, 0.52724E+03, 0.53566E+03, 0.54412E+03, 0.55262E+03, 0.56117E+03, 0.56976E+03, 0.57840E+03, 0.58708E+03 &
, 0.59580E+03, 0.60457E+03, 0.61338E+03, 0.62223E+03, 0.63113E+03, 0.64007E+03, 0.64905E+03, 0.65807E+03, 0.66714E+03, 0.67625E+03 &
, 0.68540E+03, 0.69459E+03, 0.70382E+03, 0.71310E+03, 0.72241E+03, 0.73177E+03, 0.74117E+03, 0.75061E+03, 0.76009E+03, 0.76961E+03 &
, 0.77917E+03, 0.78877E+03, 0.79841E+03, 0.80809E+03, 0.81781E+03, 0.82757E+03, 0.83737E+03, 0.84721E+03, 0.85709E+03, 0.86701E+03 &
, 0.87697E+03, 0.88696E+03, 0.89700E+03, 0.90707E+03, 0.91719E+03, 0.92734E+03, 0.93753E+03, 0.94776E+03, 0.95802E+03, 0.96833E+03 &
, 0.97867E+03, 0.98906E+03, 0.99948E+03, 0.10099E+04, 0.10204E+04, 0.10310E+04, 0.10415E+04, 0.10521E+04, 0.10628E+04, 0.10735E+04 &
, 0.10842E+04, 0.10950E+04, 0.11058E+04, 0.11166E+04, 0.11275E+04, 0.11384E+04, 0.11493E+04, 0.11603E+04, 0.11713E+04, 0.11824E+04 &
, 0.11935E+04, 0.12046E+04, 0.12158E+04, 0.12270E+04, 0.12382E+04, 0.12495E+04, 0.12608E+04, 0.12722E+04, 0.12836E+04, 0.12950E+04 &
, 0.13064E+04, 0.13179E+04, 0.13295E+04, 0.13410E+04, 0.13527E+04, 0.13643E+04, 0.13760E+04, 0.13877E+04, 0.13995E+04, 0.14113E+04 &
, 0.14231E+04, 0.14349E+04, 0.14469E+04, 0.14588E+04, 0.14708E+04, 0.14828E+04, 0.14948E+04, 0.15069E+04, 0.15190E+04, 0.15312E+04 &
, 0.15434E+04, 0.15556E+04, 0.15679E+04, 0.15802E+04, 0.15925E+04, 0.16049E+04, 0.16173E+04, 0.16298E+04, 0.16423E+04, 0.16548E+04 &
, 0.16673E+04, 0.16799E+04, 0.16926E+04, 0.17052E+04, 0.17180E+04, 0.17307E+04, 0.17435E+04, 0.17563E+04, 0.17692E+04, 0.17820E+04 &
, 0.17950E+04, 0.18079E+04, 0.18209E+04, 0.18340E+04, 0.18471E+04, 0.18602E+04, 0.18733E+04, 0.18865E+04, 0.18997E+04, 0.19130E+04 &
, 0.19263E+04, 0.19396E+04, 0.19530E+04, 0.19664E+04, 0.19799E+04, 0.19933E+04, 0.20069E+04, 0.20204E+04, 0.20340E+04, 0.20477E+04 &
, 0.20613E+04, 0.20751E+04, 0.20888E+04, 0.21026E+04, 0.21164E+04, 0.21303E+04, 0.21442E+04, 0.21581E+04, 0.21721E+04, 0.21861E+04 &
, 0.22002E+04, 0.22142E+04, 0.22284E+04, 0.22425E+04, 0.22567E+04, 0.22710E+04, 0.22853E+04, 0.22996E+04, 0.23139E+04, 0.23283E+04 &
, 0.23428E+04, 0.23573E+04, 0.23718E+04, 0.23863E+04, 0.24009E+04, 0.24155E+04, 0.24302E+04, 0.24449E+04, 0.24597E+04, 0.24744E+04 &
, 0.24893E+04, 0.25041E+04, 0.25190E+04, 0.25340E+04, 0.25490E+04, 0.25640E+04, 0.25790E+04, 0.25941E+04, 0.26093E+04, 0.26245E+04 &
, 0.26397E+04, 0.26549E+04, 0.26702E+04, 0.26856E+04, 0.27010E+04, 0.27164E+04, 0.27318E+04, 0.27473E+04, 0.27629E+04, 0.27784E+04 &
, 0.27941E+04, 0.28097E+04, 0.28254E+04, 0.28412E+04, 0.28569E+04, 0.28728E+04, 0.28886E+04, 0.29045E+04, 0.29205E+04, 0.29365E+04 &
, 0.29525E+04, 0.29686E+04, 0.29847E+04, 0.30008E+04, 0.30170E+04, 0.30333E+04, 0.30495E+04, 0.30658E+04, 0.30822E+04, 0.30986E+04 &
, 0.31150E+04, 0.31315E+04, 0.31481E+04, 0.31646E+04, 0.31812E+04, 0.31979E+04, 0.32146E+04, 0.32313E+04, 0.32481E+04, 0.32649E+04 &
, 0.32818E+04, 0.32987E+04, 0.33156E+04, 0.33326E+04, 0.33497E+04, 0.33667E+04, 0.33839E+04, 0.34010E+04, 0.34182E+04, 0.34355E+04 &
, 0.34528E+04, 0.34701E+04, 0.34875E+04, 0.35049E+04, 0.35224E+04, 0.35399E+04, 0.35574E+04, 0.35750E+04, 0.35927E+04, 0.36104E+04 &
, 0.36281E+04, 0.36459E+04, 0.36637E+04, 0.36815E+04, 0.36995E+04, 0.37174E+04, 0.37354E+04, 0.37534E+04, 0.37715E+04, 0.37897E+04 &
, 0.38078E+04, 0.38261E+04, 0.38443E+04, 0.38626E+04, 0.38810E+04, 0.38994E+04, 0.39178E+04, 0.39363E+04, 0.39549E+04, 0.39735E+04 &
, 0.39921E+04, 0.40108E+04, 0.40295E+04, 0.40483E+04, 0.40671E+04, 0.40859E+04, 0.41048E+04, 0.41238E+04, 0.41428E+04, 0.41619E+04 &
, 0.41809E+04, 0.42001E+04, 0.42193E+04, 0.42385E+04, 0.42578E+04, 0.42771E+04, 0.42965E+04, 0.43159E+04, 0.43354E+04, 0.43549E+04 &
, 0.43745E+04, 0.43941E+04, 0.44138E+04, 0.44335E+04, 0.44532E+04, 0.44730E+04, 0.44929E+04, 0.45128E+04, 0.45327E+04, 0.45527E+04 &
, 0.45728E+04, 0.45929E+04, 0.46130E+04, 0.46332E+04, 0.46535E+04, 0.46738E+04, 0.46941E+04, 0.47145E+04, 0.47349E+04, 0.47554E+04 &
, 0.47760E+04, 0.47965E+04, 0.48172E+04, 0.48379E+04, 0.48586E+04, 0.48794E+04, 0.49002E+04, 0.49211E+04, 0.49421E+04, 0.49631E+04 &
, 0.49841E+04, 0.50052E+04, 0.50263E+04, 0.50475E+04, 0.50688E+04, 0.50901E+04, 0.51114E+04, 0.51328E+04, 0.51543E+04, 0.51758E+04 &
, 0.51973E+04, 0.52189E+04, 0.52406E+04, 0.52623E+04, 0.52840E+04, 0.53059E+04, 0.53277E+04, 0.53496E+04, 0.53716E+04, 0.53936E+04 &
, 0.54157E+04, 0.54378E+04, 0.54600E+04, 0.54823E+04, 0.55045E+04, 0.55269E+04, 0.55493E+04, 0.55717E+04, 0.55942E+04, 0.56168E+04 &
, 0.56394E+04, 0.56620E+04, 0.56848E+04, 0.57075E+04, 0.57304E+04, 0.57532E+04, 0.57762E+04, 0.57991E+04, 0.58222E+04, 0.58453E+04 &
, 0.58684E+04, 0.58916E+04, 0.59149E+04, 0.59382E+04, 0.59616E+04, 0.59850E+04, 0.60085E+04, 0.60320E+04, 0.60556E+04, 0.60793E+04 &
, 0.61030E+04, 0.61268E+04, 0.61506E+04, 0.61745E+04, 0.61984E+04, 0.62224E+04, 0.62464E+04, 0.62705E+04, 0.62947E+04, 0.63189E+04 &
, 0.63432E+04, 0.63675E+04, 0.63919E+04, 0.64163E+04, 0.64408E+04, 0.64654E+04, 0.64900E+04, 0.65147E+04, 0.65394E+04, 0.65642E+04 &
, 0.65891E+04, 0.66140E+04, 0.66389E+04, 0.66640E+04, 0.66890E+04, 0.67142E+04, 0.67394E+04, 0.67647E+04, 0.67900E+04, 0.68154E+04 &
, 0.68408E+04, 0.68663E+04, 0.68919E+04, 0.69175E+04, 0.69432E+04, 0.69689E+04, 0.69947E+04, 0.70206E+04, 0.70465E+04, 0.70725E+04 &
, 0.70985E+04, 0.71246E+04, 0.71508E+04, 0.71770E+04, 0.72033E+04, 0.72296E+04, 0.72560E+04, 0.72825E+04, 0.73090E+04, 0.73356E+04 &
, 0.73623E+04, 0.73890E+04, 0.74158E+04, 0.74426E+04, 0.74695E+04, 0.74965E+04, 0.75235E+04, 0.75506E+04, 0.75778E+04, 0.76050E+04 /
 !            1           4
 DATA(QofT(           4 ,J),J=1,501)/   17.6679992675781       &
, 0.20013E+02, 0.22479E+02, 0.25055E+02, 0.27733E+02, 0.30509E+02, 0.33377E+02, 0.36333E+02, 0.39374E+02, 0.42498E+02, 0.45702E+02 &
, 0.48985E+02, 0.52343E+02, 0.55776E+02, 0.59282E+02, 0.62860E+02, 0.66507E+02, 0.70223E+02, 0.74006E+02, 0.77856E+02, 0.81771E+02 &
, 0.85750E+02, 0.89792E+02, 0.93896E+02, 0.98061E+02, 0.10229E+03, 0.10657E+03, 0.11092E+03, 0.11532E+03, 0.11978E+03, 0.12429E+03 &
, 0.12886E+03, 0.13349E+03, 0.13817E+03, 0.14291E+03, 0.14769E+03, 0.15253E+03, 0.15743E+03, 0.16237E+03, 0.16737E+03, 0.17241E+03 &
, 0.17751E+03, 0.18266E+03, 0.18785E+03, 0.19310E+03, 0.19839E+03, 0.20373E+03, 0.20912E+03, 0.21455E+03, 0.22003E+03, 0.22556E+03 &
, 0.23114E+03, 0.23675E+03, 0.24242E+03, 0.24813E+03, 0.25388E+03, 0.25968E+03, 0.26552E+03, 0.27141E+03, 0.27734E+03, 0.28331E+03 &
, 0.28932E+03, 0.29538E+03, 0.30148E+03, 0.30762E+03, 0.31380E+03, 0.32003E+03, 0.32629E+03, 0.33260E+03, 0.33894E+03, 0.34533E+03 &
, 0.35176E+03, 0.35822E+03, 0.36473E+03, 0.37127E+03, 0.37786E+03, 0.38448E+03, 0.39114E+03, 0.39784E+03, 0.40458E+03, 0.41136E+03 &
, 0.41817E+03, 0.42503E+03, 0.43192E+03, 0.43885E+03, 0.44581E+03, 0.45281E+03, 0.45985E+03, 0.46693E+03, 0.47404E+03, 0.48119E+03 &
, 0.48837E+03, 0.49559E+03, 0.50285E+03, 0.51014E+03, 0.51747E+03, 0.52483E+03, 0.53223E+03, 0.53966E+03, 0.54713E+03, 0.55464E+03 &
, 0.56218E+03, 0.56975E+03, 0.57736E+03, 0.58500E+03, 0.59268E+03, 0.60039E+03, 0.60813E+03, 0.61591E+03, 0.62373E+03, 0.63157E+03 &
, 0.63945E+03, 0.64737E+03, 0.65532E+03, 0.66330E+03, 0.67131E+03, 0.67936E+03, 0.68745E+03, 0.69556E+03, 0.70371E+03, 0.71189E+03 &
, 0.72010E+03, 0.72835E+03, 0.73663E+03, 0.74494E+03, 0.75329E+03, 0.76167E+03, 0.77008E+03, 0.77852E+03, 0.78700E+03, 0.79551E+03 &
, 0.80405E+03, 0.81262E+03, 0.82123E+03, 0.82987E+03, 0.83854E+03, 0.84724E+03, 0.85598E+03, 0.86474E+03, 0.87354E+03, 0.88237E+03 &
, 0.89124E+03, 0.90013E+03, 0.90906E+03, 0.91802E+03, 0.92701E+03, 0.93604E+03, 0.94509E+03, 0.95418E+03, 0.96330E+03, 0.97245E+03 &
, 0.98164E+03, 0.99086E+03, 0.10001E+04, 0.10094E+04, 0.10187E+04, 0.10280E+04, 0.10374E+04, 0.10468E+04, 0.10563E+04, 0.10657E+04 &
, 0.10752E+04, 0.10848E+04, 0.10944E+04, 0.11040E+04, 0.11136E+04, 0.11232E+04, 0.11329E+04, 0.11427E+04, 0.11524E+04, 0.11622E+04 &
, 0.11721E+04, 0.11819E+04, 0.11918E+04, 0.12017E+04, 0.12117E+04, 0.12217E+04, 0.12317E+04, 0.12417E+04, 0.12518E+04, 0.12620E+04 &
, 0.12721E+04, 0.12823E+04, 0.12925E+04, 0.13027E+04, 0.13130E+04, 0.13233E+04, 0.13337E+04, 0.13441E+04, 0.13545E+04, 0.13649E+04 &
, 0.13754E+04, 0.13859E+04, 0.13965E+04, 0.14070E+04, 0.14176E+04, 0.14283E+04, 0.14390E+04, 0.14497E+04, 0.14604E+04, 0.14712E+04 &
, 0.14820E+04, 0.14929E+04, 0.15037E+04, 0.15146E+04, 0.15256E+04, 0.15366E+04, 0.15476E+04, 0.15586E+04, 0.15697E+04, 0.15808E+04 &
, 0.15920E+04, 0.16031E+04, 0.16144E+04, 0.16256E+04, 0.16369E+04, 0.16482E+04, 0.16596E+04, 0.16710E+04, 0.16824E+04, 0.16938E+04 &
, 0.17053E+04, 0.17169E+04, 0.17284E+04, 0.17400E+04, 0.17516E+04, 0.17633E+04, 0.17750E+04, 0.17867E+04, 0.17985E+04, 0.18103E+04 &
, 0.18222E+04, 0.18340E+04, 0.18460E+04, 0.18579E+04, 0.18699E+04, 0.18819E+04, 0.18940E+04, 0.19061E+04, 0.19182E+04, 0.19303E+04 &
, 0.19425E+04, 0.19548E+04, 0.19671E+04, 0.19794E+04, 0.19917E+04, 0.20041E+04, 0.20165E+04, 0.20290E+04, 0.20414E+04, 0.20540E+04 &
, 0.20665E+04, 0.20791E+04, 0.20918E+04, 0.21045E+04, 0.21172E+04, 0.21299E+04, 0.21427E+04, 0.21555E+04, 0.21684E+04, 0.21813E+04 &
, 0.21942E+04, 0.22072E+04, 0.22202E+04, 0.22333E+04, 0.22464E+04, 0.22595E+04, 0.22727E+04, 0.22859E+04, 0.22991E+04, 0.23124E+04 &
, 0.23257E+04, 0.23391E+04, 0.23525E+04, 0.23659E+04, 0.23794E+04, 0.23929E+04, 0.24064E+04, 0.24200E+04, 0.24337E+04, 0.24473E+04 &
, 0.24610E+04, 0.24748E+04, 0.24886E+04, 0.25024E+04, 0.25163E+04, 0.25302E+04, 0.25441E+04, 0.25581E+04, 0.25721E+04, 0.25862E+04 &
, 0.26003E+04, 0.26145E+04, 0.26287E+04, 0.26429E+04, 0.26572E+04, 0.26715E+04, 0.26858E+04, 0.27002E+04, 0.27147E+04, 0.27291E+04 &
, 0.27437E+04, 0.27582E+04, 0.27728E+04, 0.27875E+04, 0.28021E+04, 0.28169E+04, 0.28316E+04, 0.28465E+04, 0.28613E+04, 0.28762E+04 &
, 0.28911E+04, 0.29061E+04, 0.29211E+04, 0.29362E+04, 0.29513E+04, 0.29665E+04, 0.29817E+04, 0.29969E+04, 0.30122E+04, 0.30275E+04 &
, 0.30429E+04, 0.30583E+04, 0.30737E+04, 0.30892E+04, 0.31048E+04, 0.31204E+04, 0.31360E+04, 0.31517E+04, 0.31674E+04, 0.31831E+04 &
, 0.31990E+04, 0.32148E+04, 0.32307E+04, 0.32466E+04, 0.32626E+04, 0.32787E+04, 0.32947E+04, 0.33109E+04, 0.33270E+04, 0.33432E+04 &
, 0.33595E+04, 0.33758E+04, 0.33922E+04, 0.34085E+04, 0.34250E+04, 0.34415E+04, 0.34580E+04, 0.34746E+04, 0.34912E+04, 0.35079E+04 &
, 0.35246E+04, 0.35414E+04, 0.35582E+04, 0.35751E+04, 0.35920E+04, 0.36089E+04, 0.36259E+04, 0.36430E+04, 0.36601E+04, 0.36772E+04 &
, 0.36944E+04, 0.37116E+04, 0.37289E+04, 0.37463E+04, 0.37637E+04, 0.37811E+04, 0.37986E+04, 0.38161E+04, 0.38337E+04, 0.38513E+04 &
, 0.38690E+04, 0.38867E+04, 0.39045E+04, 0.39223E+04, 0.39402E+04, 0.39581E+04, 0.39761E+04, 0.39941E+04, 0.40122E+04, 0.40303E+04 &
, 0.40485E+04, 0.40667E+04, 0.40850E+04, 0.41033E+04, 0.41217E+04, 0.41401E+04, 0.41586E+04, 0.41772E+04, 0.41957E+04, 0.42144E+04 &
, 0.42331E+04, 0.42518E+04, 0.42706E+04, 0.42894E+04, 0.43083E+04, 0.43273E+04, 0.43462E+04, 0.43653E+04, 0.43844E+04, 0.44036E+04 &
, 0.44228E+04, 0.44420E+04, 0.44613E+04, 0.44807E+04, 0.45001E+04, 0.45196E+04, 0.45391E+04, 0.45587E+04, 0.45783E+04, 0.45980E+04 &
, 0.46178E+04, 0.46375E+04, 0.46574E+04, 0.46773E+04, 0.46973E+04, 0.47173E+04, 0.47373E+04, 0.47575E+04, 0.47776E+04, 0.47979E+04 &
, 0.48182E+04, 0.48385E+04, 0.48589E+04, 0.48794E+04, 0.48999E+04, 0.49204E+04, 0.49411E+04, 0.49618E+04, 0.49825E+04, 0.50033E+04 &
, 0.50241E+04, 0.50450E+04, 0.50660E+04, 0.50870E+04, 0.51081E+04, 0.51292E+04, 0.51504E+04, 0.51717E+04, 0.51930E+04, 0.52144E+04 &
, 0.52358E+04, 0.52573E+04, 0.52788E+04, 0.53004E+04, 0.53221E+04, 0.53438E+04, 0.53655E+04, 0.53874E+04, 0.54093E+04, 0.54312E+04 &
, 0.54532E+04, 0.54753E+04, 0.54974E+04, 0.55196E+04, 0.55419E+04, 0.55642E+04, 0.55866E+04, 0.56090E+04, 0.56315E+04, 0.56541E+04 &
, 0.56767E+04, 0.56993E+04, 0.57221E+04, 0.57449E+04, 0.57677E+04, 0.57907E+04, 0.58136E+04, 0.58367E+04, 0.58598E+04, 0.58830E+04 &
, 0.59062E+04, 0.59295E+04, 0.59528E+04, 0.59763E+04, 0.59997E+04, 0.60233E+04, 0.60469E+04, 0.60706E+04, 0.60943E+04, 0.61181E+04 &
, 0.61420E+04, 0.61659E+04, 0.61899E+04, 0.62139E+04, 0.62381E+04, 0.62622E+04, 0.62865E+04, 0.63108E+04, 0.63352E+04, 0.63596E+04 &
, 0.63841E+04, 0.64087E+04, 0.64334E+04, 0.64581E+04, 0.64828E+04, 0.65077E+04, 0.65326E+04, 0.65575E+04, 0.65826E+04, 0.66077E+04 /
 !            1           5
 DATA(QofT(           5 ,J),J=1,501)/   17.8610000610352       &
, 0.20238E+02, 0.22736E+02, 0.25345E+02, 0.28058E+02, 0.30868E+02, 0.33772E+02, 0.36765E+02, 0.39844E+02, 0.43007E+02, 0.46252E+02 &
, 0.49575E+02, 0.52975E+02, 0.56451E+02, 0.60001E+02, 0.63623E+02, 0.67315E+02, 0.71078E+02, 0.74908E+02, 0.78806E+02, 0.82770E+02 &
, 0.86798E+02, 0.90891E+02, 0.95046E+02, 0.99263E+02, 0.10354E+03, 0.10788E+03, 0.11228E+03, 0.11674E+03, 0.12125E+03, 0.12582E+03 &
, 0.13045E+03, 0.13513E+03, 0.13987E+03, 0.14467E+03, 0.14952E+03, 0.15442E+03, 0.15937E+03, 0.16438E+03, 0.16944E+03, 0.17455E+03 &
, 0.17971E+03, 0.18492E+03, 0.19018E+03, 0.19549E+03, 0.20085E+03, 0.20625E+03, 0.21171E+03, 0.21721E+03, 0.22276E+03, 0.22836E+03 &
, 0.23400E+03, 0.23969E+03, 0.24543E+03, 0.25121E+03, 0.25704E+03, 0.26291E+03, 0.26882E+03, 0.27478E+03, 0.28079E+03, 0.28683E+03 &
, 0.29292E+03, 0.29905E+03, 0.30523E+03, 0.31145E+03, 0.31771E+03, 0.32401E+03, 0.33035E+03, 0.33674E+03, 0.34316E+03, 0.34963E+03 &
, 0.35614E+03, 0.36268E+03, 0.36927E+03, 0.37590E+03, 0.38257E+03, 0.38927E+03, 0.39602E+03, 0.40280E+03, 0.40963E+03, 0.41649E+03 &
, 0.42339E+03, 0.43033E+03, 0.43730E+03, 0.44432E+03, 0.45137E+03, 0.45846E+03, 0.46559E+03, 0.47275E+03, 0.47995E+03, 0.48719E+03 &
, 0.49447E+03, 0.50178E+03, 0.50913E+03, 0.51651E+03, 0.52393E+03, 0.53139E+03, 0.53888E+03, 0.54640E+03, 0.55397E+03, 0.56157E+03 &
, 0.56920E+03, 0.57687E+03, 0.58457E+03, 0.59231E+03, 0.60008E+03, 0.60789E+03, 0.61573E+03, 0.62361E+03, 0.63152E+03, 0.63947E+03 &
, 0.64745E+03, 0.65546E+03, 0.66351E+03, 0.67159E+03, 0.67971E+03, 0.68786E+03, 0.69604E+03, 0.70426E+03, 0.71251E+03, 0.72080E+03 &
, 0.72911E+03, 0.73746E+03, 0.74585E+03, 0.75427E+03, 0.76272E+03, 0.77120E+03, 0.77972E+03, 0.78827E+03, 0.79685E+03, 0.80547E+03 &
, 0.81411E+03, 0.82280E+03, 0.83151E+03, 0.84026E+03, 0.84904E+03, 0.85785E+03, 0.86669E+03, 0.87557E+03, 0.88448E+03, 0.89343E+03 &
, 0.90240E+03, 0.91141E+03, 0.92045E+03, 0.92952E+03, 0.93863E+03, 0.94777E+03, 0.95694E+03, 0.96614E+03, 0.97538E+03, 0.98465E+03 &
, 0.99395E+03, 0.10033E+04, 0.10126E+04, 0.10220E+04, 0.10315E+04, 0.10409E+04, 0.10504E+04, 0.10600E+04, 0.10695E+04, 0.10791E+04 &
, 0.10887E+04, 0.10984E+04, 0.11081E+04, 0.11178E+04, 0.11276E+04, 0.11374E+04, 0.11472E+04, 0.11570E+04, 0.11669E+04, 0.11768E+04 &
, 0.11868E+04, 0.11968E+04, 0.12068E+04, 0.12168E+04, 0.12269E+04, 0.12370E+04, 0.12472E+04, 0.12574E+04, 0.12676E+04, 0.12778E+04 &
, 0.12881E+04, 0.12984E+04, 0.13088E+04, 0.13192E+04, 0.13296E+04, 0.13400E+04, 0.13505E+04, 0.13610E+04, 0.13716E+04, 0.13821E+04 &
, 0.13928E+04, 0.14034E+04, 0.14141E+04, 0.14248E+04, 0.14355E+04, 0.14463E+04, 0.14571E+04, 0.14680E+04, 0.14789E+04, 0.14898E+04 &
, 0.15007E+04, 0.15117E+04, 0.15228E+04, 0.15338E+04, 0.15449E+04, 0.15560E+04, 0.15672E+04, 0.15784E+04, 0.15896E+04, 0.16009E+04 &
, 0.16122E+04, 0.16235E+04, 0.16348E+04, 0.16462E+04, 0.16577E+04, 0.16692E+04, 0.16807E+04, 0.16922E+04, 0.17038E+04, 0.17154E+04 &
, 0.17270E+04, 0.17387E+04, 0.17504E+04, 0.17622E+04, 0.17740E+04, 0.17858E+04, 0.17977E+04, 0.18095E+04, 0.18215E+04, 0.18334E+04 &
, 0.18454E+04, 0.18575E+04, 0.18696E+04, 0.18817E+04, 0.18938E+04, 0.19060E+04, 0.19182E+04, 0.19305E+04, 0.19428E+04, 0.19551E+04 &
, 0.19675E+04, 0.19799E+04, 0.19923E+04, 0.20048E+04, 0.20173E+04, 0.20299E+04, 0.20425E+04, 0.20551E+04, 0.20678E+04, 0.20805E+04 &
, 0.20932E+04, 0.21060E+04, 0.21188E+04, 0.21317E+04, 0.21445E+04, 0.21575E+04, 0.21704E+04, 0.21834E+04, 0.21965E+04, 0.22096E+04 &
, 0.22227E+04, 0.22359E+04, 0.22491E+04, 0.22623E+04, 0.22756E+04, 0.22889E+04, 0.23022E+04, 0.23156E+04, 0.23291E+04, 0.23425E+04 &
, 0.23561E+04, 0.23696E+04, 0.23832E+04, 0.23968E+04, 0.24105E+04, 0.24242E+04, 0.24380E+04, 0.24518E+04, 0.24656E+04, 0.24795E+04 &
, 0.24934E+04, 0.25073E+04, 0.25213E+04, 0.25354E+04, 0.25494E+04, 0.25635E+04, 0.25777E+04, 0.25919E+04, 0.26061E+04, 0.26204E+04 &
, 0.26347E+04, 0.26491E+04, 0.26635E+04, 0.26780E+04, 0.26924E+04, 0.27070E+04, 0.27215E+04, 0.27362E+04, 0.27508E+04, 0.27655E+04 &
, 0.27802E+04, 0.27950E+04, 0.28099E+04, 0.28247E+04, 0.28396E+04, 0.28546E+04, 0.28696E+04, 0.28846E+04, 0.28997E+04, 0.29148E+04 &
, 0.29300E+04, 0.29452E+04, 0.29605E+04, 0.29758E+04, 0.29911E+04, 0.30065E+04, 0.30219E+04, 0.30374E+04, 0.30529E+04, 0.30685E+04 &
, 0.30841E+04, 0.30998E+04, 0.31155E+04, 0.31312E+04, 0.31470E+04, 0.31628E+04, 0.31787E+04, 0.31946E+04, 0.32106E+04, 0.32266E+04 &
, 0.32427E+04, 0.32588E+04, 0.32750E+04, 0.32912E+04, 0.33074E+04, 0.33237E+04, 0.33400E+04, 0.33564E+04, 0.33728E+04, 0.33893E+04 &
, 0.34058E+04, 0.34224E+04, 0.34390E+04, 0.34557E+04, 0.34724E+04, 0.34892E+04, 0.35060E+04, 0.35228E+04, 0.35397E+04, 0.35567E+04 &
, 0.35737E+04, 0.35907E+04, 0.36078E+04, 0.36250E+04, 0.36421E+04, 0.36594E+04, 0.36767E+04, 0.36940E+04, 0.37114E+04, 0.37288E+04 &
, 0.37463E+04, 0.37638E+04, 0.37814E+04, 0.37991E+04, 0.38167E+04, 0.38345E+04, 0.38522E+04, 0.38701E+04, 0.38880E+04, 0.39059E+04 &
, 0.39239E+04, 0.39419E+04, 0.39600E+04, 0.39781E+04, 0.39963E+04, 0.40145E+04, 0.40328E+04, 0.40511E+04, 0.40695E+04, 0.40880E+04 &
, 0.41065E+04, 0.41250E+04, 0.41436E+04, 0.41622E+04, 0.41809E+04, 0.41997E+04, 0.42185E+04, 0.42374E+04, 0.42563E+04, 0.42752E+04 &
, 0.42942E+04, 0.43133E+04, 0.43324E+04, 0.43516E+04, 0.43708E+04, 0.43901E+04, 0.44094E+04, 0.44288E+04, 0.44483E+04, 0.44677E+04 &
, 0.44873E+04, 0.45069E+04, 0.45266E+04, 0.45463E+04, 0.45660E+04, 0.45859E+04, 0.46057E+04, 0.46257E+04, 0.46456E+04, 0.46657E+04 &
, 0.46858E+04, 0.47059E+04, 0.47261E+04, 0.47464E+04, 0.47667E+04, 0.47871E+04, 0.48075E+04, 0.48280E+04, 0.48486E+04, 0.48692E+04 &
, 0.48898E+04, 0.49105E+04, 0.49313E+04, 0.49521E+04, 0.49730E+04, 0.49940E+04, 0.50150E+04, 0.50360E+04, 0.50571E+04, 0.50783E+04 &
, 0.50996E+04, 0.51208E+04, 0.51422E+04, 0.51636E+04, 0.51851E+04, 0.52066E+04, 0.52282E+04, 0.52498E+04, 0.52715E+04, 0.52933E+04 &
, 0.53151E+04, 0.53370E+04, 0.53589E+04, 0.53809E+04, 0.54030E+04, 0.54251E+04, 0.54473E+04, 0.54695E+04, 0.54918E+04, 0.55142E+04 &
, 0.55366E+04, 0.55591E+04, 0.55817E+04, 0.56043E+04, 0.56269E+04, 0.56497E+04, 0.56725E+04, 0.56953E+04, 0.57182E+04, 0.57412E+04 &
, 0.57643E+04, 0.57874E+04, 0.58105E+04, 0.58338E+04, 0.58571E+04, 0.58804E+04, 0.59038E+04, 0.59273E+04, 0.59509E+04, 0.59745E+04 &
, 0.59981E+04, 0.60219E+04, 0.60457E+04, 0.60695E+04, 0.60935E+04, 0.61175E+04, 0.61415E+04, 0.61656E+04, 0.61898E+04, 0.62141E+04 &
, 0.62384E+04, 0.62628E+04, 0.62872E+04, 0.63118E+04, 0.63363E+04, 0.63610E+04, 0.63857E+04, 0.64105E+04, 0.64353E+04, 0.64602E+04 &
, 0.64852E+04, 0.65103E+04, 0.65354E+04, 0.65606E+04, 0.65858E+04, 0.66111E+04, 0.66365E+04, 0.66620E+04, 0.66875E+04, 0.67131E+04 /
 !            1           6
 DATA(QofT(           6 ,J),J=1,501)/   106.620002746582       &
, 0.12079E+03, 0.13568E+03, 0.15124E+03, 0.16742E+03, 0.18419E+03, 0.20151E+03, 0.21936E+03, 0.23773E+03, 0.25660E+03, 0.27595E+03 &
, 0.29577E+03, 0.31606E+03, 0.33679E+03, 0.35796E+03, 0.37957E+03, 0.40160E+03, 0.42404E+03, 0.44689E+03, 0.47014E+03, 0.49378E+03 &
, 0.51781E+03, 0.54222E+03, 0.56701E+03, 0.59216E+03, 0.61768E+03, 0.64356E+03, 0.66980E+03, 0.69638E+03, 0.72331E+03, 0.75058E+03 &
, 0.77819E+03, 0.80613E+03, 0.83440E+03, 0.86300E+03, 0.89192E+03, 0.92116E+03, 0.95071E+03, 0.98057E+03, 0.10107E+04, 0.10412E+04 &
, 0.10720E+04, 0.11031E+04, 0.11345E+04, 0.11661E+04, 0.11981E+04, 0.12304E+04, 0.12629E+04, 0.12957E+04, 0.13288E+04, 0.13622E+04 &
, 0.13959E+04, 0.14298E+04, 0.14640E+04, 0.14985E+04, 0.15333E+04, 0.15683E+04, 0.16036E+04, 0.16391E+04, 0.16749E+04, 0.17110E+04 &
, 0.17474E+04, 0.17839E+04, 0.18208E+04, 0.18579E+04, 0.18952E+04, 0.19328E+04, 0.19707E+04, 0.20088E+04, 0.20471E+04, 0.20857E+04 &
, 0.21245E+04, 0.21636E+04, 0.22029E+04, 0.22424E+04, 0.22822E+04, 0.23222E+04, 0.23625E+04, 0.24030E+04, 0.24437E+04, 0.24846E+04 &
, 0.25258E+04, 0.25672E+04, 0.26089E+04, 0.26507E+04, 0.26928E+04, 0.27351E+04, 0.27777E+04, 0.28204E+04, 0.28634E+04, 0.29066E+04 &
, 0.29500E+04, 0.29937E+04, 0.30376E+04, 0.30816E+04, 0.31259E+04, 0.31704E+04, 0.32152E+04, 0.32601E+04, 0.33053E+04, 0.33506E+04 &
, 0.33962E+04, 0.34420E+04, 0.34880E+04, 0.35342E+04, 0.35807E+04, 0.36273E+04, 0.36741E+04, 0.37212E+04, 0.37684E+04, 0.38159E+04 &
, 0.38636E+04, 0.39114E+04, 0.39595E+04, 0.40078E+04, 0.40563E+04, 0.41050E+04, 0.41539E+04, 0.42029E+04, 0.42522E+04, 0.43017E+04 &
, 0.43514E+04, 0.44013E+04, 0.44514E+04, 0.45017E+04, 0.45522E+04, 0.46029E+04, 0.46538E+04, 0.47049E+04, 0.47562E+04, 0.48077E+04 &
, 0.48594E+04, 0.49113E+04, 0.49634E+04, 0.50157E+04, 0.50682E+04, 0.51208E+04, 0.51737E+04, 0.52268E+04, 0.52801E+04, 0.53335E+04 &
, 0.53872E+04, 0.54411E+04, 0.54951E+04, 0.55494E+04, 0.56038E+04, 0.56585E+04, 0.57133E+04, 0.57683E+04, 0.58236E+04, 0.58790E+04 &
, 0.59346E+04, 0.59904E+04, 0.60464E+04, 0.61026E+04, 0.61590E+04, 0.62156E+04, 0.62724E+04, 0.63294E+04, 0.63866E+04, 0.64440E+04 &
, 0.65016E+04, 0.65594E+04, 0.66173E+04, 0.66755E+04, 0.67339E+04, 0.67924E+04, 0.68512E+04, 0.69101E+04, 0.69693E+04, 0.70286E+04 &
, 0.70882E+04, 0.71479E+04, 0.72079E+04, 0.72680E+04, 0.73283E+04, 0.73889E+04, 0.74496E+04, 0.75105E+04, 0.75717E+04, 0.76330E+04 &
, 0.76945E+04, 0.77563E+04, 0.78182E+04, 0.78803E+04, 0.79426E+04, 0.80052E+04, 0.80679E+04, 0.81308E+04, 0.81940E+04, 0.82573E+04 &
, 0.83208E+04, 0.83846E+04, 0.84485E+04, 0.85126E+04, 0.85770E+04, 0.86415E+04, 0.87063E+04, 0.87712E+04, 0.88364E+04, 0.89017E+04 &
, 0.89673E+04, 0.90331E+04, 0.90991E+04, 0.91652E+04, 0.92316E+04, 0.92982E+04, 0.93650E+04, 0.94320E+04, 0.94993E+04, 0.95667E+04 &
, 0.96343E+04, 0.97022E+04, 0.97702E+04, 0.98385E+04, 0.99070E+04, 0.99756E+04, 0.10045E+05, 0.10114E+05, 0.10183E+05, 0.10252E+05 &
, 0.10322E+05, 0.10392E+05, 0.10462E+05, 0.10533E+05, 0.10603E+05, 0.10674E+05, 0.10745E+05, 0.10816E+05, 0.10888E+05, 0.10959E+05 &
, 0.11031E+05, 0.11103E+05, 0.11176E+05, 0.11248E+05, 0.11321E+05, 0.11394E+05, 0.11467E+05, 0.11541E+05, 0.11614E+05, 0.11688E+05 &
, 0.11762E+05, 0.11837E+05, 0.11911E+05, 0.11986E+05, 0.12061E+05, 0.12136E+05, 0.12211E+05, 0.12287E+05, 0.12363E+05, 0.12439E+05 &
, 0.12515E+05, 0.12592E+05, 0.12669E+05, 0.12746E+05, 0.12823E+05, 0.12900E+05, 0.12978E+05, 0.13056E+05, 0.13134E+05, 0.13212E+05 &
, 0.13291E+05, 0.13370E+05, 0.13449E+05, 0.13528E+05, 0.13608E+05, 0.13688E+05, 0.13768E+05, 0.13848E+05, 0.13928E+05, 0.14009E+05 &
, 0.14090E+05, 0.14171E+05, 0.14253E+05, 0.14334E+05, 0.14416E+05, 0.14498E+05, 0.14581E+05, 0.14663E+05, 0.14746E+05, 0.14829E+05 &
, 0.14913E+05, 0.14996E+05, 0.15080E+05, 0.15164E+05, 0.15248E+05, 0.15333E+05, 0.15418E+05, 0.15503E+05, 0.15588E+05, 0.15674E+05 &
, 0.15759E+05, 0.15845E+05, 0.15932E+05, 0.16018E+05, 0.16105E+05, 0.16192E+05, 0.16279E+05, 0.16367E+05, 0.16455E+05, 0.16543E+05 &
, 0.16631E+05, 0.16720E+05, 0.16809E+05, 0.16898E+05, 0.16987E+05, 0.17077E+05, 0.17166E+05, 0.17256E+05, 0.17347E+05, 0.17437E+05 &
, 0.17528E+05, 0.17620E+05, 0.17711E+05, 0.17803E+05, 0.17895E+05, 0.17987E+05, 0.18079E+05, 0.18172E+05, 0.18265E+05, 0.18358E+05 &
, 0.18452E+05, 0.18545E+05, 0.18639E+05, 0.18734E+05, 0.18828E+05, 0.18923E+05, 0.19018E+05, 0.19114E+05, 0.19209E+05, 0.19305E+05 &
, 0.19402E+05, 0.19498E+05, 0.19595E+05, 0.19692E+05, 0.19789E+05, 0.19887E+05, 0.19985E+05, 0.20083E+05, 0.20181E+05, 0.20280E+05 &
, 0.20379E+05, 0.20478E+05, 0.20578E+05, 0.20678E+05, 0.20778E+05, 0.20878E+05, 0.20979E+05, 0.21080E+05, 0.21181E+05, 0.21283E+05 &
, 0.21384E+05, 0.21486E+05, 0.21589E+05, 0.21692E+05, 0.21795E+05, 0.21898E+05, 0.22001E+05, 0.22105E+05, 0.22209E+05, 0.22314E+05 &
, 0.22418E+05, 0.22523E+05, 0.22629E+05, 0.22734E+05, 0.22840E+05, 0.22946E+05, 0.23053E+05, 0.23160E+05, 0.23267E+05, 0.23374E+05 &
, 0.23482E+05, 0.23590E+05, 0.23698E+05, 0.23807E+05, 0.23916E+05, 0.24025E+05, 0.24134E+05, 0.24244E+05, 0.24354E+05, 0.24465E+05 &
, 0.24576E+05, 0.24687E+05, 0.24798E+05, 0.24910E+05, 0.25022E+05, 0.25134E+05, 0.25247E+05, 0.25359E+05, 0.25473E+05, 0.25586E+05 &
, 0.25700E+05, 0.25814E+05, 0.25929E+05, 0.26044E+05, 0.26159E+05, 0.26274E+05, 0.26390E+05, 0.26506E+05, 0.26622E+05, 0.26739E+05 &
, 0.26856E+05, 0.26974E+05, 0.27091E+05, 0.27209E+05, 0.27328E+05, 0.27446E+05, 0.27565E+05, 0.27685E+05, 0.27804E+05, 0.27924E+05 &
, 0.28045E+05, 0.28165E+05, 0.28286E+05, 0.28408E+05, 0.28529E+05, 0.28651E+05, 0.28774E+05, 0.28896E+05, 0.29019E+05, 0.29143E+05 &
, 0.29266E+05, 0.29390E+05, 0.29515E+05, 0.29639E+05, 0.29764E+05, 0.29890E+05, 0.30015E+05, 0.30141E+05, 0.30268E+05, 0.30395E+05 &
, 0.30522E+05, 0.30649E+05, 0.30777E+05, 0.30905E+05, 0.31034E+05, 0.31162E+05, 0.31292E+05, 0.31421E+05, 0.31551E+05, 0.31681E+05 &
, 0.31812E+05, 0.31943E+05, 0.32074E+05, 0.32206E+05, 0.32338E+05, 0.32470E+05, 0.32603E+05, 0.32736E+05, 0.32870E+05, 0.33003E+05 &
, 0.33138E+05, 0.33272E+05, 0.33407E+05, 0.33542E+05, 0.33678E+05, 0.33814E+05, 0.33950E+05, 0.34087E+05, 0.34224E+05, 0.34362E+05 &
, 0.34500E+05, 0.34638E+05, 0.34776E+05, 0.34915E+05, 0.35055E+05, 0.35194E+05, 0.35335E+05, 0.35475E+05, 0.35616E+05, 0.35757E+05 &
, 0.35899E+05, 0.36041E+05, 0.36183E+05, 0.36326E+05, 0.36469E+05, 0.36612E+05, 0.36756E+05, 0.36901E+05, 0.37045E+05, 0.37190E+05 &
, 0.37336E+05, 0.37482E+05, 0.37628E+05, 0.37775E+05, 0.37922E+05, 0.38069E+05, 0.38217E+05, 0.38365E+05, 0.38513E+05, 0.38662E+05 &
, 0.38812E+05, 0.38962E+05, 0.39112E+05, 0.39262E+05, 0.39413E+05, 0.39565E+05, 0.39716E+05, 0.39869E+05, 0.40021E+05, 0.40174E+05 /
 !            1           7
 DATA(QofT(           7 ,J),J=1,501)/   20.2059993743896       &
, 0.22969E+02, 0.25872E+02, 0.28906E+02, 0.32063E+02, 0.35336E+02, 0.38720E+02, 0.42210E+02, 0.45802E+02, 0.49494E+02, 0.53282E+02 &
, 0.57163E+02, 0.61136E+02, 0.65198E+02, 0.69346E+02, 0.73581E+02, 0.77898E+02, 0.82298E+02, 0.86778E+02, 0.91337E+02, 0.95974E+02 &
, 0.10069E+03, 0.10548E+03, 0.11034E+03, 0.11527E+03, 0.12028E+03, 0.12536E+03, 0.13051E+03, 0.13573E+03, 0.14101E+03, 0.14637E+03 &
, 0.15179E+03, 0.15727E+03, 0.16282E+03, 0.16844E+03, 0.17412E+03, 0.17986E+03, 0.18566E+03, 0.19152E+03, 0.19745E+03, 0.20344E+03 &
, 0.20948E+03, 0.21559E+03, 0.22175E+03, 0.22797E+03, 0.23425E+03, 0.24059E+03, 0.24698E+03, 0.25343E+03, 0.25994E+03, 0.26650E+03 &
, 0.27311E+03, 0.27978E+03, 0.28650E+03, 0.29328E+03, 0.30011E+03, 0.30699E+03, 0.31392E+03, 0.32091E+03, 0.32795E+03, 0.33504E+03 &
, 0.34218E+03, 0.34937E+03, 0.35661E+03, 0.36390E+03, 0.37124E+03, 0.37863E+03, 0.38606E+03, 0.39355E+03, 0.40109E+03, 0.40867E+03 &
, 0.41630E+03, 0.42398E+03, 0.43171E+03, 0.43948E+03, 0.44730E+03, 0.45517E+03, 0.46308E+03, 0.47104E+03, 0.47905E+03, 0.48710E+03 &
, 0.49519E+03, 0.50334E+03, 0.51152E+03, 0.51975E+03, 0.52803E+03, 0.53635E+03, 0.54472E+03, 0.55313E+03, 0.56158E+03, 0.57008E+03 &
, 0.57862E+03, 0.58721E+03, 0.59584E+03, 0.60451E+03, 0.61322E+03, 0.62198E+03, 0.63079E+03, 0.63963E+03, 0.64852E+03, 0.65745E+03 &
, 0.66642E+03, 0.67544E+03, 0.68449E+03, 0.69359E+03, 0.70273E+03, 0.71192E+03, 0.72114E+03, 0.73041E+03, 0.73972E+03, 0.74907E+03 &
, 0.75847E+03, 0.76790E+03, 0.77738E+03, 0.78690E+03, 0.79646E+03, 0.80606E+03, 0.81571E+03, 0.82539E+03, 0.83512E+03, 0.84489E+03 &
, 0.85470E+03, 0.86455E+03, 0.87444E+03, 0.88437E+03, 0.89435E+03, 0.90436E+03, 0.91442E+03, 0.92452E+03, 0.93466E+03, 0.94484E+03 &
, 0.95507E+03, 0.96533E+03, 0.97564E+03, 0.98598E+03, 0.99637E+03, 0.10068E+04, 0.10173E+04, 0.10278E+04, 0.10383E+04, 0.10489E+04 &
, 0.10596E+04, 0.10703E+04, 0.10810E+04, 0.10917E+04, 0.11025E+04, 0.11134E+04, 0.11243E+04, 0.11352E+04, 0.11462E+04, 0.11572E+04 &
, 0.11682E+04, 0.11793E+04, 0.11905E+04, 0.12017E+04, 0.12129E+04, 0.12242E+04, 0.12355E+04, 0.12468E+04, 0.12582E+04, 0.12697E+04 &
, 0.12811E+04, 0.12927E+04, 0.13042E+04, 0.13158E+04, 0.13275E+04, 0.13392E+04, 0.13509E+04, 0.13627E+04, 0.13745E+04, 0.13864E+04 &
, 0.13983E+04, 0.14102E+04, 0.14222E+04, 0.14343E+04, 0.14463E+04, 0.14585E+04, 0.14706E+04, 0.14829E+04, 0.14951E+04, 0.15074E+04 &
, 0.15198E+04, 0.15322E+04, 0.15446E+04, 0.15571E+04, 0.15696E+04, 0.15822E+04, 0.15948E+04, 0.16074E+04, 0.16201E+04, 0.16329E+04 &
, 0.16457E+04, 0.16585E+04, 0.16714E+04, 0.16843E+04, 0.16973E+04, 0.17103E+04, 0.17234E+04, 0.17365E+04, 0.17497E+04, 0.17629E+04 &
, 0.17761E+04, 0.17894E+04, 0.18028E+04, 0.18161E+04, 0.18296E+04, 0.18431E+04, 0.18566E+04, 0.18702E+04, 0.18838E+04, 0.18975E+04 &
, 0.19112E+04, 0.19249E+04, 0.19388E+04, 0.19526E+04, 0.19665E+04, 0.19805E+04, 0.19945E+04, 0.20085E+04, 0.20226E+04, 0.20368E+04 &
, 0.20510E+04, 0.20652E+04, 0.20795E+04, 0.20939E+04, 0.21083E+04, 0.21227E+04, 0.21372E+04, 0.21517E+04, 0.21663E+04, 0.21809E+04 &
, 0.21956E+04, 0.22104E+04, 0.22252E+04, 0.22400E+04, 0.22549E+04, 0.22698E+04, 0.22848E+04, 0.22998E+04, 0.23149E+04, 0.23301E+04 &
, 0.23453E+04, 0.23605E+04, 0.23758E+04, 0.23911E+04, 0.24065E+04, 0.24220E+04, 0.24375E+04, 0.24530E+04, 0.24686E+04, 0.24843E+04 &
, 0.25000E+04, 0.25158E+04, 0.25316E+04, 0.25475E+04, 0.25634E+04, 0.25793E+04, 0.25954E+04, 0.26114E+04, 0.26276E+04, 0.26438E+04 &
, 0.26600E+04, 0.26763E+04, 0.26926E+04, 0.27090E+04, 0.27255E+04, 0.27420E+04, 0.27586E+04, 0.27752E+04, 0.27919E+04, 0.28086E+04 &
, 0.28254E+04, 0.28422E+04, 0.28591E+04, 0.28761E+04, 0.28931E+04, 0.29101E+04, 0.29273E+04, 0.29444E+04, 0.29617E+04, 0.29790E+04 &
, 0.29963E+04, 0.30137E+04, 0.30312E+04, 0.30487E+04, 0.30663E+04, 0.30839E+04, 0.31016E+04, 0.31193E+04, 0.31371E+04, 0.31550E+04 &
, 0.31729E+04, 0.31909E+04, 0.32089E+04, 0.32270E+04, 0.32452E+04, 0.32634E+04, 0.32817E+04, 0.33000E+04, 0.33184E+04, 0.33369E+04 &
, 0.33554E+04, 0.33740E+04, 0.33926E+04, 0.34113E+04, 0.34300E+04, 0.34489E+04, 0.34677E+04, 0.34867E+04, 0.35057E+04, 0.35247E+04 &
, 0.35439E+04, 0.35630E+04, 0.35823E+04, 0.36016E+04, 0.36210E+04, 0.36404E+04, 0.36599E+04, 0.36795E+04, 0.36991E+04, 0.37188E+04 &
, 0.37385E+04, 0.37583E+04, 0.37782E+04, 0.37982E+04, 0.38182E+04, 0.38382E+04, 0.38584E+04, 0.38786E+04, 0.38988E+04, 0.39192E+04 &
, 0.39396E+04, 0.39600E+04, 0.39806E+04, 0.40012E+04, 0.40218E+04, 0.40425E+04, 0.40633E+04, 0.40842E+04, 0.41051E+04, 0.41261E+04 &
, 0.41472E+04, 0.41683E+04, 0.41895E+04, 0.42107E+04, 0.42321E+04, 0.42535E+04, 0.42749E+04, 0.42965E+04, 0.43181E+04, 0.43397E+04 &
, 0.43615E+04, 0.43833E+04, 0.44052E+04, 0.44271E+04, 0.44491E+04, 0.44712E+04, 0.44934E+04, 0.45156E+04, 0.45379E+04, 0.45603E+04 &
, 0.45827E+04, 0.46052E+04, 0.46278E+04, 0.46505E+04, 0.46732E+04, 0.46960E+04, 0.47189E+04, 0.47418E+04, 0.47648E+04, 0.47879E+04 &
, 0.48111E+04, 0.48343E+04, 0.48576E+04, 0.48810E+04, 0.49045E+04, 0.49280E+04, 0.49516E+04, 0.49753E+04, 0.49990E+04, 0.50229E+04 &
, 0.50468E+04, 0.50708E+04, 0.50948E+04, 0.51189E+04, 0.51432E+04, 0.51674E+04, 0.51918E+04, 0.52162E+04, 0.52407E+04, 0.52653E+04 &
, 0.52900E+04, 0.53147E+04, 0.53396E+04, 0.53645E+04, 0.53894E+04, 0.54145E+04, 0.54396E+04, 0.54648E+04, 0.54901E+04, 0.55155E+04 &
, 0.55410E+04, 0.55665E+04, 0.55921E+04, 0.56178E+04, 0.56436E+04, 0.56694E+04, 0.56953E+04, 0.57214E+04, 0.57474E+04, 0.57736E+04 &
, 0.57999E+04, 0.58262E+04, 0.58526E+04, 0.58791E+04, 0.59057E+04, 0.59324E+04, 0.59591E+04, 0.59860E+04, 0.60129E+04, 0.60399E+04 &
, 0.60670E+04, 0.60941E+04, 0.61214E+04, 0.61487E+04, 0.61761E+04, 0.62036E+04, 0.62312E+04, 0.62589E+04, 0.62867E+04, 0.63145E+04 &
, 0.63424E+04, 0.63705E+04, 0.63986E+04, 0.64267E+04, 0.64550E+04, 0.64834E+04, 0.65118E+04, 0.65404E+04, 0.65690E+04, 0.65977E+04 &
, 0.66265E+04, 0.66554E+04, 0.66844E+04, 0.67135E+04, 0.67426E+04, 0.67719E+04, 0.68012E+04, 0.68307E+04, 0.68602E+04, 0.68898E+04 &
, 0.69195E+04, 0.69493E+04, 0.69792E+04, 0.70091E+04, 0.70392E+04, 0.70694E+04, 0.70996E+04, 0.71299E+04, 0.71604E+04, 0.71909E+04 &
, 0.72215E+04, 0.72522E+04, 0.72831E+04, 0.73140E+04, 0.73449E+04, 0.73760E+04, 0.74072E+04, 0.74385E+04, 0.74699E+04, 0.75013E+04 &
, 0.75329E+04, 0.75646E+04, 0.75963E+04, 0.76282E+04, 0.76601E+04, 0.76922E+04, 0.77243E+04, 0.77565E+04, 0.77889E+04, 0.78213E+04 &
, 0.78538E+04, 0.78865E+04, 0.79192E+04, 0.79520E+04, 0.79850E+04, 0.80180E+04, 0.80511E+04, 0.80843E+04, 0.81177E+04, 0.81511E+04 &
, 0.81846E+04, 0.82182E+04, 0.82520E+04, 0.82858E+04, 0.83197E+04, 0.83538E+04, 0.83879E+04, 0.84221E+04, 0.84565E+04, 0.84909E+04 /
 !            1           8
 DATA(QofT(           8 ,J),J=1,501)/   20.4890003204346       &
, 0.23295E+02, 0.26242E+02, 0.29322E+02, 0.32526E+02, 0.35848E+02, 0.39283E+02, 0.42825E+02, 0.46471E+02, 0.50218E+02, 0.54063E+02 &
, 0.58002E+02, 0.62035E+02, 0.66157E+02, 0.70368E+02, 0.74666E+02, 0.79048E+02, 0.83514E+02, 0.88062E+02, 0.92689E+02, 0.97396E+02 &
, 0.10218E+03, 0.10704E+03, 0.11198E+03, 0.11699E+03, 0.12207E+03, 0.12722E+03, 0.13245E+03, 0.13775E+03, 0.14311E+03, 0.14855E+03 &
, 0.15405E+03, 0.15962E+03, 0.16525E+03, 0.17095E+03, 0.17671E+03, 0.18254E+03, 0.18843E+03, 0.19438E+03, 0.20040E+03, 0.20647E+03 &
, 0.21261E+03, 0.21881E+03, 0.22507E+03, 0.23138E+03, 0.23775E+03, 0.24419E+03, 0.25068E+03, 0.25722E+03, 0.26382E+03, 0.27048E+03 &
, 0.27720E+03, 0.28397E+03, 0.29079E+03, 0.29767E+03, 0.30460E+03, 0.31159E+03, 0.31863E+03, 0.32572E+03, 0.33286E+03, 0.34006E+03 &
, 0.34730E+03, 0.35460E+03, 0.36195E+03, 0.36935E+03, 0.37680E+03, 0.38431E+03, 0.39186E+03, 0.39946E+03, 0.40710E+03, 0.41480E+03 &
, 0.42255E+03, 0.43034E+03, 0.43819E+03, 0.44608E+03, 0.45402E+03, 0.46200E+03, 0.47003E+03, 0.47811E+03, 0.48624E+03, 0.49441E+03 &
, 0.50263E+03, 0.51090E+03, 0.51921E+03, 0.52757E+03, 0.53597E+03, 0.54441E+03, 0.55291E+03, 0.56144E+03, 0.57003E+03, 0.57865E+03 &
, 0.58732E+03, 0.59604E+03, 0.60480E+03, 0.61360E+03, 0.62245E+03, 0.63134E+03, 0.64028E+03, 0.64926E+03, 0.65828E+03, 0.66735E+03 &
, 0.67646E+03, 0.68561E+03, 0.69481E+03, 0.70404E+03, 0.71333E+03, 0.72265E+03, 0.73202E+03, 0.74143E+03, 0.75088E+03, 0.76037E+03 &
, 0.76991E+03, 0.77949E+03, 0.78911E+03, 0.79878E+03, 0.80849E+03, 0.81824E+03, 0.82803E+03, 0.83786E+03, 0.84774E+03, 0.85766E+03 &
, 0.86762E+03, 0.87762E+03, 0.88767E+03, 0.89775E+03, 0.90788E+03, 0.91805E+03, 0.92827E+03, 0.93852E+03, 0.94882E+03, 0.95916E+03 &
, 0.96954E+03, 0.97997E+03, 0.99043E+03, 0.10009E+04, 0.10115E+04, 0.10221E+04, 0.10327E+04, 0.10434E+04, 0.10541E+04, 0.10649E+04 &
, 0.10757E+04, 0.10865E+04, 0.10974E+04, 0.11084E+04, 0.11193E+04, 0.11303E+04, 0.11414E+04, 0.11525E+04, 0.11637E+04, 0.11748E+04 &
, 0.11861E+04, 0.11973E+04, 0.12087E+04, 0.12200E+04, 0.12314E+04, 0.12429E+04, 0.12544E+04, 0.12659E+04, 0.12775E+04, 0.12891E+04 &
, 0.13007E+04, 0.13125E+04, 0.13242E+04, 0.13360E+04, 0.13478E+04, 0.13597E+04, 0.13716E+04, 0.13836E+04, 0.13956E+04, 0.14077E+04 &
, 0.14198E+04, 0.14319E+04, 0.14441E+04, 0.14563E+04, 0.14686E+04, 0.14810E+04, 0.14933E+04, 0.15057E+04, 0.15182E+04, 0.15307E+04 &
, 0.15432E+04, 0.15558E+04, 0.15685E+04, 0.15812E+04, 0.15939E+04, 0.16067E+04, 0.16195E+04, 0.16324E+04, 0.16453E+04, 0.16582E+04 &
, 0.16712E+04, 0.16843E+04, 0.16974E+04, 0.17105E+04, 0.17237E+04, 0.17369E+04, 0.17502E+04, 0.17636E+04, 0.17769E+04, 0.17904E+04 &
, 0.18038E+04, 0.18173E+04, 0.18309E+04, 0.18445E+04, 0.18582E+04, 0.18719E+04, 0.18857E+04, 0.18995E+04, 0.19133E+04, 0.19272E+04 &
, 0.19412E+04, 0.19552E+04, 0.19692E+04, 0.19833E+04, 0.19974E+04, 0.20116E+04, 0.20259E+04, 0.20402E+04, 0.20545E+04, 0.20689E+04 &
, 0.20833E+04, 0.20978E+04, 0.21124E+04, 0.21269E+04, 0.21416E+04, 0.21563E+04, 0.21710E+04, 0.21858E+04, 0.22006E+04, 0.22155E+04 &
, 0.22305E+04, 0.22455E+04, 0.22605E+04, 0.22756E+04, 0.22907E+04, 0.23059E+04, 0.23212E+04, 0.23365E+04, 0.23518E+04, 0.23672E+04 &
, 0.23827E+04, 0.23982E+04, 0.24138E+04, 0.24294E+04, 0.24450E+04, 0.24608E+04, 0.24765E+04, 0.24924E+04, 0.25082E+04, 0.25242E+04 &
, 0.25401E+04, 0.25562E+04, 0.25723E+04, 0.25884E+04, 0.26046E+04, 0.26209E+04, 0.26372E+04, 0.26536E+04, 0.26700E+04, 0.26864E+04 &
, 0.27030E+04, 0.27196E+04, 0.27362E+04, 0.27529E+04, 0.27696E+04, 0.27864E+04, 0.28033E+04, 0.28202E+04, 0.28372E+04, 0.28542E+04 &
, 0.28713E+04, 0.28884E+04, 0.29056E+04, 0.29229E+04, 0.29402E+04, 0.29576E+04, 0.29750E+04, 0.29925E+04, 0.30100E+04, 0.30276E+04 &
, 0.30453E+04, 0.30630E+04, 0.30808E+04, 0.30986E+04, 0.31165E+04, 0.31345E+04, 0.31525E+04, 0.31706E+04, 0.31887E+04, 0.32069E+04 &
, 0.32251E+04, 0.32434E+04, 0.32618E+04, 0.32802E+04, 0.32987E+04, 0.33173E+04, 0.33359E+04, 0.33546E+04, 0.33733E+04, 0.33921E+04 &
, 0.34110E+04, 0.34299E+04, 0.34488E+04, 0.34679E+04, 0.34870E+04, 0.35062E+04, 0.35254E+04, 0.35447E+04, 0.35640E+04, 0.35834E+04 &
, 0.36029E+04, 0.36225E+04, 0.36421E+04, 0.36618E+04, 0.36815E+04, 0.37013E+04, 0.37211E+04, 0.37411E+04, 0.37611E+04, 0.37811E+04 &
, 0.38012E+04, 0.38214E+04, 0.38417E+04, 0.38620E+04, 0.38824E+04, 0.39028E+04, 0.39234E+04, 0.39439E+04, 0.39646E+04, 0.39853E+04 &
, 0.40061E+04, 0.40269E+04, 0.40478E+04, 0.40688E+04, 0.40899E+04, 0.41110E+04, 0.41322E+04, 0.41534E+04, 0.41748E+04, 0.41962E+04 &
, 0.42176E+04, 0.42391E+04, 0.42607E+04, 0.42824E+04, 0.43042E+04, 0.43260E+04, 0.43478E+04, 0.43698E+04, 0.43918E+04, 0.44139E+04 &
, 0.44361E+04, 0.44583E+04, 0.44806E+04, 0.45030E+04, 0.45254E+04, 0.45479E+04, 0.45705E+04, 0.45932E+04, 0.46159E+04, 0.46387E+04 &
, 0.46616E+04, 0.46845E+04, 0.47076E+04, 0.47307E+04, 0.47538E+04, 0.47771E+04, 0.48004E+04, 0.48238E+04, 0.48472E+04, 0.48708E+04 &
, 0.48944E+04, 0.49181E+04, 0.49419E+04, 0.49657E+04, 0.49896E+04, 0.50136E+04, 0.50377E+04, 0.50618E+04, 0.50860E+04, 0.51103E+04 &
, 0.51347E+04, 0.51592E+04, 0.51837E+04, 0.52083E+04, 0.52330E+04, 0.52578E+04, 0.52826E+04, 0.53075E+04, 0.53325E+04, 0.53576E+04 &
, 0.53828E+04, 0.54080E+04, 0.54333E+04, 0.54587E+04, 0.54842E+04, 0.55097E+04, 0.55354E+04, 0.55611E+04, 0.55869E+04, 0.56128E+04 &
, 0.56387E+04, 0.56648E+04, 0.56909E+04, 0.57171E+04, 0.57434E+04, 0.57697E+04, 0.57962E+04, 0.58227E+04, 0.58493E+04, 0.58760E+04 &
, 0.59028E+04, 0.59297E+04, 0.59566E+04, 0.59837E+04, 0.60108E+04, 0.60380E+04, 0.60653E+04, 0.60927E+04, 0.61201E+04, 0.61477E+04 &
, 0.61753E+04, 0.62030E+04, 0.62308E+04, 0.62587E+04, 0.62867E+04, 0.63147E+04, 0.63429E+04, 0.63711E+04, 0.63994E+04, 0.64279E+04 &
, 0.64564E+04, 0.64850E+04, 0.65136E+04, 0.65424E+04, 0.65713E+04, 0.66002E+04, 0.66292E+04, 0.66584E+04, 0.66876E+04, 0.67169E+04 &
, 0.67463E+04, 0.67758E+04, 0.68053E+04, 0.68350E+04, 0.68648E+04, 0.68946E+04, 0.69246E+04, 0.69546E+04, 0.69847E+04, 0.70149E+04 &
, 0.70452E+04, 0.70757E+04, 0.71062E+04, 0.71367E+04, 0.71674E+04, 0.71982E+04, 0.72291E+04, 0.72601E+04, 0.72911E+04, 0.73223E+04 &
, 0.73535E+04, 0.73849E+04, 0.74163E+04, 0.74479E+04, 0.74795E+04, 0.75112E+04, 0.75431E+04, 0.75750E+04, 0.76070E+04, 0.76392E+04 &
, 0.76714E+04, 0.77037E+04, 0.77361E+04, 0.77686E+04, 0.78012E+04, 0.78340E+04, 0.78668E+04, 0.78997E+04, 0.79327E+04, 0.79658E+04 &
, 0.79990E+04, 0.80323E+04, 0.80658E+04, 0.80993E+04, 0.81329E+04, 0.81666E+04, 0.82004E+04, 0.82344E+04, 0.82684E+04, 0.83025E+04 &
, 0.83367E+04, 0.83711E+04, 0.84055E+04, 0.84401E+04, 0.84747E+04, 0.85095E+04, 0.85443E+04, 0.85793E+04, 0.86143E+04, 0.86495E+04 /
 !            1           9
 DATA(QofT(           9 ,J),J=1,501)/   122.129997253418       &
, 0.13884E+03, 0.15640E+03, 0.17475E+03, 0.19384E+03, 0.21364E+03, 0.23410E+03, 0.25520E+03, 0.27693E+03, 0.29925E+03, 0.32216E+03 &
, 0.34563E+03, 0.36966E+03, 0.39422E+03, 0.41931E+03, 0.44492E+03, 0.47103E+03, 0.49763E+03, 0.52473E+03, 0.55230E+03, 0.58034E+03 &
, 0.60884E+03, 0.63780E+03, 0.66721E+03, 0.69706E+03, 0.72734E+03, 0.75805E+03, 0.78919E+03, 0.82074E+03, 0.85271E+03, 0.88508E+03 &
, 0.91786E+03, 0.95103E+03, 0.98460E+03, 0.10186E+04, 0.10529E+04, 0.10876E+04, 0.11227E+04, 0.11582E+04, 0.11940E+04, 0.12302E+04 &
, 0.12668E+04, 0.13037E+04, 0.13410E+04, 0.13786E+04, 0.14166E+04, 0.14549E+04, 0.14936E+04, 0.15326E+04, 0.15719E+04, 0.16116E+04 &
, 0.16516E+04, 0.16919E+04, 0.17325E+04, 0.17735E+04, 0.18148E+04, 0.18564E+04, 0.18984E+04, 0.19406E+04, 0.19832E+04, 0.20261E+04 &
, 0.20692E+04, 0.21127E+04, 0.21565E+04, 0.22006E+04, 0.22450E+04, 0.22897E+04, 0.23347E+04, 0.23799E+04, 0.24255E+04, 0.24714E+04 &
, 0.25175E+04, 0.25640E+04, 0.26107E+04, 0.26577E+04, 0.27050E+04, 0.27526E+04, 0.28004E+04, 0.28486E+04, 0.28970E+04, 0.29457E+04 &
, 0.29946E+04, 0.30439E+04, 0.30934E+04, 0.31432E+04, 0.31932E+04, 0.32436E+04, 0.32942E+04, 0.33450E+04, 0.33961E+04, 0.34475E+04 &
, 0.34992E+04, 0.35511E+04, 0.36033E+04, 0.36558E+04, 0.37085E+04, 0.37614E+04, 0.38147E+04, 0.38682E+04, 0.39219E+04, 0.39759E+04 &
, 0.40302E+04, 0.40847E+04, 0.41395E+04, 0.41945E+04, 0.42498E+04, 0.43054E+04, 0.43612E+04, 0.44172E+04, 0.44735E+04, 0.45301E+04 &
, 0.45869E+04, 0.46440E+04, 0.47013E+04, 0.47589E+04, 0.48167E+04, 0.48748E+04, 0.49331E+04, 0.49917E+04, 0.50505E+04, 0.51096E+04 &
, 0.51690E+04, 0.52285E+04, 0.52884E+04, 0.53485E+04, 0.54088E+04, 0.54694E+04, 0.55302E+04, 0.55913E+04, 0.56527E+04, 0.57143E+04 &
, 0.57761E+04, 0.58382E+04, 0.59005E+04, 0.59631E+04, 0.60260E+04, 0.60891E+04, 0.61524E+04, 0.62160E+04, 0.62799E+04, 0.63440E+04 &
, 0.64083E+04, 0.64729E+04, 0.65378E+04, 0.66029E+04, 0.66683E+04, 0.67339E+04, 0.67997E+04, 0.68659E+04, 0.69323E+04, 0.69989E+04 &
, 0.70658E+04, 0.71329E+04, 0.72003E+04, 0.72680E+04, 0.73359E+04, 0.74040E+04, 0.74725E+04, 0.75411E+04, 0.76101E+04, 0.76793E+04 &
, 0.77487E+04, 0.78184E+04, 0.78884E+04, 0.79587E+04, 0.80291E+04, 0.80999E+04, 0.81709E+04, 0.82422E+04, 0.83137E+04, 0.83855E+04 &
, 0.84576E+04, 0.85299E+04, 0.86025E+04, 0.86753E+04, 0.87485E+04, 0.88218E+04, 0.88955E+04, 0.89694E+04, 0.90436E+04, 0.91180E+04 &
, 0.91928E+04, 0.92677E+04, 0.93430E+04, 0.94185E+04, 0.94943E+04, 0.95704E+04, 0.96467E+04, 0.97233E+04, 0.98002E+04, 0.98774E+04 &
, 0.99548E+04, 0.10032E+05, 0.10110E+05, 0.10189E+05, 0.10267E+05, 0.10346E+05, 0.10425E+05, 0.10505E+05, 0.10584E+05, 0.10664E+05 &
, 0.10744E+05, 0.10825E+05, 0.10906E+05, 0.10987E+05, 0.11068E+05, 0.11149E+05, 0.11231E+05, 0.11314E+05, 0.11396E+05, 0.11479E+05 &
, 0.11562E+05, 0.11645E+05, 0.11729E+05, 0.11813E+05, 0.11897E+05, 0.11981E+05, 0.12066E+05, 0.12151E+05, 0.12237E+05, 0.12322E+05 &
, 0.12408E+05, 0.12494E+05, 0.12581E+05, 0.12668E+05, 0.12755E+05, 0.12842E+05, 0.12930E+05, 0.13018E+05, 0.13106E+05, 0.13195E+05 &
, 0.13284E+05, 0.13373E+05, 0.13463E+05, 0.13552E+05, 0.13643E+05, 0.13733E+05, 0.13824E+05, 0.13915E+05, 0.14006E+05, 0.14098E+05 &
, 0.14190E+05, 0.14282E+05, 0.14375E+05, 0.14468E+05, 0.14561E+05, 0.14654E+05, 0.14748E+05, 0.14843E+05, 0.14937E+05, 0.15032E+05 &
, 0.15127E+05, 0.15222E+05, 0.15318E+05, 0.15414E+05, 0.15511E+05, 0.15607E+05, 0.15704E+05, 0.15802E+05, 0.15899E+05, 0.15997E+05 &
, 0.16096E+05, 0.16194E+05, 0.16293E+05, 0.16393E+05, 0.16492E+05, 0.16592E+05, 0.16693E+05, 0.16793E+05, 0.16894E+05, 0.16996E+05 &
, 0.17097E+05, 0.17199E+05, 0.17302E+05, 0.17404E+05, 0.17507E+05, 0.17611E+05, 0.17714E+05, 0.17819E+05, 0.17923E+05, 0.18028E+05 &
, 0.18133E+05, 0.18238E+05, 0.18344E+05, 0.18450E+05, 0.18556E+05, 0.18663E+05, 0.18770E+05, 0.18878E+05, 0.18986E+05, 0.19094E+05 &
, 0.19202E+05, 0.19311E+05, 0.19421E+05, 0.19530E+05, 0.19640E+05, 0.19751E+05, 0.19861E+05, 0.19972E+05, 0.20084E+05, 0.20196E+05 &
, 0.20308E+05, 0.20420E+05, 0.20533E+05, 0.20647E+05, 0.20760E+05, 0.20874E+05, 0.20989E+05, 0.21103E+05, 0.21218E+05, 0.21334E+05 &
, 0.21450E+05, 0.21566E+05, 0.21683E+05, 0.21800E+05, 0.21917E+05, 0.22035E+05, 0.22153E+05, 0.22271E+05, 0.22390E+05, 0.22510E+05 &
, 0.22629E+05, 0.22749E+05, 0.22870E+05, 0.22991E+05, 0.23112E+05, 0.23233E+05, 0.23355E+05, 0.23478E+05, 0.23601E+05, 0.23724E+05 &
, 0.23847E+05, 0.23971E+05, 0.24096E+05, 0.24221E+05, 0.24346E+05, 0.24471E+05, 0.24597E+05, 0.24724E+05, 0.24851E+05, 0.24978E+05 &
, 0.25105E+05, 0.25233E+05, 0.25362E+05, 0.25491E+05, 0.25620E+05, 0.25750E+05, 0.25880E+05, 0.26010E+05, 0.26141E+05, 0.26272E+05 &
, 0.26404E+05, 0.26536E+05, 0.26669E+05, 0.26802E+05, 0.26936E+05, 0.27069E+05, 0.27204E+05, 0.27338E+05, 0.27474E+05, 0.27609E+05 &
, 0.27745E+05, 0.27882E+05, 0.28018E+05, 0.28156E+05, 0.28294E+05, 0.28432E+05, 0.28570E+05, 0.28709E+05, 0.28849E+05, 0.28989E+05 &
, 0.29129E+05, 0.29270E+05, 0.29411E+05, 0.29553E+05, 0.29695E+05, 0.29838E+05, 0.29981E+05, 0.30125E+05, 0.30269E+05, 0.30413E+05 &
, 0.30558E+05, 0.30703E+05, 0.30849E+05, 0.30995E+05, 0.31142E+05, 0.31289E+05, 0.31437E+05, 0.31585E+05, 0.31734E+05, 0.31883E+05 &
, 0.32032E+05, 0.32182E+05, 0.32333E+05, 0.32484E+05, 0.32635E+05, 0.32787E+05, 0.32939E+05, 0.33092E+05, 0.33246E+05, 0.33399E+05 &
, 0.33554E+05, 0.33709E+05, 0.33864E+05, 0.34020E+05, 0.34176E+05, 0.34332E+05, 0.34490E+05, 0.34647E+05, 0.34806E+05, 0.34964E+05 &
, 0.35123E+05, 0.35283E+05, 0.35443E+05, 0.35604E+05, 0.35765E+05, 0.35927E+05, 0.36089E+05, 0.36252E+05, 0.36415E+05, 0.36579E+05 &
, 0.36743E+05, 0.36907E+05, 0.37073E+05, 0.37238E+05, 0.37405E+05, 0.37571E+05, 0.37739E+05, 0.37906E+05, 0.38075E+05, 0.38244E+05 &
, 0.38413E+05, 0.38583E+05, 0.38753E+05, 0.38924E+05, 0.39096E+05, 0.39268E+05, 0.39440E+05, 0.39613E+05, 0.39787E+05, 0.39961E+05 &
, 0.40136E+05, 0.40311E+05, 0.40487E+05, 0.40663E+05, 0.40840E+05, 0.41017E+05, 0.41195E+05, 0.41373E+05, 0.41552E+05, 0.41732E+05 &
, 0.41912E+05, 0.42093E+05, 0.42274E+05, 0.42456E+05, 0.42638E+05, 0.42821E+05, 0.43004E+05, 0.43188E+05, 0.43373E+05, 0.43558E+05 &
, 0.43744E+05, 0.43930E+05, 0.44117E+05, 0.44304E+05, 0.44492E+05, 0.44681E+05, 0.44870E+05, 0.45060E+05, 0.45250E+05, 0.45441E+05 &
, 0.45632E+05, 0.45824E+05, 0.46017E+05, 0.46210E+05, 0.46404E+05, 0.46598E+05, 0.46793E+05, 0.46988E+05, 0.47184E+05, 0.47381E+05 &
, 0.47578E+05, 0.47776E+05, 0.47975E+05, 0.48174E+05, 0.48374E+05, 0.48574E+05, 0.48775E+05, 0.48976E+05, 0.49178E+05, 0.49381E+05 &
, 0.49584E+05, 0.49788E+05, 0.49993E+05, 0.50198E+05, 0.50404E+05, 0.50610E+05, 0.50817E+05, 0.51025E+05, 0.51233E+05, 0.51442E+05 /
 !            2           1
 DATA(QofT(          10 ,J),J=1,501)/   17.9790000915527       &
, 0.19760E+02, 0.21542E+02, 0.23323E+02, 0.25104E+02, 0.26885E+02, 0.28666E+02, 0.30448E+02, 0.32229E+02, 0.34010E+02, 0.35792E+02 &
, 0.37573E+02, 0.39354E+02, 0.41135E+02, 0.42917E+02, 0.44698E+02, 0.46479E+02, 0.48261E+02, 0.50042E+02, 0.51823E+02, 0.53605E+02 &
, 0.55386E+02, 0.57168E+02, 0.58949E+02, 0.60731E+02, 0.62512E+02, 0.64294E+02, 0.66075E+02, 0.67857E+02, 0.69638E+02, 0.71420E+02 &
, 0.73202E+02, 0.74984E+02, 0.76766E+02, 0.78548E+02, 0.80330E+02, 0.82113E+02, 0.83896E+02, 0.85679E+02, 0.87462E+02, 0.89246E+02 &
, 0.91031E+02, 0.92815E+02, 0.94601E+02, 0.96387E+02, 0.98174E+02, 0.99961E+02, 0.10175E+03, 0.10354E+03, 0.10533E+03, 0.10712E+03 &
, 0.10891E+03, 0.11071E+03, 0.11250E+03, 0.11430E+03, 0.11610E+03, 0.11790E+03, 0.11971E+03, 0.12151E+03, 0.12332E+03, 0.12513E+03 &
, 0.12694E+03, 0.12876E+03, 0.13057E+03, 0.13239E+03, 0.13422E+03, 0.13605E+03, 0.13788E+03, 0.13971E+03, 0.14155E+03, 0.14339E+03 &
, 0.14524E+03, 0.14709E+03, 0.14894E+03, 0.15080E+03, 0.15266E+03, 0.15453E+03, 0.15641E+03, 0.15829E+03, 0.16017E+03, 0.16206E+03 &
, 0.16395E+03, 0.16586E+03, 0.16776E+03, 0.16968E+03, 0.17159E+03, 0.17352E+03, 0.17545E+03, 0.17739E+03, 0.17934E+03, 0.18129E+03 &
, 0.18325E+03, 0.18522E+03, 0.18719E+03, 0.18918E+03, 0.19117E+03, 0.19317E+03, 0.19517E+03, 0.19719E+03, 0.19921E+03, 0.20124E+03 &
, 0.20328E+03, 0.20533E+03, 0.20739E+03, 0.20946E+03, 0.21153E+03, 0.21362E+03, 0.21572E+03, 0.21782E+03, 0.21994E+03, 0.22206E+03 &
, 0.22419E+03, 0.22634E+03, 0.22850E+03, 0.23066E+03, 0.23284E+03, 0.23502E+03, 0.23722E+03, 0.23943E+03, 0.24165E+03, 0.24388E+03 &
, 0.24612E+03, 0.24838E+03, 0.25064E+03, 0.25292E+03, 0.25521E+03, 0.25751E+03, 0.25982E+03, 0.26215E+03, 0.26449E+03, 0.26684E+03 &
, 0.26920E+03, 0.27157E+03, 0.27396E+03, 0.27636E+03, 0.27877E+03, 0.28120E+03, 0.28364E+03, 0.28609E+03, 0.28856E+03, 0.29104E+03 &
, 0.29353E+03, 0.29604E+03, 0.29856E+03, 0.30110E+03, 0.30365E+03, 0.30621E+03, 0.30879E+03, 0.31138E+03, 0.31399E+03, 0.31661E+03 &
, 0.31925E+03, 0.32190E+03, 0.32456E+03, 0.32724E+03, 0.32994E+03, 0.33265E+03, 0.33538E+03, 0.33812E+03, 0.34088E+03, 0.34365E+03 &
, 0.34644E+03, 0.34925E+03, 0.35207E+03, 0.35491E+03, 0.35776E+03, 0.36063E+03, 0.36352E+03, 0.36642E+03, 0.36934E+03, 0.37228E+03 &
, 0.37523E+03, 0.37820E+03, 0.38119E+03, 0.38419E+03, 0.38721E+03, 0.39025E+03, 0.39331E+03, 0.39638E+03, 0.39948E+03, 0.40259E+03 &
, 0.40571E+03, 0.40886E+03, 0.41202E+03, 0.41520E+03, 0.41840E+03, 0.42162E+03, 0.42486E+03, 0.42811E+03, 0.43139E+03, 0.43468E+03 &
, 0.43799E+03, 0.44132E+03, 0.44467E+03, 0.44804E+03, 0.45143E+03, 0.45484E+03, 0.45827E+03, 0.46171E+03, 0.46518E+03, 0.46867E+03 &
, 0.47218E+03, 0.47570E+03, 0.47925E+03, 0.48282E+03, 0.48640E+03, 0.49001E+03, 0.49364E+03, 0.49729E+03, 0.50096E+03, 0.50465E+03 &
, 0.50836E+03, 0.51210E+03, 0.51585E+03, 0.51963E+03, 0.52343E+03, 0.52725E+03, 0.53109E+03, 0.53495E+03, 0.53883E+03, 0.54274E+03 &
, 0.54667E+03, 0.55062E+03, 0.55459E+03, 0.55859E+03, 0.56261E+03, 0.56665E+03, 0.57071E+03, 0.57480E+03, 0.57891E+03, 0.58304E+03 &
, 0.58719E+03, 0.59137E+03, 0.59557E+03, 0.59980E+03, 0.60405E+03, 0.60832E+03, 0.61262E+03, 0.61694E+03, 0.62129E+03, 0.62566E+03 &
, 0.63005E+03, 0.63447E+03, 0.63891E+03, 0.64338E+03, 0.64787E+03, 0.65239E+03, 0.65693E+03, 0.66150E+03, 0.66609E+03, 0.67071E+03 &
, 0.67536E+03, 0.68003E+03, 0.68472E+03, 0.68944E+03, 0.69419E+03, 0.69896E+03, 0.70376E+03, 0.70859E+03, 0.71344E+03, 0.71832E+03 &
, 0.72322E+03, 0.72815E+03, 0.73311E+03, 0.73810E+03, 0.74311E+03, 0.74815E+03, 0.75322E+03, 0.75832E+03, 0.76344E+03, 0.76859E+03 &
, 0.77377E+03, 0.77898E+03, 0.78421E+03, 0.78948E+03, 0.79477E+03, 0.80009E+03, 0.80544E+03, 0.81082E+03, 0.81622E+03, 0.82166E+03 &
, 0.82712E+03, 0.83262E+03, 0.83814E+03, 0.84370E+03, 0.84928E+03, 0.85489E+03, 0.86054E+03, 0.86621E+03, 0.87191E+03, 0.87765E+03 &
, 0.88341E+03, 0.88920E+03, 0.89503E+03, 0.90088E+03, 0.90677E+03, 0.91269E+03, 0.91864E+03, 0.92462E+03, 0.93063E+03, 0.93668E+03 &
, 0.94275E+03, 0.94886E+03, 0.95500E+03, 0.96117E+03, 0.96738E+03, 0.97361E+03, 0.97988E+03, 0.98619E+03, 0.99252E+03, 0.99889E+03 &
, 0.10053E+04, 0.10117E+04, 0.10182E+04, 0.10247E+04, 0.10312E+04, 0.10378E+04, 0.10444E+04, 0.10510E+04, 0.10577E+04, 0.10644E+04 &
, 0.10712E+04, 0.10779E+04, 0.10847E+04, 0.10916E+04, 0.10985E+04, 0.11054E+04, 0.11123E+04, 0.11193E+04, 0.11263E+04, 0.11334E+04 &
, 0.11405E+04, 0.11476E+04, 0.11548E+04, 0.11620E+04, 0.11692E+04, 0.11765E+04, 0.11838E+04, 0.11912E+04, 0.11986E+04, 0.12060E+04 &
, 0.12135E+04, 0.12210E+04, 0.12285E+04, 0.12361E+04, 0.12437E+04, 0.12513E+04, 0.12590E+04, 0.12668E+04, 0.12745E+04, 0.12823E+04 &
, 0.12902E+04, 0.12981E+04, 0.13060E+04, 0.13140E+04, 0.13220E+04, 0.13300E+04, 0.13381E+04, 0.13462E+04, 0.13544E+04, 0.13626E+04 &
, 0.13708E+04, 0.13791E+04, 0.13874E+04, 0.13958E+04, 0.14042E+04, 0.14127E+04, 0.14212E+04, 0.14297E+04, 0.14383E+04, 0.14469E+04 &
, 0.14555E+04, 0.14642E+04, 0.14730E+04, 0.14818E+04, 0.14906E+04, 0.14995E+04, 0.15084E+04, 0.15174E+04, 0.15264E+04, 0.15354E+04 &
, 0.15445E+04, 0.15536E+04, 0.15628E+04, 0.15720E+04, 0.15813E+04, 0.15906E+04, 0.16000E+04, 0.16094E+04, 0.16188E+04, 0.16283E+04 &
, 0.16379E+04, 0.16474E+04, 0.16571E+04, 0.16668E+04, 0.16765E+04, 0.16862E+04, 0.16961E+04, 0.17059E+04, 0.17158E+04, 0.17258E+04 &
, 0.17358E+04, 0.17458E+04, 0.17559E+04, 0.17661E+04, 0.17763E+04, 0.17865E+04, 0.17968E+04, 0.18072E+04, 0.18176E+04, 0.18280E+04 &
, 0.18385E+04, 0.18490E+04, 0.18596E+04, 0.18702E+04, 0.18809E+04, 0.18917E+04, 0.19024E+04, 0.19133E+04, 0.19242E+04, 0.19351E+04 &
, 0.19461E+04, 0.19571E+04, 0.19682E+04, 0.19794E+04, 0.19906E+04, 0.20018E+04, 0.20131E+04, 0.20245E+04, 0.20359E+04, 0.20473E+04 &
, 0.20588E+04, 0.20704E+04, 0.20820E+04, 0.20937E+04, 0.21054E+04, 0.21172E+04, 0.21290E+04, 0.21409E+04, 0.21528E+04, 0.21648E+04 &
, 0.21769E+04, 0.21890E+04, 0.22011E+04, 0.22133E+04, 0.22256E+04, 0.22379E+04, 0.22503E+04, 0.22627E+04, 0.22752E+04, 0.22878E+04 &
, 0.23004E+04, 0.23131E+04, 0.23258E+04, 0.23386E+04, 0.23514E+04, 0.23643E+04, 0.23772E+04, 0.23902E+04, 0.24033E+04, 0.24164E+04 &
, 0.24296E+04, 0.24429E+04, 0.24562E+04, 0.24695E+04, 0.24829E+04, 0.24964E+04, 0.25100E+04, 0.25236E+04, 0.25372E+04, 0.25510E+04 &
, 0.25647E+04, 0.25786E+04, 0.25925E+04, 0.26065E+04, 0.26205E+04, 0.26346E+04, 0.26487E+04, 0.26629E+04, 0.26772E+04, 0.26916E+04 &
, 0.27060E+04, 0.27204E+04, 0.27350E+04, 0.27496E+04, 0.27642E+04, 0.27789E+04, 0.27937E+04, 0.28086E+04, 0.28235E+04, 0.28385E+04 &
, 0.28535E+04, 0.28686E+04, 0.28838E+04, 0.28990E+04, 0.29144E+04, 0.29297E+04, 0.29452E+04, 0.29607E+04, 0.29763E+04, 0.29919E+04 /
 !            2           2
 DATA(QofT(          11 ,J),J=1,501)/   35.9570007324219       &
, 0.39519E+02, 0.43081E+02, 0.46644E+02, 0.50206E+02, 0.53768E+02, 0.57330E+02, 0.60893E+02, 0.64455E+02, 0.68017E+02, 0.71580E+02 &
, 0.75142E+02, 0.78705E+02, 0.82267E+02, 0.85830E+02, 0.89392E+02, 0.92955E+02, 0.96517E+02, 0.10008E+03, 0.10364E+03, 0.10720E+03 &
, 0.11077E+03, 0.11433E+03, 0.11789E+03, 0.12146E+03, 0.12502E+03, 0.12858E+03, 0.13214E+03, 0.13571E+03, 0.13927E+03, 0.14283E+03 &
, 0.14640E+03, 0.14996E+03, 0.15353E+03, 0.15709E+03, 0.16066E+03, 0.16422E+03, 0.16779E+03, 0.17135E+03, 0.17492E+03, 0.17849E+03 &
, 0.18206E+03, 0.18563E+03, 0.18921E+03, 0.19278E+03, 0.19636E+03, 0.19993E+03, 0.20351E+03, 0.20710E+03, 0.21068E+03, 0.21427E+03 &
, 0.21786E+03, 0.22145E+03, 0.22505E+03, 0.22865E+03, 0.23226E+03, 0.23587E+03, 0.23948E+03, 0.24310E+03, 0.24673E+03, 0.25036E+03 &
, 0.25399E+03, 0.25764E+03, 0.26128E+03, 0.26494E+03, 0.26860E+03, 0.27227E+03, 0.27595E+03, 0.27964E+03, 0.28333E+03, 0.28703E+03 &
, 0.29074E+03, 0.29446E+03, 0.29820E+03, 0.30194E+03, 0.30569E+03, 0.30945E+03, 0.31322E+03, 0.31701E+03, 0.32081E+03, 0.32461E+03 &
, 0.32844E+03, 0.33227E+03, 0.33612E+03, 0.33998E+03, 0.34385E+03, 0.34774E+03, 0.35164E+03, 0.35556E+03, 0.35949E+03, 0.36344E+03 &
, 0.36740E+03, 0.37138E+03, 0.37538E+03, 0.37939E+03, 0.38342E+03, 0.38746E+03, 0.39153E+03, 0.39561E+03, 0.39971E+03, 0.40383E+03 &
, 0.40797E+03, 0.41212E+03, 0.41630E+03, 0.42049E+03, 0.42471E+03, 0.42894E+03, 0.43320E+03, 0.43747E+03, 0.44177E+03, 0.44609E+03 &
, 0.45042E+03, 0.45479E+03, 0.45917E+03, 0.46357E+03, 0.46800E+03, 0.47245E+03, 0.47693E+03, 0.48142E+03, 0.48594E+03, 0.49049E+03 &
, 0.49506E+03, 0.49965E+03, 0.50427E+03, 0.50891E+03, 0.51358E+03, 0.51827E+03, 0.52299E+03, 0.52773E+03, 0.53250E+03, 0.53730E+03 &
, 0.54212E+03, 0.54697E+03, 0.55185E+03, 0.55675E+03, 0.56168E+03, 0.56664E+03, 0.57163E+03, 0.57664E+03, 0.58169E+03, 0.58676E+03 &
, 0.59186E+03, 0.59699E+03, 0.60215E+03, 0.60734E+03, 0.61255E+03, 0.61780E+03, 0.62308E+03, 0.62839E+03, 0.63373E+03, 0.63910E+03 &
, 0.64450E+03, 0.64993E+03, 0.65540E+03, 0.66089E+03, 0.66642E+03, 0.67198E+03, 0.67757E+03, 0.68319E+03, 0.68885E+03, 0.69454E+03 &
, 0.70026E+03, 0.70602E+03, 0.71181E+03, 0.71764E+03, 0.72349E+03, 0.72939E+03, 0.73531E+03, 0.74127E+03, 0.74727E+03, 0.75330E+03 &
, 0.75937E+03, 0.76547E+03, 0.77161E+03, 0.77779E+03, 0.78400E+03, 0.79024E+03, 0.79653E+03, 0.80285E+03, 0.80921E+03, 0.81560E+03 &
, 0.82203E+03, 0.82850E+03, 0.83501E+03, 0.84156E+03, 0.84814E+03, 0.85477E+03, 0.86143E+03, 0.86813E+03, 0.87487E+03, 0.88165E+03 &
, 0.88847E+03, 0.89533E+03, 0.90223E+03, 0.90917E+03, 0.91615E+03, 0.92317E+03, 0.93023E+03, 0.93733E+03, 0.94448E+03, 0.95166E+03 &
, 0.95889E+03, 0.96616E+03, 0.97348E+03, 0.98083E+03, 0.98823E+03, 0.99567E+03, 0.10032E+04, 0.10107E+04, 0.10183E+04, 0.10259E+04 &
, 0.10335E+04, 0.10412E+04, 0.10490E+04, 0.10568E+04, 0.10646E+04, 0.10725E+04, 0.10804E+04, 0.10884E+04, 0.10964E+04, 0.11045E+04 &
, 0.11126E+04, 0.11207E+04, 0.11290E+04, 0.11372E+04, 0.11455E+04, 0.11538E+04, 0.11622E+04, 0.11707E+04, 0.11792E+04, 0.11877E+04 &
, 0.11963E+04, 0.12049E+04, 0.12136E+04, 0.12224E+04, 0.12312E+04, 0.12400E+04, 0.12489E+04, 0.12578E+04, 0.12668E+04, 0.12758E+04 &
, 0.12849E+04, 0.12941E+04, 0.13033E+04, 0.13125E+04, 0.13218E+04, 0.13311E+04, 0.13405E+04, 0.13500E+04, 0.13595E+04, 0.13691E+04 &
, 0.13787E+04, 0.13883E+04, 0.13981E+04, 0.14078E+04, 0.14177E+04, 0.14275E+04, 0.14375E+04, 0.14475E+04, 0.14575E+04, 0.14676E+04 &
, 0.14778E+04, 0.14880E+04, 0.14983E+04, 0.15086E+04, 0.15190E+04, 0.15294E+04, 0.15399E+04, 0.15505E+04, 0.15611E+04, 0.15718E+04 &
, 0.15825E+04, 0.15933E+04, 0.16042E+04, 0.16151E+04, 0.16261E+04, 0.16371E+04, 0.16482E+04, 0.16594E+04, 0.16706E+04, 0.16818E+04 &
, 0.16932E+04, 0.17046E+04, 0.17160E+04, 0.17276E+04, 0.17392E+04, 0.17508E+04, 0.17625E+04, 0.17743E+04, 0.17861E+04, 0.17980E+04 &
, 0.18100E+04, 0.18220E+04, 0.18341E+04, 0.18463E+04, 0.18585E+04, 0.18708E+04, 0.18831E+04, 0.18956E+04, 0.19081E+04, 0.19206E+04 &
, 0.19332E+04, 0.19459E+04, 0.19587E+04, 0.19715E+04, 0.19844E+04, 0.19973E+04, 0.20104E+04, 0.20235E+04, 0.20366E+04, 0.20499E+04 &
, 0.20632E+04, 0.20766E+04, 0.20900E+04, 0.21035E+04, 0.21171E+04, 0.21308E+04, 0.21445E+04, 0.21583E+04, 0.21722E+04, 0.21861E+04 &
, 0.22001E+04, 0.22142E+04, 0.22284E+04, 0.22427E+04, 0.22570E+04, 0.22714E+04, 0.22858E+04, 0.23004E+04, 0.23150E+04, 0.23297E+04 &
, 0.23444E+04, 0.23593E+04, 0.23742E+04, 0.23892E+04, 0.24043E+04, 0.24194E+04, 0.24346E+04, 0.24499E+04, 0.24653E+04, 0.24808E+04 &
, 0.24963E+04, 0.25119E+04, 0.25276E+04, 0.25434E+04, 0.25593E+04, 0.25752E+04, 0.25912E+04, 0.26074E+04, 0.26235E+04, 0.26398E+04 &
, 0.26562E+04, 0.26726E+04, 0.26891E+04, 0.27057E+04, 0.27224E+04, 0.27392E+04, 0.27560E+04, 0.27729E+04, 0.27900E+04, 0.28071E+04 &
, 0.28243E+04, 0.28415E+04, 0.28589E+04, 0.28763E+04, 0.28939E+04, 0.29115E+04, 0.29292E+04, 0.29470E+04, 0.29649E+04, 0.29829E+04 &
, 0.30010E+04, 0.30191E+04, 0.30374E+04, 0.30557E+04, 0.30741E+04, 0.30926E+04, 0.31112E+04, 0.31299E+04, 0.31487E+04, 0.31676E+04 &
, 0.31866E+04, 0.32057E+04, 0.32248E+04, 0.32441E+04, 0.32634E+04, 0.32829E+04, 0.33024E+04, 0.33221E+04, 0.33418E+04, 0.33616E+04 &
, 0.33815E+04, 0.34016E+04, 0.34217E+04, 0.34419E+04, 0.34622E+04, 0.34826E+04, 0.35031E+04, 0.35237E+04, 0.35444E+04, 0.35652E+04 &
, 0.35861E+04, 0.36072E+04, 0.36283E+04, 0.36495E+04, 0.36708E+04, 0.36922E+04, 0.37137E+04, 0.37353E+04, 0.37570E+04, 0.37789E+04 &
, 0.38008E+04, 0.38228E+04, 0.38450E+04, 0.38672E+04, 0.38895E+04, 0.39120E+04, 0.39346E+04, 0.39572E+04, 0.39800E+04, 0.40029E+04 &
, 0.40259E+04, 0.40490E+04, 0.40722E+04, 0.40955E+04, 0.41189E+04, 0.41424E+04, 0.41661E+04, 0.41898E+04, 0.42137E+04, 0.42377E+04 &
, 0.42617E+04, 0.42859E+04, 0.43103E+04, 0.43347E+04, 0.43592E+04, 0.43839E+04, 0.44086E+04, 0.44335E+04, 0.44585E+04, 0.44836E+04 &
, 0.45088E+04, 0.45342E+04, 0.45596E+04, 0.45852E+04, 0.46109E+04, 0.46367E+04, 0.46627E+04, 0.46887E+04, 0.47149E+04, 0.47412E+04 &
, 0.47676E+04, 0.47941E+04, 0.48208E+04, 0.48475E+04, 0.48744E+04, 0.49014E+04, 0.49286E+04, 0.49559E+04, 0.49832E+04, 0.50107E+04 &
, 0.50384E+04, 0.50661E+04, 0.50940E+04, 0.51220E+04, 0.51502E+04, 0.51784E+04, 0.52068E+04, 0.52353E+04, 0.52640E+04, 0.52928E+04 &
, 0.53217E+04, 0.53507E+04, 0.53799E+04, 0.54092E+04, 0.54386E+04, 0.54681E+04, 0.54978E+04, 0.55277E+04, 0.55576E+04, 0.55877E+04 &
, 0.56179E+04, 0.56483E+04, 0.56787E+04, 0.57094E+04, 0.57401E+04, 0.57710E+04, 0.58020E+04, 0.58332E+04, 0.58645E+04, 0.58960E+04 &
, 0.59275E+04, 0.59593E+04, 0.59911E+04, 0.60231E+04, 0.60552E+04, 0.60875E+04, 0.61199E+04, 0.61525E+04, 0.61852E+04, 0.62180E+04 /
 !            2           3
 DATA(QofT(          12 ,J),J=1,501)/   38.0909996032715       &
, 0.41866E+02, 0.45642E+02, 0.49417E+02, 0.53193E+02, 0.56969E+02, 0.60744E+02, 0.64520E+02, 0.68296E+02, 0.72071E+02, 0.75847E+02 &
, 0.79623E+02, 0.83399E+02, 0.87175E+02, 0.90951E+02, 0.94727E+02, 0.98502E+02, 0.10228E+03, 0.10605E+03, 0.10983E+03, 0.11361E+03 &
, 0.11738E+03, 0.12116E+03, 0.12493E+03, 0.12871E+03, 0.13249E+03, 0.13626E+03, 0.14004E+03, 0.14382E+03, 0.14759E+03, 0.15137E+03 &
, 0.15515E+03, 0.15892E+03, 0.16270E+03, 0.16648E+03, 0.17026E+03, 0.17404E+03, 0.17781E+03, 0.18159E+03, 0.18538E+03, 0.18916E+03 &
, 0.19294E+03, 0.19672E+03, 0.20051E+03, 0.20430E+03, 0.20808E+03, 0.21187E+03, 0.21567E+03, 0.21946E+03, 0.22326E+03, 0.22706E+03 &
, 0.23086E+03, 0.23466E+03, 0.23847E+03, 0.24228E+03, 0.24610E+03, 0.24992E+03, 0.25374E+03, 0.25757E+03, 0.26141E+03, 0.26524E+03 &
, 0.26909E+03, 0.27294E+03, 0.27680E+03, 0.28066E+03, 0.28453E+03, 0.28841E+03, 0.29229E+03, 0.29619E+03, 0.30009E+03, 0.30400E+03 &
, 0.30791E+03, 0.31184E+03, 0.31578E+03, 0.31973E+03, 0.32368E+03, 0.32765E+03, 0.33163E+03, 0.33562E+03, 0.33962E+03, 0.34363E+03 &
, 0.34766E+03, 0.35170E+03, 0.35575E+03, 0.35981E+03, 0.36389E+03, 0.36799E+03, 0.37209E+03, 0.37621E+03, 0.38035E+03, 0.38450E+03 &
, 0.38867E+03, 0.39285E+03, 0.39705E+03, 0.40127E+03, 0.40550E+03, 0.40975E+03, 0.41402E+03, 0.41831E+03, 0.42261E+03, 0.42693E+03 &
, 0.43128E+03, 0.43564E+03, 0.44002E+03, 0.44442E+03, 0.44884E+03, 0.45328E+03, 0.45774E+03, 0.46222E+03, 0.46673E+03, 0.47125E+03 &
, 0.47580E+03, 0.48037E+03, 0.48496E+03, 0.48957E+03, 0.49421E+03, 0.49887E+03, 0.50355E+03, 0.50826E+03, 0.51299E+03, 0.51775E+03 &
, 0.52253E+03, 0.52733E+03, 0.53216E+03, 0.53702E+03, 0.54190E+03, 0.54681E+03, 0.55174E+03, 0.55670E+03, 0.56169E+03, 0.56670E+03 &
, 0.57174E+03, 0.57681E+03, 0.58191E+03, 0.58703E+03, 0.59218E+03, 0.59736E+03, 0.60257E+03, 0.60781E+03, 0.61307E+03, 0.61837E+03 &
, 0.62370E+03, 0.62905E+03, 0.63444E+03, 0.63986E+03, 0.64530E+03, 0.65078E+03, 0.65629E+03, 0.66183E+03, 0.66740E+03, 0.67300E+03 &
, 0.67864E+03, 0.68431E+03, 0.69001E+03, 0.69574E+03, 0.70151E+03, 0.70731E+03, 0.71314E+03, 0.71900E+03, 0.72490E+03, 0.73084E+03 &
, 0.73681E+03, 0.74281E+03, 0.74885E+03, 0.75492E+03, 0.76103E+03, 0.76717E+03, 0.77335E+03, 0.77957E+03, 0.78582E+03, 0.79211E+03 &
, 0.79843E+03, 0.80479E+03, 0.81119E+03, 0.81762E+03, 0.82410E+03, 0.83061E+03, 0.83716E+03, 0.84374E+03, 0.85037E+03, 0.85703E+03 &
, 0.86374E+03, 0.87048E+03, 0.87726E+03, 0.88408E+03, 0.89094E+03, 0.89784E+03, 0.90478E+03, 0.91176E+03, 0.91878E+03, 0.92585E+03 &
, 0.93295E+03, 0.94010E+03, 0.94728E+03, 0.95451E+03, 0.96178E+03, 0.96910E+03, 0.97645E+03, 0.98385E+03, 0.99129E+03, 0.99878E+03 &
, 0.10063E+04, 0.10139E+04, 0.10215E+04, 0.10292E+04, 0.10369E+04, 0.10446E+04, 0.10524E+04, 0.10602E+04, 0.10681E+04, 0.10761E+04 &
, 0.10840E+04, 0.10920E+04, 0.11001E+04, 0.11082E+04, 0.11164E+04, 0.11246E+04, 0.11329E+04, 0.11412E+04, 0.11495E+04, 0.11579E+04 &
, 0.11663E+04, 0.11748E+04, 0.11834E+04, 0.11920E+04, 0.12006E+04, 0.12093E+04, 0.12180E+04, 0.12268E+04, 0.12357E+04, 0.12446E+04 &
, 0.12535E+04, 0.12625E+04, 0.12715E+04, 0.12806E+04, 0.12898E+04, 0.12990E+04, 0.13082E+04, 0.13175E+04, 0.13269E+04, 0.13363E+04 &
, 0.13457E+04, 0.13553E+04, 0.13648E+04, 0.13745E+04, 0.13841E+04, 0.13939E+04, 0.14036E+04, 0.14135E+04, 0.14234E+04, 0.14333E+04 &
, 0.14433E+04, 0.14534E+04, 0.14635E+04, 0.14737E+04, 0.14839E+04, 0.14942E+04, 0.15045E+04, 0.15149E+04, 0.15254E+04, 0.15359E+04 &
, 0.15465E+04, 0.15571E+04, 0.15678E+04, 0.15786E+04, 0.15894E+04, 0.16002E+04, 0.16112E+04, 0.16222E+04, 0.16332E+04, 0.16443E+04 &
, 0.16555E+04, 0.16667E+04, 0.16780E+04, 0.16894E+04, 0.17008E+04, 0.17123E+04, 0.17238E+04, 0.17354E+04, 0.17471E+04, 0.17588E+04 &
, 0.17706E+04, 0.17825E+04, 0.17944E+04, 0.18064E+04, 0.18184E+04, 0.18306E+04, 0.18427E+04, 0.18550E+04, 0.18673E+04, 0.18797E+04 &
, 0.18921E+04, 0.19046E+04, 0.19172E+04, 0.19299E+04, 0.19426E+04, 0.19554E+04, 0.19682E+04, 0.19811E+04, 0.19941E+04, 0.20072E+04 &
, 0.20203E+04, 0.20335E+04, 0.20468E+04, 0.20601E+04, 0.20735E+04, 0.20870E+04, 0.21006E+04, 0.21142E+04, 0.21279E+04, 0.21417E+04 &
, 0.21555E+04, 0.21694E+04, 0.21834E+04, 0.21974E+04, 0.22116E+04, 0.22258E+04, 0.22401E+04, 0.22544E+04, 0.22689E+04, 0.22834E+04 &
, 0.22979E+04, 0.23126E+04, 0.23273E+04, 0.23421E+04, 0.23570E+04, 0.23720E+04, 0.23870E+04, 0.24021E+04, 0.24173E+04, 0.24326E+04 &
, 0.24480E+04, 0.24634E+04, 0.24789E+04, 0.24945E+04, 0.25102E+04, 0.25259E+04, 0.25418E+04, 0.25577E+04, 0.25737E+04, 0.25897E+04 &
, 0.26059E+04, 0.26221E+04, 0.26385E+04, 0.26549E+04, 0.26714E+04, 0.26879E+04, 0.27046E+04, 0.27213E+04, 0.27381E+04, 0.27551E+04 &
, 0.27721E+04, 0.27891E+04, 0.28063E+04, 0.28236E+04, 0.28409E+04, 0.28583E+04, 0.28758E+04, 0.28934E+04, 0.29111E+04, 0.29289E+04 &
, 0.29468E+04, 0.29647E+04, 0.29828E+04, 0.30009E+04, 0.30192E+04, 0.30375E+04, 0.30559E+04, 0.30744E+04, 0.30930E+04, 0.31117E+04 &
, 0.31304E+04, 0.31493E+04, 0.31683E+04, 0.31873E+04, 0.32065E+04, 0.32257E+04, 0.32450E+04, 0.32645E+04, 0.32840E+04, 0.33036E+04 &
, 0.33233E+04, 0.33432E+04, 0.33631E+04, 0.33831E+04, 0.34032E+04, 0.34234E+04, 0.34437E+04, 0.34641E+04, 0.34846E+04, 0.35052E+04 &
, 0.35259E+04, 0.35467E+04, 0.35676E+04, 0.35886E+04, 0.36097E+04, 0.36309E+04, 0.36522E+04, 0.36736E+04, 0.36951E+04, 0.37167E+04 &
, 0.37385E+04, 0.37603E+04, 0.37822E+04, 0.38042E+04, 0.38264E+04, 0.38486E+04, 0.38710E+04, 0.38934E+04, 0.39160E+04, 0.39387E+04 &
, 0.39614E+04, 0.39843E+04, 0.40073E+04, 0.40304E+04, 0.40536E+04, 0.40769E+04, 0.41004E+04, 0.41239E+04, 0.41475E+04, 0.41713E+04 &
, 0.41952E+04, 0.42192E+04, 0.42433E+04, 0.42675E+04, 0.42918E+04, 0.43162E+04, 0.43408E+04, 0.43655E+04, 0.43902E+04, 0.44151E+04 &
, 0.44401E+04, 0.44653E+04, 0.44905E+04, 0.45159E+04, 0.45414E+04, 0.45669E+04, 0.45927E+04, 0.46185E+04, 0.46444E+04, 0.46705E+04 &
, 0.46967E+04, 0.47230E+04, 0.47495E+04, 0.47760E+04, 0.48027E+04, 0.48295E+04, 0.48564E+04, 0.48835E+04, 0.49106E+04, 0.49379E+04 &
, 0.49653E+04, 0.49929E+04, 0.50205E+04, 0.50483E+04, 0.50762E+04, 0.51043E+04, 0.51325E+04, 0.51608E+04, 0.51892E+04, 0.52177E+04 &
, 0.52464E+04, 0.52752E+04, 0.53042E+04, 0.53333E+04, 0.53625E+04, 0.53918E+04, 0.54213E+04, 0.54509E+04, 0.54806E+04, 0.55105E+04 &
, 0.55405E+04, 0.55706E+04, 0.56009E+04, 0.56313E+04, 0.56618E+04, 0.56925E+04, 0.57233E+04, 0.57542E+04, 0.57853E+04, 0.58165E+04 &
, 0.58479E+04, 0.58794E+04, 0.59110E+04, 0.59428E+04, 0.59747E+04, 0.60067E+04, 0.60389E+04, 0.60713E+04, 0.61038E+04, 0.61364E+04 &
, 0.61691E+04, 0.62021E+04, 0.62351E+04, 0.62683E+04, 0.63016E+04, 0.63351E+04, 0.63688E+04, 0.64025E+04, 0.64365E+04, 0.64705E+04 /
 !            2           4
 DATA(QofT(          13 ,J),J=1,501)/   222.309997558594       &
, 0.24433E+03, 0.26636E+03, 0.28839E+03, 0.31042E+03, 0.33245E+03, 0.35448E+03, 0.37651E+03, 0.39854E+03, 0.42057E+03, 0.44261E+03 &
, 0.46464E+03, 0.48667E+03, 0.50870E+03, 0.53073E+03, 0.55276E+03, 0.57479E+03, 0.59682E+03, 0.61886E+03, 0.64089E+03, 0.66292E+03 &
, 0.68495E+03, 0.70698E+03, 0.72902E+03, 0.75105E+03, 0.77308E+03, 0.79512E+03, 0.81715E+03, 0.83919E+03, 0.86122E+03, 0.88326E+03 &
, 0.90529E+03, 0.92733E+03, 0.94937E+03, 0.97142E+03, 0.99346E+03, 0.10155E+04, 0.10376E+04, 0.10596E+04, 0.10817E+04, 0.11037E+04 &
, 0.11258E+04, 0.11479E+04, 0.11700E+04, 0.11921E+04, 0.12142E+04, 0.12363E+04, 0.12584E+04, 0.12805E+04, 0.13027E+04, 0.13248E+04 &
, 0.13470E+04, 0.13692E+04, 0.13914E+04, 0.14137E+04, 0.14359E+04, 0.14582E+04, 0.14805E+04, 0.15029E+04, 0.15252E+04, 0.15476E+04 &
, 0.15700E+04, 0.15925E+04, 0.16150E+04, 0.16375E+04, 0.16601E+04, 0.16827E+04, 0.17054E+04, 0.17281E+04, 0.17508E+04, 0.17736E+04 &
, 0.17965E+04, 0.18194E+04, 0.18423E+04, 0.18653E+04, 0.18884E+04, 0.19115E+04, 0.19347E+04, 0.19580E+04, 0.19813E+04, 0.20047E+04 &
, 0.20282E+04, 0.20517E+04, 0.20753E+04, 0.20990E+04, 0.21228E+04, 0.21466E+04, 0.21706E+04, 0.21946E+04, 0.22187E+04, 0.22429E+04 &
, 0.22672E+04, 0.22915E+04, 0.23160E+04, 0.23406E+04, 0.23652E+04, 0.23900E+04, 0.24148E+04, 0.24398E+04, 0.24649E+04, 0.24901E+04 &
, 0.25154E+04, 0.25408E+04, 0.25663E+04, 0.25919E+04, 0.26176E+04, 0.26435E+04, 0.26695E+04, 0.26956E+04, 0.27218E+04, 0.27481E+04 &
, 0.27746E+04, 0.28012E+04, 0.28279E+04, 0.28548E+04, 0.28818E+04, 0.29089E+04, 0.29361E+04, 0.29635E+04, 0.29911E+04, 0.30188E+04 &
, 0.30466E+04, 0.30745E+04, 0.31026E+04, 0.31309E+04, 0.31593E+04, 0.31878E+04, 0.32165E+04, 0.32454E+04, 0.32744E+04, 0.33036E+04 &
, 0.33329E+04, 0.33623E+04, 0.33920E+04, 0.34218E+04, 0.34517E+04, 0.34819E+04, 0.35122E+04, 0.35426E+04, 0.35732E+04, 0.36040E+04 &
, 0.36350E+04, 0.36661E+04, 0.36974E+04, 0.37289E+04, 0.37606E+04, 0.37924E+04, 0.38245E+04, 0.38567E+04, 0.38890E+04, 0.39216E+04 &
, 0.39543E+04, 0.39873E+04, 0.40204E+04, 0.40537E+04, 0.40872E+04, 0.41209E+04, 0.41548E+04, 0.41889E+04, 0.42232E+04, 0.42576E+04 &
, 0.42923E+04, 0.43272E+04, 0.43622E+04, 0.43975E+04, 0.44330E+04, 0.44687E+04, 0.45046E+04, 0.45406E+04, 0.45769E+04, 0.46135E+04 &
, 0.46502E+04, 0.46871E+04, 0.47243E+04, 0.47616E+04, 0.47992E+04, 0.48370E+04, 0.48750E+04, 0.49132E+04, 0.49517E+04, 0.49904E+04 &
, 0.50293E+04, 0.50684E+04, 0.51078E+04, 0.51474E+04, 0.51872E+04, 0.52272E+04, 0.52675E+04, 0.53080E+04, 0.53487E+04, 0.53897E+04 &
, 0.54309E+04, 0.54724E+04, 0.55141E+04, 0.55560E+04, 0.55982E+04, 0.56406E+04, 0.56833E+04, 0.57262E+04, 0.57694E+04, 0.58128E+04 &
, 0.58564E+04, 0.59003E+04, 0.59445E+04, 0.59889E+04, 0.60336E+04, 0.60785E+04, 0.61237E+04, 0.61692E+04, 0.62149E+04, 0.62608E+04 &
, 0.63071E+04, 0.63536E+04, 0.64004E+04, 0.64474E+04, 0.64947E+04, 0.65423E+04, 0.65901E+04, 0.66382E+04, 0.66866E+04, 0.67353E+04 &
, 0.67842E+04, 0.68335E+04, 0.68830E+04, 0.69328E+04, 0.69828E+04, 0.70332E+04, 0.70838E+04, 0.71348E+04, 0.71860E+04, 0.72375E+04 &
, 0.72893E+04, 0.73414E+04, 0.73938E+04, 0.74464E+04, 0.74994E+04, 0.75527E+04, 0.76063E+04, 0.76601E+04, 0.77143E+04, 0.77688E+04 &
, 0.78236E+04, 0.78787E+04, 0.79341E+04, 0.79898E+04, 0.80458E+04, 0.81022E+04, 0.81588E+04, 0.82158E+04, 0.82731E+04, 0.83307E+04 &
, 0.83886E+04, 0.84468E+04, 0.85054E+04, 0.85643E+04, 0.86235E+04, 0.86830E+04, 0.87429E+04, 0.88031E+04, 0.88637E+04, 0.89245E+04 &
, 0.89857E+04, 0.90473E+04, 0.91092E+04, 0.91714E+04, 0.92339E+04, 0.92968E+04, 0.93601E+04, 0.94237E+04, 0.94876E+04, 0.95519E+04 &
, 0.96166E+04, 0.96816E+04, 0.97469E+04, 0.98126E+04, 0.98787E+04, 0.99451E+04, 0.10012E+05, 0.10079E+05, 0.10146E+05, 0.10214E+05 &
, 0.10283E+05, 0.10351E+05, 0.10420E+05, 0.10489E+05, 0.10559E+05, 0.10629E+05, 0.10700E+05, 0.10771E+05, 0.10842E+05, 0.10913E+05 &
, 0.10985E+05, 0.11058E+05, 0.11131E+05, 0.11204E+05, 0.11277E+05, 0.11351E+05, 0.11425E+05, 0.11500E+05, 0.11575E+05, 0.11651E+05 &
, 0.11727E+05, 0.11803E+05, 0.11880E+05, 0.11957E+05, 0.12034E+05, 0.12112E+05, 0.12191E+05, 0.12269E+05, 0.12348E+05, 0.12428E+05 &
, 0.12508E+05, 0.12588E+05, 0.12669E+05, 0.12751E+05, 0.12832E+05, 0.12914E+05, 0.12997E+05, 0.13080E+05, 0.13163E+05, 0.13247E+05 &
, 0.13331E+05, 0.13416E+05, 0.13501E+05, 0.13587E+05, 0.13673E+05, 0.13759E+05, 0.13846E+05, 0.13933E+05, 0.14021E+05, 0.14109E+05 &
, 0.14198E+05, 0.14287E+05, 0.14377E+05, 0.14467E+05, 0.14557E+05, 0.14648E+05, 0.14740E+05, 0.14832E+05, 0.14924E+05, 0.15017E+05 &
, 0.15110E+05, 0.15204E+05, 0.15298E+05, 0.15393E+05, 0.15488E+05, 0.15584E+05, 0.15680E+05, 0.15777E+05, 0.15874E+05, 0.15972E+05 &
, 0.16070E+05, 0.16169E+05, 0.16268E+05, 0.16367E+05, 0.16467E+05, 0.16568E+05, 0.16669E+05, 0.16771E+05, 0.16873E+05, 0.16976E+05 &
, 0.17079E+05, 0.17183E+05, 0.17287E+05, 0.17391E+05, 0.17497E+05, 0.17602E+05, 0.17709E+05, 0.17815E+05, 0.17923E+05, 0.18031E+05 &
, 0.18139E+05, 0.18248E+05, 0.18357E+05, 0.18467E+05, 0.18578E+05, 0.18689E+05, 0.18800E+05, 0.18913E+05, 0.19025E+05, 0.19139E+05 &
, 0.19252E+05, 0.19367E+05, 0.19482E+05, 0.19597E+05, 0.19713E+05, 0.19830E+05, 0.19947E+05, 0.20065E+05, 0.20183E+05, 0.20302E+05 &
, 0.20421E+05, 0.20541E+05, 0.20662E+05, 0.20783E+05, 0.20905E+05, 0.21027E+05, 0.21150E+05, 0.21273E+05, 0.21397E+05, 0.21522E+05 &
, 0.21647E+05, 0.21773E+05, 0.21900E+05, 0.22027E+05, 0.22155E+05, 0.22283E+05, 0.22412E+05, 0.22541E+05, 0.22671E+05, 0.22802E+05 &
, 0.22934E+05, 0.23065E+05, 0.23198E+05, 0.23331E+05, 0.23465E+05, 0.23600E+05, 0.23735E+05, 0.23871E+05, 0.24007E+05, 0.24144E+05 &
, 0.24282E+05, 0.24420E+05, 0.24559E+05, 0.24698E+05, 0.24839E+05, 0.24980E+05, 0.25121E+05, 0.25263E+05, 0.25406E+05, 0.25550E+05 &
, 0.25694E+05, 0.25839E+05, 0.25984E+05, 0.26131E+05, 0.26278E+05, 0.26425E+05, 0.26573E+05, 0.26722E+05, 0.26872E+05, 0.27022E+05 &
, 0.27173E+05, 0.27325E+05, 0.27477E+05, 0.27630E+05, 0.27784E+05, 0.27938E+05, 0.28094E+05, 0.28250E+05, 0.28406E+05, 0.28563E+05 &
, 0.28721E+05, 0.28880E+05, 0.29040E+05, 0.29200E+05, 0.29361E+05, 0.29522E+05, 0.29685E+05, 0.29848E+05, 0.30012E+05, 0.30176E+05 &
, 0.30341E+05, 0.30507E+05, 0.30674E+05, 0.30842E+05, 0.31010E+05, 0.31179E+05, 0.31349E+05, 0.31519E+05, 0.31691E+05, 0.31863E+05 &
, 0.32036E+05, 0.32209E+05, 0.32384E+05, 0.32559E+05, 0.32735E+05, 0.32911E+05, 0.33089E+05, 0.33267E+05, 0.33446E+05, 0.33626E+05 &
, 0.33807E+05, 0.33988E+05, 0.34170E+05, 0.34353E+05, 0.34537E+05, 0.34722E+05, 0.34907E+05, 0.35094E+05, 0.35281E+05, 0.35469E+05 &
, 0.35657E+05, 0.35847E+05, 0.36037E+05, 0.36229E+05, 0.36421E+05, 0.36613E+05, 0.36807E+05, 0.37002E+05, 0.37197E+05, 0.37393E+05 /
 !            2           5
 DATA(QofT(          14 ,J),J=1,501)/   76.1809997558594       &
, 0.83732E+02, 0.91283E+02, 0.98834E+02, 0.10639E+03, 0.11394E+03, 0.12149E+03, 0.12904E+03, 0.13659E+03, 0.14414E+03, 0.15169E+03 &
, 0.15925E+03, 0.16680E+03, 0.17435E+03, 0.18190E+03, 0.18945E+03, 0.19700E+03, 0.20456E+03, 0.21211E+03, 0.21966E+03, 0.22721E+03 &
, 0.23476E+03, 0.24232E+03, 0.24987E+03, 0.25742E+03, 0.26497E+03, 0.27252E+03, 0.28008E+03, 0.28763E+03, 0.29518E+03, 0.30274E+03 &
, 0.31029E+03, 0.31785E+03, 0.32540E+03, 0.33296E+03, 0.34052E+03, 0.34808E+03, 0.35564E+03, 0.36320E+03, 0.37076E+03, 0.37833E+03 &
, 0.38590E+03, 0.39347E+03, 0.40104E+03, 0.40862E+03, 0.41620E+03, 0.42379E+03, 0.43138E+03, 0.43898E+03, 0.44658E+03, 0.45419E+03 &
, 0.46181E+03, 0.46943E+03, 0.47706E+03, 0.48470E+03, 0.49234E+03, 0.50000E+03, 0.50767E+03, 0.51535E+03, 0.52304E+03, 0.53074E+03 &
, 0.53845E+03, 0.54618E+03, 0.55392E+03, 0.56168E+03, 0.56945E+03, 0.57724E+03, 0.58505E+03, 0.59287E+03, 0.60071E+03, 0.60858E+03 &
, 0.61646E+03, 0.62436E+03, 0.63228E+03, 0.64023E+03, 0.64819E+03, 0.65619E+03, 0.66420E+03, 0.67224E+03, 0.68031E+03, 0.68840E+03 &
, 0.69652E+03, 0.70467E+03, 0.71285E+03, 0.72105E+03, 0.72929E+03, 0.73755E+03, 0.74585E+03, 0.75418E+03, 0.76254E+03, 0.77094E+03 &
, 0.77937E+03, 0.78783E+03, 0.79633E+03, 0.80487E+03, 0.81344E+03, 0.82205E+03, 0.83070E+03, 0.83939E+03, 0.84811E+03, 0.85688E+03 &
, 0.86569E+03, 0.87453E+03, 0.88342E+03, 0.89236E+03, 0.90133E+03, 0.91035E+03, 0.91942E+03, 0.92853E+03, 0.93768E+03, 0.94688E+03 &
, 0.95613E+03, 0.96542E+03, 0.97476E+03, 0.98415E+03, 0.99359E+03, 0.10031E+04, 0.10126E+04, 0.10222E+04, 0.10319E+04, 0.10415E+04 &
, 0.10513E+04, 0.10611E+04, 0.10709E+04, 0.10808E+04, 0.10908E+04, 0.11008E+04, 0.11109E+04, 0.11210E+04, 0.11312E+04, 0.11414E+04 &
, 0.11517E+04, 0.11621E+04, 0.11725E+04, 0.11830E+04, 0.11935E+04, 0.12041E+04, 0.12147E+04, 0.12255E+04, 0.12362E+04, 0.12471E+04 &
, 0.12580E+04, 0.12689E+04, 0.12800E+04, 0.12911E+04, 0.13022E+04, 0.13134E+04, 0.13247E+04, 0.13361E+04, 0.13475E+04, 0.13590E+04 &
, 0.13705E+04, 0.13821E+04, 0.13938E+04, 0.14056E+04, 0.14174E+04, 0.14293E+04, 0.14413E+04, 0.14533E+04, 0.14654E+04, 0.14776E+04 &
, 0.14899E+04, 0.15022E+04, 0.15146E+04, 0.15271E+04, 0.15396E+04, 0.15522E+04, 0.15649E+04, 0.15777E+04, 0.15905E+04, 0.16035E+04 &
, 0.16165E+04, 0.16295E+04, 0.16427E+04, 0.16559E+04, 0.16692E+04, 0.16826E+04, 0.16961E+04, 0.17097E+04, 0.17233E+04, 0.17370E+04 &
, 0.17508E+04, 0.17647E+04, 0.17786E+04, 0.17927E+04, 0.18068E+04, 0.18210E+04, 0.18353E+04, 0.18497E+04, 0.18642E+04, 0.18787E+04 &
, 0.18933E+04, 0.19081E+04, 0.19229E+04, 0.19378E+04, 0.19528E+04, 0.19678E+04, 0.19830E+04, 0.19983E+04, 0.20136E+04, 0.20290E+04 &
, 0.20446E+04, 0.20602E+04, 0.20759E+04, 0.20917E+04, 0.21076E+04, 0.21236E+04, 0.21397E+04, 0.21558E+04, 0.21721E+04, 0.21885E+04 &
, 0.22050E+04, 0.22215E+04, 0.22382E+04, 0.22549E+04, 0.22718E+04, 0.22887E+04, 0.23058E+04, 0.23229E+04, 0.23402E+04, 0.23575E+04 &
, 0.23750E+04, 0.23925E+04, 0.24102E+04, 0.24280E+04, 0.24458E+04, 0.24638E+04, 0.24818E+04, 0.25000E+04, 0.25183E+04, 0.25367E+04 &
, 0.25552E+04, 0.25738E+04, 0.25925E+04, 0.26113E+04, 0.26302E+04, 0.26492E+04, 0.26684E+04, 0.26876E+04, 0.27070E+04, 0.27264E+04 &
, 0.27460E+04, 0.27657E+04, 0.27855E+04, 0.28054E+04, 0.28255E+04, 0.28456E+04, 0.28659E+04, 0.28862E+04, 0.29067E+04, 0.29273E+04 &
, 0.29481E+04, 0.29689E+04, 0.29899E+04, 0.30109E+04, 0.30321E+04, 0.30535E+04, 0.30749E+04, 0.30964E+04, 0.31181E+04, 0.31399E+04 &
, 0.31618E+04, 0.31839E+04, 0.32060E+04, 0.32283E+04, 0.32508E+04, 0.32733E+04, 0.32960E+04, 0.33188E+04, 0.33417E+04, 0.33647E+04 &
, 0.33879E+04, 0.34112E+04, 0.34346E+04, 0.34582E+04, 0.34819E+04, 0.35057E+04, 0.35297E+04, 0.35538E+04, 0.35780E+04, 0.36023E+04 &
, 0.36268E+04, 0.36514E+04, 0.36762E+04, 0.37011E+04, 0.37261E+04, 0.37513E+04, 0.37766E+04, 0.38020E+04, 0.38276E+04, 0.38533E+04 &
, 0.38792E+04, 0.39052E+04, 0.39313E+04, 0.39576E+04, 0.39840E+04, 0.40106E+04, 0.40373E+04, 0.40642E+04, 0.40912E+04, 0.41183E+04 &
, 0.41456E+04, 0.41731E+04, 0.42006E+04, 0.42284E+04, 0.42563E+04, 0.42843E+04, 0.43125E+04, 0.43408E+04, 0.43693E+04, 0.43979E+04 &
, 0.44267E+04, 0.44557E+04, 0.44848E+04, 0.45140E+04, 0.45434E+04, 0.45730E+04, 0.46027E+04, 0.46326E+04, 0.46626E+04, 0.46928E+04 &
, 0.47231E+04, 0.47537E+04, 0.47843E+04, 0.48152E+04, 0.48462E+04, 0.48773E+04, 0.49086E+04, 0.49401E+04, 0.49718E+04, 0.50036E+04 &
, 0.50356E+04, 0.50677E+04, 0.51000E+04, 0.51325E+04, 0.51652E+04, 0.51980E+04, 0.52310E+04, 0.52641E+04, 0.52975E+04, 0.53310E+04 &
, 0.53646E+04, 0.53985E+04, 0.54325E+04, 0.54667E+04, 0.55011E+04, 0.55356E+04, 0.55704E+04, 0.56053E+04, 0.56404E+04, 0.56756E+04 &
, 0.57111E+04, 0.57467E+04, 0.57825E+04, 0.58185E+04, 0.58547E+04, 0.58910E+04, 0.59276E+04, 0.59643E+04, 0.60012E+04, 0.60383E+04 &
, 0.60756E+04, 0.61131E+04, 0.61508E+04, 0.61886E+04, 0.62267E+04, 0.62649E+04, 0.63033E+04, 0.63419E+04, 0.63808E+04, 0.64198E+04 &
, 0.64590E+04, 0.64984E+04, 0.65380E+04, 0.65778E+04, 0.66177E+04, 0.66579E+04, 0.66983E+04, 0.67389E+04, 0.67797E+04, 0.68207E+04 &
, 0.68619E+04, 0.69033E+04, 0.69449E+04, 0.69867E+04, 0.70287E+04, 0.70710E+04, 0.71134E+04, 0.71560E+04, 0.71989E+04, 0.72419E+04 &
, 0.72852E+04, 0.73287E+04, 0.73724E+04, 0.74163E+04, 0.74604E+04, 0.75048E+04, 0.75493E+04, 0.75941E+04, 0.76391E+04, 0.76843E+04 &
, 0.77297E+04, 0.77753E+04, 0.78212E+04, 0.78673E+04, 0.79136E+04, 0.79601E+04, 0.80069E+04, 0.80539E+04, 0.81011E+04, 0.81485E+04 &
, 0.81962E+04, 0.82441E+04, 0.82922E+04, 0.83405E+04, 0.83891E+04, 0.84379E+04, 0.84870E+04, 0.85362E+04, 0.85858E+04, 0.86355E+04 &
, 0.86855E+04, 0.87357E+04, 0.87862E+04, 0.88369E+04, 0.88878E+04, 0.89390E+04, 0.89904E+04, 0.90421E+04, 0.90940E+04, 0.91461E+04 &
, 0.91985E+04, 0.92512E+04, 0.93041E+04, 0.93572E+04, 0.94106E+04, 0.94642E+04, 0.95181E+04, 0.95722E+04, 0.96266E+04, 0.96813E+04 &
, 0.97362E+04, 0.97913E+04, 0.98467E+04, 0.99024E+04, 0.99583E+04, 0.10014E+05, 0.10071E+05, 0.10128E+05, 0.10185E+05, 0.10242E+05 &
, 0.10299E+05, 0.10357E+05, 0.10415E+05, 0.10473E+05, 0.10532E+05, 0.10591E+05, 0.10650E+05, 0.10709E+05, 0.10769E+05, 0.10829E+05 &
, 0.10889E+05, 0.10949E+05, 0.11010E+05, 0.11071E+05, 0.11132E+05, 0.11194E+05, 0.11256E+05, 0.11318E+05, 0.11380E+05, 0.11443E+05 &
, 0.11506E+05, 0.11569E+05, 0.11633E+05, 0.11697E+05, 0.11761E+05, 0.11825E+05, 0.11890E+05, 0.11955E+05, 0.12020E+05, 0.12086E+05 &
, 0.12151E+05, 0.12218E+05, 0.12284E+05, 0.12351E+05, 0.12418E+05, 0.12485E+05, 0.12553E+05, 0.12621E+05, 0.12689E+05, 0.12757E+05 &
, 0.12826E+05, 0.12895E+05, 0.12965E+05, 0.13035E+05, 0.13105E+05, 0.13175E+05, 0.13246E+05, 0.13317E+05, 0.13388E+05, 0.13460E+05 /
 !            2           6
 DATA(QofT(          15 ,J),J=1,501)/   444.589996337891       &
, 0.48865E+03, 0.53271E+03, 0.57676E+03, 0.62082E+03, 0.66488E+03, 0.70894E+03, 0.75300E+03, 0.79705E+03, 0.84111E+03, 0.88517E+03 &
, 0.92923E+03, 0.97329E+03, 0.10174E+04, 0.10614E+04, 0.11055E+04, 0.11495E+04, 0.11936E+04, 0.12377E+04, 0.12817E+04, 0.13258E+04 &
, 0.13698E+04, 0.14139E+04, 0.14580E+04, 0.15020E+04, 0.15461E+04, 0.15902E+04, 0.16342E+04, 0.16783E+04, 0.17224E+04, 0.17664E+04 &
, 0.18105E+04, 0.18546E+04, 0.18987E+04, 0.19428E+04, 0.19869E+04, 0.20310E+04, 0.20751E+04, 0.21192E+04, 0.21633E+04, 0.22075E+04 &
, 0.22516E+04, 0.22958E+04, 0.23400E+04, 0.23842E+04, 0.24284E+04, 0.24727E+04, 0.25170E+04, 0.25613E+04, 0.26057E+04, 0.26500E+04 &
, 0.26945E+04, 0.27389E+04, 0.27834E+04, 0.28280E+04, 0.28726E+04, 0.29173E+04, 0.29620E+04, 0.30068E+04, 0.30516E+04, 0.30965E+04 &
, 0.31415E+04, 0.31866E+04, 0.32317E+04, 0.32770E+04, 0.33223E+04, 0.33677E+04, 0.34132E+04, 0.34588E+04, 0.35046E+04, 0.35504E+04 &
, 0.35963E+04, 0.36424E+04, 0.36886E+04, 0.37349E+04, 0.37814E+04, 0.38279E+04, 0.38747E+04, 0.39215E+04, 0.39685E+04, 0.40157E+04 &
, 0.40630E+04, 0.41105E+04, 0.41581E+04, 0.42060E+04, 0.42539E+04, 0.43021E+04, 0.43504E+04, 0.43990E+04, 0.44477E+04, 0.44966E+04 &
, 0.45457E+04, 0.45950E+04, 0.46445E+04, 0.46942E+04, 0.47442E+04, 0.47943E+04, 0.48447E+04, 0.48953E+04, 0.49461E+04, 0.49971E+04 &
, 0.50484E+04, 0.50999E+04, 0.51517E+04, 0.52037E+04, 0.52559E+04, 0.53084E+04, 0.53612E+04, 0.54142E+04, 0.54675E+04, 0.55210E+04 &
, 0.55749E+04, 0.56290E+04, 0.56833E+04, 0.57380E+04, 0.57929E+04, 0.58481E+04, 0.59036E+04, 0.59594E+04, 0.60155E+04, 0.60719E+04 &
, 0.61285E+04, 0.61855E+04, 0.62428E+04, 0.63005E+04, 0.63584E+04, 0.64166E+04, 0.64752E+04, 0.65341E+04, 0.65933E+04, 0.66528E+04 &
, 0.67127E+04, 0.67729E+04, 0.68334E+04, 0.68943E+04, 0.69555E+04, 0.70171E+04, 0.70790E+04, 0.71413E+04, 0.72040E+04, 0.72669E+04 &
, 0.73303E+04, 0.73940E+04, 0.74581E+04, 0.75225E+04, 0.75874E+04, 0.76526E+04, 0.77181E+04, 0.77841E+04, 0.78504E+04, 0.79171E+04 &
, 0.79843E+04, 0.80518E+04, 0.81197E+04, 0.81880E+04, 0.82566E+04, 0.83257E+04, 0.83952E+04, 0.84652E+04, 0.85355E+04, 0.86062E+04 &
, 0.86774E+04, 0.87489E+04, 0.88209E+04, 0.88933E+04, 0.89662E+04, 0.90395E+04, 0.91132E+04, 0.91873E+04, 0.92619E+04, 0.93369E+04 &
, 0.94124E+04, 0.94883E+04, 0.95646E+04, 0.96414E+04, 0.97187E+04, 0.97964E+04, 0.98746E+04, 0.99532E+04, 0.10032E+05, 0.10112E+05 &
, 0.10192E+05, 0.10272E+05, 0.10353E+05, 0.10435E+05, 0.10517E+05, 0.10599E+05, 0.10682E+05, 0.10766E+05, 0.10849E+05, 0.10934E+05 &
, 0.11019E+05, 0.11104E+05, 0.11190E+05, 0.11276E+05, 0.11363E+05, 0.11451E+05, 0.11539E+05, 0.11627E+05, 0.11716E+05, 0.11806E+05 &
, 0.11896E+05, 0.11986E+05, 0.12077E+05, 0.12169E+05, 0.12261E+05, 0.12354E+05, 0.12447E+05, 0.12541E+05, 0.12635E+05, 0.12730E+05 &
, 0.12825E+05, 0.12921E+05, 0.13018E+05, 0.13115E+05, 0.13212E+05, 0.13311E+05, 0.13410E+05, 0.13509E+05, 0.13609E+05, 0.13709E+05 &
, 0.13810E+05, 0.13912E+05, 0.14014E+05, 0.14117E+05, 0.14221E+05, 0.14325E+05, 0.14429E+05, 0.14535E+05, 0.14641E+05, 0.14747E+05 &
, 0.14854E+05, 0.14962E+05, 0.15070E+05, 0.15179E+05, 0.15289E+05, 0.15399E+05, 0.15510E+05, 0.15621E+05, 0.15733E+05, 0.15846E+05 &
, 0.15959E+05, 0.16073E+05, 0.16188E+05, 0.16303E+05, 0.16419E+05, 0.16536E+05, 0.16653E+05, 0.16771E+05, 0.16889E+05, 0.17009E+05 &
, 0.17129E+05, 0.17249E+05, 0.17370E+05, 0.17492E+05, 0.17615E+05, 0.17738E+05, 0.17862E+05, 0.17987E+05, 0.18113E+05, 0.18239E+05 &
, 0.18366E+05, 0.18493E+05, 0.18621E+05, 0.18750E+05, 0.18880E+05, 0.19010E+05, 0.19141E+05, 0.19273E+05, 0.19406E+05, 0.19539E+05 &
, 0.19673E+05, 0.19808E+05, 0.19944E+05, 0.20080E+05, 0.20217E+05, 0.20355E+05, 0.20493E+05, 0.20632E+05, 0.20772E+05, 0.20913E+05 &
, 0.21055E+05, 0.21197E+05, 0.21340E+05, 0.21484E+05, 0.21629E+05, 0.21775E+05, 0.21921E+05, 0.22068E+05, 0.22216E+05, 0.22364E+05 &
, 0.22514E+05, 0.22664E+05, 0.22815E+05, 0.22967E+05, 0.23120E+05, 0.23274E+05, 0.23428E+05, 0.23583E+05, 0.23739E+05, 0.23896E+05 &
, 0.24054E+05, 0.24212E+05, 0.24372E+05, 0.24532E+05, 0.24693E+05, 0.24855E+05, 0.25018E+05, 0.25182E+05, 0.25346E+05, 0.25512E+05 &
, 0.25678E+05, 0.25845E+05, 0.26013E+05, 0.26182E+05, 0.26352E+05, 0.26523E+05, 0.26695E+05, 0.26867E+05, 0.27041E+05, 0.27215E+05 &
, 0.27390E+05, 0.27567E+05, 0.27744E+05, 0.27922E+05, 0.28101E+05, 0.28281E+05, 0.28462E+05, 0.28644E+05, 0.28826E+05, 0.29010E+05 &
, 0.29195E+05, 0.29380E+05, 0.29567E+05, 0.29754E+05, 0.29943E+05, 0.30132E+05, 0.30323E+05, 0.30514E+05, 0.30707E+05, 0.30900E+05 &
, 0.31095E+05, 0.31290E+05, 0.31487E+05, 0.31684E+05, 0.31882E+05, 0.32082E+05, 0.32282E+05, 0.32484E+05, 0.32686E+05, 0.32890E+05 &
, 0.33095E+05, 0.33300E+05, 0.33507E+05, 0.33715E+05, 0.33923E+05, 0.34133E+05, 0.34344E+05, 0.34556E+05, 0.34769E+05, 0.34983E+05 &
, 0.35198E+05, 0.35415E+05, 0.35632E+05, 0.35850E+05, 0.36070E+05, 0.36291E+05, 0.36512E+05, 0.36735E+05, 0.36959E+05, 0.37184E+05 &
, 0.37410E+05, 0.37638E+05, 0.37866E+05, 0.38096E+05, 0.38326E+05, 0.38558E+05, 0.38791E+05, 0.39025E+05, 0.39261E+05, 0.39497E+05 &
, 0.39735E+05, 0.39974E+05, 0.40214E+05, 0.40455E+05, 0.40697E+05, 0.40941E+05, 0.41185E+05, 0.41431E+05, 0.41678E+05, 0.41927E+05 &
, 0.42176E+05, 0.42427E+05, 0.42679E+05, 0.42932E+05, 0.43187E+05, 0.43442E+05, 0.43699E+05, 0.43957E+05, 0.44217E+05, 0.44477E+05 &
, 0.44739E+05, 0.45002E+05, 0.45267E+05, 0.45533E+05, 0.45800E+05, 0.46068E+05, 0.46337E+05, 0.46608E+05, 0.46880E+05, 0.47154E+05 &
, 0.47429E+05, 0.47705E+05, 0.47982E+05, 0.48261E+05, 0.48541E+05, 0.48822E+05, 0.49105E+05, 0.49389E+05, 0.49674E+05, 0.49961E+05 &
, 0.50249E+05, 0.50539E+05, 0.50829E+05, 0.51122E+05, 0.51415E+05, 0.51710E+05, 0.52007E+05, 0.52304E+05, 0.52603E+05, 0.52904E+05 &
, 0.53206E+05, 0.53509E+05, 0.53814E+05, 0.54120E+05, 0.54428E+05, 0.54737E+05, 0.55047E+05, 0.55359E+05, 0.55673E+05, 0.55987E+05 &
, 0.56304E+05, 0.56622E+05, 0.56941E+05, 0.57261E+05, 0.57584E+05, 0.57907E+05, 0.58233E+05, 0.58559E+05, 0.58887E+05, 0.59217E+05 &
, 0.59548E+05, 0.59881E+05, 0.60215E+05, 0.60551E+05, 0.60888E+05, 0.61227E+05, 0.61567E+05, 0.61909E+05, 0.62253E+05, 0.62598E+05 &
, 0.62944E+05, 0.63293E+05, 0.63642E+05, 0.63994E+05, 0.64347E+05, 0.64701E+05, 0.65057E+05, 0.65415E+05, 0.65774E+05, 0.66135E+05 &
, 0.66498E+05, 0.66862E+05, 0.67228E+05, 0.67596E+05, 0.67965E+05, 0.68336E+05, 0.68708E+05, 0.69082E+05, 0.69458E+05, 0.69835E+05 &
, 0.70215E+05, 0.70595E+05, 0.70978E+05, 0.71362E+05, 0.71748E+05, 0.72136E+05, 0.72525E+05, 0.72916E+05, 0.73309E+05, 0.73704E+05 &
, 0.74100E+05, 0.74498E+05, 0.74898E+05, 0.75299E+05, 0.75703E+05, 0.76108E+05, 0.76515E+05, 0.76923E+05, 0.77334E+05, 0.77746E+05 /
 !            2           7
 DATA(QofT(          16 ,J),J=1,501)/   20.2089996337891       &
, 0.22213E+02, 0.24217E+02, 0.26221E+02, 0.28226E+02, 0.30230E+02, 0.32234E+02, 0.34238E+02, 0.36243E+02, 0.38247E+02, 0.40251E+02 &
, 0.42255E+02, 0.44260E+02, 0.46264E+02, 0.48268E+02, 0.50273E+02, 0.52277E+02, 0.54281E+02, 0.56286E+02, 0.58290E+02, 0.60294E+02 &
, 0.62299E+02, 0.64303E+02, 0.66308E+02, 0.68312E+02, 0.70317E+02, 0.72321E+02, 0.74326E+02, 0.76330E+02, 0.78335E+02, 0.80340E+02 &
, 0.82345E+02, 0.84350E+02, 0.86355E+02, 0.88360E+02, 0.90366E+02, 0.92372E+02, 0.94378E+02, 0.96385E+02, 0.98392E+02, 0.10040E+03 &
, 0.10241E+03, 0.10442E+03, 0.10643E+03, 0.10844E+03, 0.11045E+03, 0.11246E+03, 0.11447E+03, 0.11649E+03, 0.11850E+03, 0.12052E+03 &
, 0.12254E+03, 0.12456E+03, 0.12658E+03, 0.12861E+03, 0.13063E+03, 0.13266E+03, 0.13469E+03, 0.13673E+03, 0.13876E+03, 0.14080E+03 &
, 0.14285E+03, 0.14489E+03, 0.14694E+03, 0.14899E+03, 0.15105E+03, 0.15311E+03, 0.15518E+03, 0.15725E+03, 0.15932E+03, 0.16140E+03 &
, 0.16348E+03, 0.16557E+03, 0.16766E+03, 0.16976E+03, 0.17186E+03, 0.17397E+03, 0.17609E+03, 0.17821E+03, 0.18034E+03, 0.18248E+03 &
, 0.18462E+03, 0.18677E+03, 0.18892E+03, 0.19109E+03, 0.19326E+03, 0.19544E+03, 0.19762E+03, 0.19982E+03, 0.20202E+03, 0.20423E+03 &
, 0.20645E+03, 0.20868E+03, 0.21091E+03, 0.21316E+03, 0.21541E+03, 0.21768E+03, 0.21995E+03, 0.22224E+03, 0.22453E+03, 0.22683E+03 &
, 0.22915E+03, 0.23147E+03, 0.23381E+03, 0.23615E+03, 0.23851E+03, 0.24088E+03, 0.24326E+03, 0.24565E+03, 0.24805E+03, 0.25047E+03 &
, 0.25289E+03, 0.25533E+03, 0.25778E+03, 0.26024E+03, 0.26272E+03, 0.26520E+03, 0.26771E+03, 0.27022E+03, 0.27275E+03, 0.27528E+03 &
, 0.27784E+03, 0.28040E+03, 0.28298E+03, 0.28558E+03, 0.28819E+03, 0.29081E+03, 0.29345E+03, 0.29610E+03, 0.29876E+03, 0.30144E+03 &
, 0.30414E+03, 0.30685E+03, 0.30957E+03, 0.31231E+03, 0.31506E+03, 0.31784E+03, 0.32062E+03, 0.32342E+03, 0.32624E+03, 0.32908E+03 &
, 0.33193E+03, 0.33479E+03, 0.33767E+03, 0.34057E+03, 0.34349E+03, 0.34642E+03, 0.34937E+03, 0.35234E+03, 0.35532E+03, 0.35832E+03 &
, 0.36134E+03, 0.36438E+03, 0.36743E+03, 0.37050E+03, 0.37359E+03, 0.37670E+03, 0.37983E+03, 0.38297E+03, 0.38613E+03, 0.38931E+03 &
, 0.39251E+03, 0.39573E+03, 0.39897E+03, 0.40223E+03, 0.40550E+03, 0.40880E+03, 0.41211E+03, 0.41545E+03, 0.41880E+03, 0.42217E+03 &
, 0.42557E+03, 0.42898E+03, 0.43242E+03, 0.43587E+03, 0.43934E+03, 0.44284E+03, 0.44636E+03, 0.44989E+03, 0.45345E+03, 0.45703E+03 &
, 0.46063E+03, 0.46425E+03, 0.46789E+03, 0.47156E+03, 0.47524E+03, 0.47895E+03, 0.48268E+03, 0.48643E+03, 0.49021E+03, 0.49400E+03 &
, 0.49782E+03, 0.50166E+03, 0.50553E+03, 0.50941E+03, 0.51332E+03, 0.51725E+03, 0.52121E+03, 0.52519E+03, 0.52919E+03, 0.53322E+03 &
, 0.53727E+03, 0.54134E+03, 0.54544E+03, 0.54956E+03, 0.55371E+03, 0.55788E+03, 0.56207E+03, 0.56629E+03, 0.57054E+03, 0.57480E+03 &
, 0.57910E+03, 0.58342E+03, 0.58776E+03, 0.59213E+03, 0.59653E+03, 0.60095E+03, 0.60539E+03, 0.60986E+03, 0.61436E+03, 0.61889E+03 &
, 0.62344E+03, 0.62802E+03, 0.63262E+03, 0.63725E+03, 0.64191E+03, 0.64659E+03, 0.65130E+03, 0.65604E+03, 0.66080E+03, 0.66560E+03 &
, 0.67042E+03, 0.67527E+03, 0.68014E+03, 0.68505E+03, 0.68998E+03, 0.69494E+03, 0.69993E+03, 0.70495E+03, 0.70999E+03, 0.71507E+03 &
, 0.72017E+03, 0.72531E+03, 0.73047E+03, 0.73566E+03, 0.74088E+03, 0.74613E+03, 0.75142E+03, 0.75673E+03, 0.76207E+03, 0.76744E+03 &
, 0.77284E+03, 0.77827E+03, 0.78373E+03, 0.78923E+03, 0.79475E+03, 0.80031E+03, 0.80590E+03, 0.81151E+03, 0.81716E+03, 0.82284E+03 &
, 0.82856E+03, 0.83430E+03, 0.84008E+03, 0.84589E+03, 0.85173E+03, 0.85761E+03, 0.86351E+03, 0.86945E+03, 0.87543E+03, 0.88143E+03 &
, 0.88747E+03, 0.89354E+03, 0.89965E+03, 0.90579E+03, 0.91196E+03, 0.91817E+03, 0.92442E+03, 0.93069E+03, 0.93700E+03, 0.94335E+03 &
, 0.94973E+03, 0.95615E+03, 0.96260E+03, 0.96908E+03, 0.97560E+03, 0.98216E+03, 0.98875E+03, 0.99538E+03, 0.10020E+04, 0.10087E+04 &
, 0.10155E+04, 0.10223E+04, 0.10291E+04, 0.10359E+04, 0.10428E+04, 0.10497E+04, 0.10567E+04, 0.10637E+04, 0.10707E+04, 0.10778E+04 &
, 0.10849E+04, 0.10920E+04, 0.10992E+04, 0.11065E+04, 0.11137E+04, 0.11210E+04, 0.11284E+04, 0.11357E+04, 0.11432E+04, 0.11506E+04 &
, 0.11581E+04, 0.11657E+04, 0.11732E+04, 0.11809E+04, 0.11885E+04, 0.11962E+04, 0.12040E+04, 0.12117E+04, 0.12196E+04, 0.12274E+04 &
, 0.12353E+04, 0.12433E+04, 0.12513E+04, 0.12593E+04, 0.12674E+04, 0.12755E+04, 0.12836E+04, 0.12918E+04, 0.13001E+04, 0.13084E+04 &
, 0.13167E+04, 0.13250E+04, 0.13335E+04, 0.13419E+04, 0.13504E+04, 0.13590E+04, 0.13676E+04, 0.13762E+04, 0.13849E+04, 0.13936E+04 &
, 0.14024E+04, 0.14112E+04, 0.14200E+04, 0.14289E+04, 0.14379E+04, 0.14469E+04, 0.14559E+04, 0.14650E+04, 0.14742E+04, 0.14833E+04 &
, 0.14926E+04, 0.15018E+04, 0.15112E+04, 0.15205E+04, 0.15299E+04, 0.15394E+04, 0.15489E+04, 0.15585E+04, 0.15681E+04, 0.15777E+04 &
, 0.15875E+04, 0.15972E+04, 0.16070E+04, 0.16169E+04, 0.16268E+04, 0.16367E+04, 0.16467E+04, 0.16568E+04, 0.16669E+04, 0.16770E+04 &
, 0.16872E+04, 0.16975E+04, 0.17078E+04, 0.17181E+04, 0.17285E+04, 0.17390E+04, 0.17495E+04, 0.17601E+04, 0.17707E+04, 0.17814E+04 &
, 0.17921E+04, 0.18029E+04, 0.18137E+04, 0.18246E+04, 0.18355E+04, 0.18465E+04, 0.18575E+04, 0.18686E+04, 0.18798E+04, 0.18910E+04 &
, 0.19022E+04, 0.19135E+04, 0.19249E+04, 0.19363E+04, 0.19478E+04, 0.19593E+04, 0.19709E+04, 0.19826E+04, 0.19943E+04, 0.20061E+04 &
, 0.20179E+04, 0.20297E+04, 0.20417E+04, 0.20537E+04, 0.20657E+04, 0.20778E+04, 0.20900E+04, 0.21022E+04, 0.21145E+04, 0.21268E+04 &
, 0.21392E+04, 0.21517E+04, 0.21642E+04, 0.21768E+04, 0.21894E+04, 0.22021E+04, 0.22148E+04, 0.22277E+04, 0.22405E+04, 0.22535E+04 &
, 0.22665E+04, 0.22795E+04, 0.22927E+04, 0.23058E+04, 0.23191E+04, 0.23324E+04, 0.23458E+04, 0.23592E+04, 0.23727E+04, 0.23863E+04 &
, 0.23999E+04, 0.24136E+04, 0.24273E+04, 0.24411E+04, 0.24550E+04, 0.24690E+04, 0.24830E+04, 0.24971E+04, 0.25112E+04, 0.25254E+04 &
, 0.25397E+04, 0.25540E+04, 0.25684E+04, 0.25829E+04, 0.25974E+04, 0.26120E+04, 0.26267E+04, 0.26415E+04, 0.26563E+04, 0.26711E+04 &
, 0.26861E+04, 0.27011E+04, 0.27162E+04, 0.27313E+04, 0.27466E+04, 0.27618E+04, 0.27772E+04, 0.27926E+04, 0.28081E+04, 0.28237E+04 &
, 0.28394E+04, 0.28551E+04, 0.28708E+04, 0.28867E+04, 0.29026E+04, 0.29186E+04, 0.29347E+04, 0.29508E+04, 0.29671E+04, 0.29834E+04 &
, 0.29997E+04, 0.30162E+04, 0.30327E+04, 0.30493E+04, 0.30659E+04, 0.30826E+04, 0.30995E+04, 0.31163E+04, 0.31333E+04, 0.31503E+04 &
, 0.31674E+04, 0.31846E+04, 0.32019E+04, 0.32192E+04, 0.32367E+04, 0.32541E+04, 0.32717E+04, 0.32894E+04, 0.33071E+04, 0.33249E+04 &
, 0.33428E+04, 0.33608E+04, 0.33788E+04, 0.33969E+04, 0.34151E+04, 0.34334E+04, 0.34518E+04, 0.34702E+04, 0.34887E+04, 0.35073E+04 /
 !            2           8
 DATA(QofT(          17 ,J),J=1,501)/   235.690002441406       &
, 0.25905E+03, 0.28242E+03, 0.30579E+03, 0.32916E+03, 0.35252E+03, 0.37589E+03, 0.39926E+03, 0.42263E+03, 0.44600E+03, 0.46937E+03 &
, 0.49274E+03, 0.51610E+03, 0.53947E+03, 0.56284E+03, 0.58621E+03, 0.60958E+03, 0.63295E+03, 0.65632E+03, 0.67969E+03, 0.70306E+03 &
, 0.72643E+03, 0.74980E+03, 0.77317E+03, 0.79655E+03, 0.81992E+03, 0.84329E+03, 0.86666E+03, 0.89004E+03, 0.91341E+03, 0.93678E+03 &
, 0.96016E+03, 0.98354E+03, 0.10069E+04, 0.10303E+04, 0.10537E+04, 0.10771E+04, 0.11005E+04, 0.11239E+04, 0.11473E+04, 0.11707E+04 &
, 0.11941E+04, 0.12175E+04, 0.12409E+04, 0.12644E+04, 0.12878E+04, 0.13113E+04, 0.13347E+04, 0.13582E+04, 0.13817E+04, 0.14052E+04 &
, 0.14288E+04, 0.14523E+04, 0.14759E+04, 0.14995E+04, 0.15231E+04, 0.15468E+04, 0.15705E+04, 0.15942E+04, 0.16179E+04, 0.16417E+04 &
, 0.16655E+04, 0.16893E+04, 0.17132E+04, 0.17371E+04, 0.17611E+04, 0.17851E+04, 0.18092E+04, 0.18333E+04, 0.18574E+04, 0.18817E+04 &
, 0.19059E+04, 0.19303E+04, 0.19546E+04, 0.19791E+04, 0.20036E+04, 0.20282E+04, 0.20528E+04, 0.20776E+04, 0.21023E+04, 0.21272E+04 &
, 0.21522E+04, 0.21772E+04, 0.22023E+04, 0.22275E+04, 0.22528E+04, 0.22781E+04, 0.23036E+04, 0.23291E+04, 0.23548E+04, 0.23805E+04 &
, 0.24063E+04, 0.24323E+04, 0.24583E+04, 0.24845E+04, 0.25107E+04, 0.25371E+04, 0.25635E+04, 0.25901E+04, 0.26168E+04, 0.26436E+04 &
, 0.26706E+04, 0.26976E+04, 0.27248E+04, 0.27521E+04, 0.27795E+04, 0.28070E+04, 0.28347E+04, 0.28625E+04, 0.28905E+04, 0.29186E+04 &
, 0.29468E+04, 0.29751E+04, 0.30036E+04, 0.30323E+04, 0.30610E+04, 0.30900E+04, 0.31190E+04, 0.31483E+04, 0.31776E+04, 0.32072E+04 &
, 0.32368E+04, 0.32667E+04, 0.32967E+04, 0.33268E+04, 0.33571E+04, 0.33876E+04, 0.34183E+04, 0.34491E+04, 0.34800E+04, 0.35112E+04 &
, 0.35425E+04, 0.35740E+04, 0.36056E+04, 0.36375E+04, 0.36695E+04, 0.37017E+04, 0.37340E+04, 0.37666E+04, 0.37993E+04, 0.38322E+04 &
, 0.38653E+04, 0.38986E+04, 0.39321E+04, 0.39658E+04, 0.39996E+04, 0.40337E+04, 0.40679E+04, 0.41024E+04, 0.41370E+04, 0.41719E+04 &
, 0.42069E+04, 0.42421E+04, 0.42776E+04, 0.43132E+04, 0.43491E+04, 0.43852E+04, 0.44215E+04, 0.44579E+04, 0.44946E+04, 0.45316E+04 &
, 0.45687E+04, 0.46060E+04, 0.46436E+04, 0.46814E+04, 0.47194E+04, 0.47576E+04, 0.47961E+04, 0.48348E+04, 0.48737E+04, 0.49128E+04 &
, 0.49522E+04, 0.49917E+04, 0.50316E+04, 0.50716E+04, 0.51119E+04, 0.51525E+04, 0.51932E+04, 0.52342E+04, 0.52755E+04, 0.53170E+04 &
, 0.53587E+04, 0.54007E+04, 0.54429E+04, 0.54854E+04, 0.55281E+04, 0.55711E+04, 0.56144E+04, 0.56578E+04, 0.57016E+04, 0.57456E+04 &
, 0.57898E+04, 0.58344E+04, 0.58791E+04, 0.59242E+04, 0.59695E+04, 0.60150E+04, 0.60609E+04, 0.61070E+04, 0.61533E+04, 0.62000E+04 &
, 0.62469E+04, 0.62941E+04, 0.63416E+04, 0.63893E+04, 0.64373E+04, 0.64856E+04, 0.65342E+04, 0.65831E+04, 0.66322E+04, 0.66817E+04 &
, 0.67314E+04, 0.67814E+04, 0.68317E+04, 0.68823E+04, 0.69332E+04, 0.69844E+04, 0.70359E+04, 0.70876E+04, 0.71397E+04, 0.71921E+04 &
, 0.72448E+04, 0.72978E+04, 0.73510E+04, 0.74046E+04, 0.74585E+04, 0.75128E+04, 0.75673E+04, 0.76221E+04, 0.76773E+04, 0.77327E+04 &
, 0.77885E+04, 0.78446E+04, 0.79011E+04, 0.79578E+04, 0.80149E+04, 0.80723E+04, 0.81300E+04, 0.81881E+04, 0.82464E+04, 0.83052E+04 &
, 0.83642E+04, 0.84236E+04, 0.84833E+04, 0.85434E+04, 0.86038E+04, 0.86645E+04, 0.87256E+04, 0.87870E+04, 0.88488E+04, 0.89109E+04 &
, 0.89734E+04, 0.90362E+04, 0.90994E+04, 0.91629E+04, 0.92268E+04, 0.92910E+04, 0.93556E+04, 0.94206E+04, 0.94859E+04, 0.95516E+04 &
, 0.96176E+04, 0.96840E+04, 0.97508E+04, 0.98180E+04, 0.98855E+04, 0.99534E+04, 0.10022E+05, 0.10090E+05, 0.10159E+05, 0.10229E+05 &
, 0.10299E+05, 0.10369E+05, 0.10439E+05, 0.10510E+05, 0.10582E+05, 0.10653E+05, 0.10726E+05, 0.10798E+05, 0.10871E+05, 0.10944E+05 &
, 0.11018E+05, 0.11092E+05, 0.11167E+05, 0.11242E+05, 0.11317E+05, 0.11393E+05, 0.11469E+05, 0.11545E+05, 0.11622E+05, 0.11700E+05 &
, 0.11778E+05, 0.11856E+05, 0.11935E+05, 0.12014E+05, 0.12093E+05, 0.12173E+05, 0.12254E+05, 0.12334E+05, 0.12416E+05, 0.12497E+05 &
, 0.12579E+05, 0.12662E+05, 0.12745E+05, 0.12828E+05, 0.12912E+05, 0.12996E+05, 0.13081E+05, 0.13166E+05, 0.13252E+05, 0.13338E+05 &
, 0.13425E+05, 0.13512E+05, 0.13599E+05, 0.13687E+05, 0.13776E+05, 0.13865E+05, 0.13954E+05, 0.14044E+05, 0.14134E+05, 0.14225E+05 &
, 0.14316E+05, 0.14408E+05, 0.14500E+05, 0.14593E+05, 0.14686E+05, 0.14779E+05, 0.14874E+05, 0.14968E+05, 0.15063E+05, 0.15159E+05 &
, 0.15255E+05, 0.15352E+05, 0.15449E+05, 0.15546E+05, 0.15644E+05, 0.15743E+05, 0.15842E+05, 0.15942E+05, 0.16042E+05, 0.16142E+05 &
, 0.16244E+05, 0.16345E+05, 0.16447E+05, 0.16550E+05, 0.16653E+05, 0.16757E+05, 0.16862E+05, 0.16966E+05, 0.17072E+05, 0.17178E+05 &
, 0.17284E+05, 0.17391E+05, 0.17499E+05, 0.17607E+05, 0.17715E+05, 0.17824E+05, 0.17934E+05, 0.18044E+05, 0.18155E+05, 0.18266E+05 &
, 0.18378E+05, 0.18491E+05, 0.18604E+05, 0.18718E+05, 0.18832E+05, 0.18946E+05, 0.19062E+05, 0.19178E+05, 0.19294E+05, 0.19411E+05 &
, 0.19529E+05, 0.19647E+05, 0.19766E+05, 0.19885E+05, 0.20005E+05, 0.20126E+05, 0.20247E+05, 0.20369E+05, 0.20491E+05, 0.20614E+05 &
, 0.20738E+05, 0.20862E+05, 0.20987E+05, 0.21112E+05, 0.21238E+05, 0.21365E+05, 0.21492E+05, 0.21620E+05, 0.21748E+05, 0.21877E+05 &
, 0.22007E+05, 0.22138E+05, 0.22269E+05, 0.22400E+05, 0.22532E+05, 0.22665E+05, 0.22799E+05, 0.22933E+05, 0.23068E+05, 0.23204E+05 &
, 0.23340E+05, 0.23476E+05, 0.23614E+05, 0.23752E+05, 0.23891E+05, 0.24030E+05, 0.24170E+05, 0.24311E+05, 0.24453E+05, 0.24595E+05 &
, 0.24738E+05, 0.24881E+05, 0.25025E+05, 0.25170E+05, 0.25316E+05, 0.25462E+05, 0.25609E+05, 0.25757E+05, 0.25905E+05, 0.26054E+05 &
, 0.26204E+05, 0.26354E+05, 0.26505E+05, 0.26657E+05, 0.26810E+05, 0.26963E+05, 0.27117E+05, 0.27272E+05, 0.27427E+05, 0.27583E+05 &
, 0.27740E+05, 0.27898E+05, 0.28056E+05, 0.28215E+05, 0.28375E+05, 0.28536E+05, 0.28697E+05, 0.28859E+05, 0.29022E+05, 0.29186E+05 &
, 0.29350E+05, 0.29515E+05, 0.29681E+05, 0.29847E+05, 0.30015E+05, 0.30183E+05, 0.30352E+05, 0.30522E+05, 0.30692E+05, 0.30863E+05 &
, 0.31035E+05, 0.31208E+05, 0.31382E+05, 0.31556E+05, 0.31732E+05, 0.31908E+05, 0.32084E+05, 0.32262E+05, 0.32440E+05, 0.32620E+05 &
, 0.32800E+05, 0.32980E+05, 0.33162E+05, 0.33345E+05, 0.33528E+05, 0.33712E+05, 0.33897E+05, 0.34083E+05, 0.34270E+05, 0.34457E+05 &
, 0.34645E+05, 0.34834E+05, 0.35025E+05, 0.35215E+05, 0.35407E+05, 0.35600E+05, 0.35793E+05, 0.35987E+05, 0.36183E+05, 0.36379E+05 &
, 0.36575E+05, 0.36773E+05, 0.36972E+05, 0.37171E+05, 0.37372E+05, 0.37573E+05, 0.37775E+05, 0.37978E+05, 0.38182E+05, 0.38387E+05 &
, 0.38593E+05, 0.38800E+05, 0.39007E+05, 0.39216E+05, 0.39425E+05, 0.39635E+05, 0.39847E+05, 0.40059E+05, 0.40272E+05, 0.40486E+05 /
 !            2           9
 DATA(QofT(          18 ,J),J=1,501)/   687.469970703125       &
, 0.75561E+03, 0.82375E+03, 0.89189E+03, 0.96004E+03, 0.10282E+04, 0.10963E+04, 0.11645E+04, 0.12326E+04, 0.13008E+04, 0.13689E+04 &
, 0.14371E+04, 0.15052E+04, 0.15734E+04, 0.16415E+04, 0.17097E+04, 0.17778E+04, 0.18460E+04, 0.19141E+04, 0.19823E+04, 0.20504E+04 &
, 0.21186E+04, 0.21867E+04, 0.22549E+04, 0.23230E+04, 0.23912E+04, 0.24593E+04, 0.25275E+04, 0.25956E+04, 0.26638E+04, 0.27320E+04 &
, 0.28001E+04, 0.28683E+04, 0.29365E+04, 0.30047E+04, 0.30729E+04, 0.31411E+04, 0.32093E+04, 0.32775E+04, 0.33457E+04, 0.34140E+04 &
, 0.34823E+04, 0.35506E+04, 0.36189E+04, 0.36872E+04, 0.37556E+04, 0.38240E+04, 0.38924E+04, 0.39609E+04, 0.40294E+04, 0.40980E+04 &
, 0.41666E+04, 0.42353E+04, 0.43040E+04, 0.43728E+04, 0.44417E+04, 0.45107E+04, 0.45797E+04, 0.46488E+04, 0.47180E+04, 0.47873E+04 &
, 0.48567E+04, 0.49262E+04, 0.49958E+04, 0.50655E+04, 0.51354E+04, 0.52054E+04, 0.52755E+04, 0.53458E+04, 0.54162E+04, 0.54867E+04 &
, 0.55575E+04, 0.56284E+04, 0.56994E+04, 0.57707E+04, 0.58421E+04, 0.59137E+04, 0.59855E+04, 0.60576E+04, 0.61298E+04, 0.62022E+04 &
, 0.62749E+04, 0.63478E+04, 0.64209E+04, 0.64943E+04, 0.65679E+04, 0.66418E+04, 0.67159E+04, 0.67903E+04, 0.68650E+04, 0.69399E+04 &
, 0.70152E+04, 0.70907E+04, 0.71665E+04, 0.72426E+04, 0.73190E+04, 0.73958E+04, 0.74728E+04, 0.75502E+04, 0.76279E+04, 0.77060E+04 &
, 0.77843E+04, 0.78631E+04, 0.79421E+04, 0.80216E+04, 0.81014E+04, 0.81816E+04, 0.82621E+04, 0.83430E+04, 0.84243E+04, 0.85060E+04 &
, 0.85881E+04, 0.86706E+04, 0.87535E+04, 0.88368E+04, 0.89205E+04, 0.90046E+04, 0.90892E+04, 0.91742E+04, 0.92596E+04, 0.93455E+04 &
, 0.94318E+04, 0.95185E+04, 0.96058E+04, 0.96934E+04, 0.97816E+04, 0.98702E+04, 0.99592E+04, 0.10049E+05, 0.10139E+05, 0.10229E+05 &
, 0.10320E+05, 0.10412E+05, 0.10504E+05, 0.10596E+05, 0.10689E+05, 0.10783E+05, 0.10877E+05, 0.10972E+05, 0.11067E+05, 0.11162E+05 &
, 0.11258E+05, 0.11355E+05, 0.11452E+05, 0.11550E+05, 0.11649E+05, 0.11748E+05, 0.11847E+05, 0.11947E+05, 0.12048E+05, 0.12149E+05 &
, 0.12251E+05, 0.12353E+05, 0.12456E+05, 0.12559E+05, 0.12664E+05, 0.12768E+05, 0.12874E+05, 0.12979E+05, 0.13086E+05, 0.13193E+05 &
, 0.13301E+05, 0.13409E+05, 0.13518E+05, 0.13628E+05, 0.13738E+05, 0.13849E+05, 0.13961E+05, 0.14073E+05, 0.14186E+05, 0.14300E+05 &
, 0.14414E+05, 0.14529E+05, 0.14644E+05, 0.14760E+05, 0.14877E+05, 0.14995E+05, 0.15113E+05, 0.15232E+05, 0.15352E+05, 0.15472E+05 &
, 0.15593E+05, 0.15715E+05, 0.15837E+05, 0.15961E+05, 0.16085E+05, 0.16209E+05, 0.16335E+05, 0.16461E+05, 0.16587E+05, 0.16715E+05 &
, 0.16843E+05, 0.16972E+05, 0.17102E+05, 0.17233E+05, 0.17364E+05, 0.17496E+05, 0.17629E+05, 0.17763E+05, 0.17897E+05, 0.18032E+05 &
, 0.18168E+05, 0.18305E+05, 0.18442E+05, 0.18581E+05, 0.18720E+05, 0.18860E+05, 0.19001E+05, 0.19142E+05, 0.19285E+05, 0.19428E+05 &
, 0.19572E+05, 0.19717E+05, 0.19863E+05, 0.20009E+05, 0.20156E+05, 0.20305E+05, 0.20454E+05, 0.20604E+05, 0.20755E+05, 0.20906E+05 &
, 0.21059E+05, 0.21212E+05, 0.21367E+05, 0.21522E+05, 0.21678E+05, 0.21835E+05, 0.21993E+05, 0.22151E+05, 0.22311E+05, 0.22472E+05 &
, 0.22633E+05, 0.22796E+05, 0.22959E+05, 0.23123E+05, 0.23288E+05, 0.23455E+05, 0.23622E+05, 0.23790E+05, 0.23959E+05, 0.24129E+05 &
, 0.24299E+05, 0.24471E+05, 0.24644E+05, 0.24818E+05, 0.24993E+05, 0.25168E+05, 0.25345E+05, 0.25523E+05, 0.25702E+05, 0.25881E+05 &
, 0.26062E+05, 0.26244E+05, 0.26427E+05, 0.26610E+05, 0.26795E+05, 0.26981E+05, 0.27168E+05, 0.27356E+05, 0.27545E+05, 0.27735E+05 &
, 0.27926E+05, 0.28118E+05, 0.28311E+05, 0.28505E+05, 0.28700E+05, 0.28897E+05, 0.29094E+05, 0.29293E+05, 0.29492E+05, 0.29693E+05 &
, 0.29895E+05, 0.30098E+05, 0.30302E+05, 0.30507E+05, 0.30713E+05, 0.30921E+05, 0.31129E+05, 0.31339E+05, 0.31550E+05, 0.31762E+05 &
, 0.31975E+05, 0.32189E+05, 0.32404E+05, 0.32621E+05, 0.32839E+05, 0.33058E+05, 0.33278E+05, 0.33499E+05, 0.33722E+05, 0.33945E+05 &
, 0.34170E+05, 0.34396E+05, 0.34624E+05, 0.34852E+05, 0.35082E+05, 0.35313E+05, 0.35545E+05, 0.35778E+05, 0.36013E+05, 0.36249E+05 &
, 0.36486E+05, 0.36725E+05, 0.36964E+05, 0.37205E+05, 0.37448E+05, 0.37691E+05, 0.37936E+05, 0.38182E+05, 0.38430E+05, 0.38678E+05 &
, 0.38928E+05, 0.39180E+05, 0.39432E+05, 0.39686E+05, 0.39942E+05, 0.40198E+05, 0.40457E+05, 0.40716E+05, 0.40977E+05, 0.41239E+05 &
, 0.41502E+05, 0.41767E+05, 0.42033E+05, 0.42301E+05, 0.42570E+05, 0.42840E+05, 0.43112E+05, 0.43385E+05, 0.43659E+05, 0.43935E+05 &
, 0.44213E+05, 0.44492E+05, 0.44772E+05, 0.45054E+05, 0.45337E+05, 0.45621E+05, 0.45908E+05, 0.46195E+05, 0.46484E+05, 0.46775E+05 &
, 0.47067E+05, 0.47360E+05, 0.47655E+05, 0.47951E+05, 0.48249E+05, 0.48549E+05, 0.48850E+05, 0.49152E+05, 0.49456E+05, 0.49762E+05 &
, 0.50069E+05, 0.50378E+05, 0.50688E+05, 0.51000E+05, 0.51313E+05, 0.51628E+05, 0.51944E+05, 0.52262E+05, 0.52582E+05, 0.52903E+05 &
, 0.53226E+05, 0.53551E+05, 0.53877E+05, 0.54205E+05, 0.54534E+05, 0.54865E+05, 0.55198E+05, 0.55532E+05, 0.55868E+05, 0.56205E+05 &
, 0.56545E+05, 0.56886E+05, 0.57228E+05, 0.57573E+05, 0.57919E+05, 0.58266E+05, 0.58616E+05, 0.58967E+05, 0.59320E+05, 0.59675E+05 &
, 0.60031E+05, 0.60389E+05, 0.60749E+05, 0.61110E+05, 0.61474E+05, 0.61839E+05, 0.62206E+05, 0.62574E+05, 0.62945E+05, 0.63317E+05 &
, 0.63691E+05, 0.64067E+05, 0.64445E+05, 0.64824E+05, 0.65206E+05, 0.65589E+05, 0.65974E+05, 0.66361E+05, 0.66749E+05, 0.67140E+05 &
, 0.67533E+05, 0.67927E+05, 0.68323E+05, 0.68721E+05, 0.69121E+05, 0.69523E+05, 0.69927E+05, 0.70333E+05, 0.70741E+05, 0.71150E+05 &
, 0.71562E+05, 0.71975E+05, 0.72391E+05, 0.72808E+05, 0.73228E+05, 0.73649E+05, 0.74073E+05, 0.74498E+05, 0.74925E+05, 0.75355E+05 &
, 0.75786E+05, 0.76220E+05, 0.76655E+05, 0.77093E+05, 0.77532E+05, 0.77974E+05, 0.78418E+05, 0.78863E+05, 0.79311E+05, 0.79761E+05 &
, 0.80213E+05, 0.80667E+05, 0.81124E+05, 0.81582E+05, 0.82042E+05, 0.82505E+05, 0.82970E+05, 0.83437E+05, 0.83906E+05, 0.84377E+05 &
, 0.84850E+05, 0.85326E+05, 0.85803E+05, 0.86283E+05, 0.86766E+05, 0.87250E+05, 0.87736E+05, 0.88225E+05, 0.88716E+05, 0.89209E+05 &
, 0.89705E+05, 0.90203E+05, 0.90703E+05, 0.91205E+05, 0.91710E+05, 0.92216E+05, 0.92726E+05, 0.93237E+05, 0.93751E+05, 0.94267E+05 &
, 0.94785E+05, 0.95306E+05, 0.95829E+05, 0.96355E+05, 0.96882E+05, 0.97413E+05, 0.97945E+05, 0.98480E+05, 0.99018E+05, 0.99557E+05 &
, 0.10010E+06, 0.10064E+06, 0.10119E+06, 0.10174E+06, 0.10229E+06, 0.10285E+06, 0.10340E+06, 0.10396E+06, 0.10452E+06, 0.10509E+06 &
, 0.10566E+06, 0.10622E+06, 0.10680E+06, 0.10737E+06, 0.10795E+06, 0.10853E+06, 0.10911E+06, 0.10969E+06, 0.11028E+06, 0.11087E+06 &
, 0.11146E+06, 0.11206E+06, 0.11265E+06, 0.11325E+06, 0.11386E+06, 0.11446E+06, 0.11507E+06, 0.11568E+06, 0.11629E+06, 0.11691E+06 /
 !            2          10
 DATA(QofT(          19 ,J),J=1,501)/   40.4160003662109       &
, 0.44424E+02, 0.48432E+02, 0.52440E+02, 0.56449E+02, 0.60457E+02, 0.64465E+02, 0.68473E+02, 0.72482E+02, 0.76490E+02, 0.80498E+02 &
, 0.84507E+02, 0.88515E+02, 0.92523E+02, 0.96532E+02, 0.10054E+03, 0.10455E+03, 0.10856E+03, 0.11257E+03, 0.11657E+03, 0.12058E+03 &
, 0.12459E+03, 0.12860E+03, 0.13261E+03, 0.13662E+03, 0.14063E+03, 0.14464E+03, 0.14864E+03, 0.15265E+03, 0.15666E+03, 0.16067E+03 &
, 0.16468E+03, 0.16869E+03, 0.17270E+03, 0.17672E+03, 0.18073E+03, 0.18474E+03, 0.18875E+03, 0.19277E+03, 0.19678E+03, 0.20080E+03 &
, 0.20482E+03, 0.20884E+03, 0.21286E+03, 0.21688E+03, 0.22091E+03, 0.22494E+03, 0.22897E+03, 0.23300E+03, 0.23704E+03, 0.24108E+03 &
, 0.24512E+03, 0.24917E+03, 0.25322E+03, 0.25728E+03, 0.26134E+03, 0.26541E+03, 0.26948E+03, 0.27356E+03, 0.27765E+03, 0.28174E+03 &
, 0.28584E+03, 0.28994E+03, 0.29406E+03, 0.29818E+03, 0.30231E+03, 0.30645E+03, 0.31060E+03, 0.31476E+03, 0.31893E+03, 0.32311E+03 &
, 0.32730E+03, 0.33150E+03, 0.33571E+03, 0.33994E+03, 0.34418E+03, 0.34843E+03, 0.35269E+03, 0.35697E+03, 0.36126E+03, 0.36557E+03 &
, 0.36989E+03, 0.37423E+03, 0.37858E+03, 0.38295E+03, 0.38733E+03, 0.39173E+03, 0.39615E+03, 0.40059E+03, 0.40504E+03, 0.40951E+03 &
, 0.41400E+03, 0.41851E+03, 0.42304E+03, 0.42759E+03, 0.43216E+03, 0.43675E+03, 0.44136E+03, 0.44599E+03, 0.45064E+03, 0.45532E+03 &
, 0.46001E+03, 0.46473E+03, 0.46947E+03, 0.47424E+03, 0.47903E+03, 0.48384E+03, 0.48867E+03, 0.49354E+03, 0.49842E+03, 0.50333E+03 &
, 0.50827E+03, 0.51323E+03, 0.51822E+03, 0.52323E+03, 0.52827E+03, 0.53334E+03, 0.53843E+03, 0.54356E+03, 0.54871E+03, 0.55388E+03 &
, 0.55909E+03, 0.56433E+03, 0.56959E+03, 0.57489E+03, 0.58021E+03, 0.58556E+03, 0.59095E+03, 0.59636E+03, 0.60181E+03, 0.60728E+03 &
, 0.61279E+03, 0.61833E+03, 0.62390E+03, 0.62950E+03, 0.63514E+03, 0.64081E+03, 0.64651E+03, 0.65224E+03, 0.65801E+03, 0.66381E+03 &
, 0.66965E+03, 0.67552E+03, 0.68142E+03, 0.68736E+03, 0.69334E+03, 0.69935E+03, 0.70539E+03, 0.71147E+03, 0.71759E+03, 0.72375E+03 &
, 0.72994E+03, 0.73617E+03, 0.74243E+03, 0.74873E+03, 0.75507E+03, 0.76145E+03, 0.76787E+03, 0.77432E+03, 0.78082E+03, 0.78735E+03 &
, 0.79392E+03, 0.80054E+03, 0.80719E+03, 0.81388E+03, 0.82061E+03, 0.82739E+03, 0.83420E+03, 0.84106E+03, 0.84795E+03, 0.85489E+03 &
, 0.86187E+03, 0.86889E+03, 0.87596E+03, 0.88306E+03, 0.89021E+03, 0.89741E+03, 0.90464E+03, 0.91192E+03, 0.91925E+03, 0.92662E+03 &
, 0.93403E+03, 0.94149E+03, 0.94899E+03, 0.95654E+03, 0.96413E+03, 0.97177E+03, 0.97946E+03, 0.98719E+03, 0.99497E+03, 0.10028E+04 &
, 0.10107E+04, 0.10186E+04, 0.10266E+04, 0.10346E+04, 0.10426E+04, 0.10507E+04, 0.10589E+04, 0.10671E+04, 0.10754E+04, 0.10837E+04 &
, 0.10920E+04, 0.11004E+04, 0.11089E+04, 0.11174E+04, 0.11260E+04, 0.11346E+04, 0.11432E+04, 0.11520E+04, 0.11607E+04, 0.11695E+04 &
, 0.11784E+04, 0.11873E+04, 0.11963E+04, 0.12053E+04, 0.12144E+04, 0.12235E+04, 0.12327E+04, 0.12420E+04, 0.12513E+04, 0.12606E+04 &
, 0.12700E+04, 0.12795E+04, 0.12890E+04, 0.12986E+04, 0.13082E+04, 0.13179E+04, 0.13277E+04, 0.13375E+04, 0.13473E+04, 0.13572E+04 &
, 0.13672E+04, 0.13773E+04, 0.13873E+04, 0.13975E+04, 0.14077E+04, 0.14180E+04, 0.14283E+04, 0.14387E+04, 0.14491E+04, 0.14597E+04 &
, 0.14702E+04, 0.14809E+04, 0.14915E+04, 0.15023E+04, 0.15131E+04, 0.15240E+04, 0.15349E+04, 0.15459E+04, 0.15570E+04, 0.15681E+04 &
, 0.15793E+04, 0.15906E+04, 0.16019E+04, 0.16133E+04, 0.16248E+04, 0.16363E+04, 0.16479E+04, 0.16595E+04, 0.16712E+04, 0.16830E+04 &
, 0.16949E+04, 0.17068E+04, 0.17188E+04, 0.17308E+04, 0.17430E+04, 0.17551E+04, 0.17674E+04, 0.17797E+04, 0.17921E+04, 0.18046E+04 &
, 0.18171E+04, 0.18297E+04, 0.18424E+04, 0.18552E+04, 0.18680E+04, 0.18809E+04, 0.18938E+04, 0.19069E+04, 0.19200E+04, 0.19332E+04 &
, 0.19464E+04, 0.19598E+04, 0.19732E+04, 0.19866E+04, 0.20002E+04, 0.20138E+04, 0.20275E+04, 0.20413E+04, 0.20551E+04, 0.20691E+04 &
, 0.20831E+04, 0.20972E+04, 0.21113E+04, 0.21256E+04, 0.21399E+04, 0.21543E+04, 0.21687E+04, 0.21833E+04, 0.21979E+04, 0.22126E+04 &
, 0.22274E+04, 0.22423E+04, 0.22572E+04, 0.22723E+04, 0.22874E+04, 0.23026E+04, 0.23179E+04, 0.23332E+04, 0.23487E+04, 0.23642E+04 &
, 0.23798E+04, 0.23955E+04, 0.24113E+04, 0.24271E+04, 0.24431E+04, 0.24591E+04, 0.24752E+04, 0.24914E+04, 0.25077E+04, 0.25241E+04 &
, 0.25406E+04, 0.25571E+04, 0.25738E+04, 0.25905E+04, 0.26073E+04, 0.26242E+04, 0.26412E+04, 0.26583E+04, 0.26755E+04, 0.26928E+04 &
, 0.27101E+04, 0.27276E+04, 0.27451E+04, 0.27627E+04, 0.27805E+04, 0.27983E+04, 0.28162E+04, 0.28342E+04, 0.28523E+04, 0.28705E+04 &
, 0.28888E+04, 0.29072E+04, 0.29257E+04, 0.29442E+04, 0.29629E+04, 0.29817E+04, 0.30005E+04, 0.30195E+04, 0.30386E+04, 0.30577E+04 &
, 0.30770E+04, 0.30963E+04, 0.31158E+04, 0.31354E+04, 0.31550E+04, 0.31748E+04, 0.31946E+04, 0.32146E+04, 0.32347E+04, 0.32548E+04 &
, 0.32751E+04, 0.32955E+04, 0.33159E+04, 0.33365E+04, 0.33572E+04, 0.33780E+04, 0.33989E+04, 0.34199E+04, 0.34410E+04, 0.34622E+04 &
, 0.34835E+04, 0.35050E+04, 0.35265E+04, 0.35481E+04, 0.35699E+04, 0.35918E+04, 0.36137E+04, 0.36358E+04, 0.36580E+04, 0.36803E+04 &
, 0.37027E+04, 0.37252E+04, 0.37479E+04, 0.37706E+04, 0.37935E+04, 0.38165E+04, 0.38396E+04, 0.38628E+04, 0.38861E+04, 0.39095E+04 &
, 0.39331E+04, 0.39568E+04, 0.39805E+04, 0.40044E+04, 0.40285E+04, 0.40526E+04, 0.40769E+04, 0.41012E+04, 0.41257E+04, 0.41503E+04 &
, 0.41751E+04, 0.41999E+04, 0.42249E+04, 0.42500E+04, 0.42752E+04, 0.43006E+04, 0.43260E+04, 0.43516E+04, 0.43773E+04, 0.44032E+04 &
, 0.44291E+04, 0.44552E+04, 0.44814E+04, 0.45078E+04, 0.45343E+04, 0.45609E+04, 0.45876E+04, 0.46144E+04, 0.46414E+04, 0.46685E+04 &
, 0.46958E+04, 0.47231E+04, 0.47506E+04, 0.47783E+04, 0.48060E+04, 0.48339E+04, 0.48620E+04, 0.48901E+04, 0.49184E+04, 0.49468E+04 &
, 0.49754E+04, 0.50041E+04, 0.50329E+04, 0.50619E+04, 0.50910E+04, 0.51203E+04, 0.51497E+04, 0.51792E+04, 0.52088E+04, 0.52386E+04 &
, 0.52686E+04, 0.52987E+04, 0.53289E+04, 0.53592E+04, 0.53897E+04, 0.54204E+04, 0.54512E+04, 0.54821E+04, 0.55132E+04, 0.55444E+04 &
, 0.55758E+04, 0.56073E+04, 0.56389E+04, 0.56707E+04, 0.57027E+04, 0.57348E+04, 0.57670E+04, 0.57994E+04, 0.58320E+04, 0.58647E+04 &
, 0.58975E+04, 0.59305E+04, 0.59637E+04, 0.59970E+04, 0.60304E+04, 0.60640E+04, 0.60978E+04, 0.61317E+04, 0.61657E+04, 0.62000E+04 &
, 0.62343E+04, 0.62689E+04, 0.63036E+04, 0.63384E+04, 0.63734E+04, 0.64086E+04, 0.64439E+04, 0.64794E+04, 0.65150E+04, 0.65508E+04 &
, 0.65868E+04, 0.66229E+04, 0.66592E+04, 0.66956E+04, 0.67323E+04, 0.67690E+04, 0.68060E+04, 0.68431E+04, 0.68804E+04, 0.69178E+04 &
, 0.69554E+04, 0.69932E+04, 0.70311E+04, 0.70692E+04, 0.71075E+04, 0.71460E+04, 0.71846E+04, 0.72234E+04, 0.72623E+04, 0.73015E+04 /
 !            2          11
 DATA(QofT(          20 ,J),J=1,501)/   471.359985351562       &
, 0.51809E+03, 0.56482E+03, 0.61156E+03, 0.65829E+03, 0.70502E+03, 0.75176E+03, 0.79850E+03, 0.84523E+03, 0.89197E+03, 0.93870E+03 &
, 0.98544E+03, 0.10322E+04, 0.10789E+04, 0.11257E+04, 0.11724E+04, 0.12191E+04, 0.12659E+04, 0.13126E+04, 0.13593E+04, 0.14061E+04 &
, 0.14528E+04, 0.14996E+04, 0.15463E+04, 0.15930E+04, 0.16398E+04, 0.16865E+04, 0.17333E+04, 0.17800E+04, 0.18268E+04, 0.18735E+04 &
, 0.19203E+04, 0.19670E+04, 0.20138E+04, 0.20606E+04, 0.21073E+04, 0.21541E+04, 0.22009E+04, 0.22477E+04, 0.22945E+04, 0.23414E+04 &
, 0.23882E+04, 0.24351E+04, 0.24820E+04, 0.25289E+04, 0.25758E+04, 0.26228E+04, 0.26698E+04, 0.27168E+04, 0.27638E+04, 0.28109E+04 &
, 0.28581E+04, 0.29053E+04, 0.29525E+04, 0.29998E+04, 0.30471E+04, 0.30945E+04, 0.31420E+04, 0.31896E+04, 0.32372E+04, 0.32849E+04 &
, 0.33326E+04, 0.33805E+04, 0.34284E+04, 0.34765E+04, 0.35246E+04, 0.35728E+04, 0.36212E+04, 0.36697E+04, 0.37182E+04, 0.37669E+04 &
, 0.38157E+04, 0.38647E+04, 0.39138E+04, 0.39630E+04, 0.40124E+04, 0.40619E+04, 0.41116E+04, 0.41614E+04, 0.42114E+04, 0.42615E+04 &
, 0.43118E+04, 0.43623E+04, 0.44130E+04, 0.44639E+04, 0.45149E+04, 0.45662E+04, 0.46176E+04, 0.46692E+04, 0.47211E+04, 0.47731E+04 &
, 0.48254E+04, 0.48779E+04, 0.49306E+04, 0.49835E+04, 0.50367E+04, 0.50901E+04, 0.51437E+04, 0.51976E+04, 0.52518E+04, 0.53061E+04 &
, 0.53608E+04, 0.54157E+04, 0.54708E+04, 0.55263E+04, 0.55820E+04, 0.56379E+04, 0.56942E+04, 0.57507E+04, 0.58075E+04, 0.58646E+04 &
, 0.59220E+04, 0.59797E+04, 0.60377E+04, 0.60960E+04, 0.61546E+04, 0.62135E+04, 0.62727E+04, 0.63323E+04, 0.63922E+04, 0.64523E+04 &
, 0.65129E+04, 0.65737E+04, 0.66349E+04, 0.66964E+04, 0.67583E+04, 0.68205E+04, 0.68830E+04, 0.69459E+04, 0.70092E+04, 0.70728E+04 &
, 0.71368E+04, 0.72012E+04, 0.72659E+04, 0.73309E+04, 0.73964E+04, 0.74622E+04, 0.75284E+04, 0.75950E+04, 0.76620E+04, 0.77294E+04 &
, 0.77971E+04, 0.78653E+04, 0.79339E+04, 0.80028E+04, 0.80722E+04, 0.81420E+04, 0.82121E+04, 0.82827E+04, 0.83537E+04, 0.84252E+04 &
, 0.84970E+04, 0.85693E+04, 0.86420E+04, 0.87152E+04, 0.87887E+04, 0.88628E+04, 0.89372E+04, 0.90121E+04, 0.90875E+04, 0.91632E+04 &
, 0.92395E+04, 0.93162E+04, 0.93934E+04, 0.94710E+04, 0.95491E+04, 0.96276E+04, 0.97067E+04, 0.97862E+04, 0.98661E+04, 0.99466E+04 &
, 0.10028E+05, 0.10109E+05, 0.10191E+05, 0.10273E+05, 0.10356E+05, 0.10440E+05, 0.10523E+05, 0.10608E+05, 0.10693E+05, 0.10778E+05 &
, 0.10864E+05, 0.10951E+05, 0.11037E+05, 0.11125E+05, 0.11213E+05, 0.11302E+05, 0.11391E+05, 0.11480E+05, 0.11570E+05, 0.11661E+05 &
, 0.11752E+05, 0.11844E+05, 0.11936E+05, 0.12029E+05, 0.12122E+05, 0.12216E+05, 0.12311E+05, 0.12406E+05, 0.12502E+05, 0.12598E+05 &
, 0.12695E+05, 0.12792E+05, 0.12890E+05, 0.12989E+05, 0.13088E+05, 0.13187E+05, 0.13288E+05, 0.13389E+05, 0.13490E+05, 0.13592E+05 &
, 0.13695E+05, 0.13798E+05, 0.13902E+05, 0.14006E+05, 0.14112E+05, 0.14217E+05, 0.14324E+05, 0.14431E+05, 0.14538E+05, 0.14647E+05 &
, 0.14755E+05, 0.14865E+05, 0.14975E+05, 0.15086E+05, 0.15197E+05, 0.15309E+05, 0.15422E+05, 0.15536E+05, 0.15650E+05, 0.15764E+05 &
, 0.15880E+05, 0.15996E+05, 0.16113E+05, 0.16230E+05, 0.16348E+05, 0.16467E+05, 0.16586E+05, 0.16707E+05, 0.16828E+05, 0.16949E+05 &
, 0.17071E+05, 0.17194E+05, 0.17318E+05, 0.17442E+05, 0.17567E+05, 0.17693E+05, 0.17820E+05, 0.17947E+05, 0.18075E+05, 0.18204E+05 &
, 0.18333E+05, 0.18463E+05, 0.18594E+05, 0.18726E+05, 0.18858E+05, 0.18991E+05, 0.19125E+05, 0.19260E+05, 0.19395E+05, 0.19532E+05 &
, 0.19669E+05, 0.19806E+05, 0.19945E+05, 0.20084E+05, 0.20224E+05, 0.20365E+05, 0.20507E+05, 0.20649E+05, 0.20792E+05, 0.20936E+05 &
, 0.21081E+05, 0.21227E+05, 0.21373E+05, 0.21521E+05, 0.21669E+05, 0.21818E+05, 0.21967E+05, 0.22118E+05, 0.22269E+05, 0.22422E+05 &
, 0.22575E+05, 0.22729E+05, 0.22884E+05, 0.23039E+05, 0.23196E+05, 0.23353E+05, 0.23511E+05, 0.23670E+05, 0.23830E+05, 0.23991E+05 &
, 0.24153E+05, 0.24316E+05, 0.24479E+05, 0.24644E+05, 0.24809E+05, 0.24975E+05, 0.25142E+05, 0.25310E+05, 0.25479E+05, 0.25649E+05 &
, 0.25820E+05, 0.25991E+05, 0.26164E+05, 0.26337E+05, 0.26512E+05, 0.26687E+05, 0.26864E+05, 0.27041E+05, 0.27219E+05, 0.27398E+05 &
, 0.27578E+05, 0.27760E+05, 0.27942E+05, 0.28125E+05, 0.28309E+05, 0.28494E+05, 0.28680E+05, 0.28867E+05, 0.29055E+05, 0.29244E+05 &
, 0.29434E+05, 0.29625E+05, 0.29817E+05, 0.30010E+05, 0.30204E+05, 0.30399E+05, 0.30595E+05, 0.30792E+05, 0.30990E+05, 0.31189E+05 &
, 0.31390E+05, 0.31591E+05, 0.31793E+05, 0.31997E+05, 0.32201E+05, 0.32407E+05, 0.32613E+05, 0.32821E+05, 0.33030E+05, 0.33239E+05 &
, 0.33450E+05, 0.33662E+05, 0.33876E+05, 0.34090E+05, 0.34305E+05, 0.34522E+05, 0.34739E+05, 0.34958E+05, 0.35178E+05, 0.35398E+05 &
, 0.35620E+05, 0.35844E+05, 0.36068E+05, 0.36294E+05, 0.36520E+05, 0.36748E+05, 0.36977E+05, 0.37207E+05, 0.37438E+05, 0.37671E+05 &
, 0.37904E+05, 0.38139E+05, 0.38375E+05, 0.38612E+05, 0.38851E+05, 0.39091E+05, 0.39331E+05, 0.39573E+05, 0.39817E+05, 0.40061E+05 &
, 0.40307E+05, 0.40554E+05, 0.40802E+05, 0.41051E+05, 0.41302E+05, 0.41554E+05, 0.41807E+05, 0.42062E+05, 0.42317E+05, 0.42574E+05 &
, 0.42833E+05, 0.43092E+05, 0.43353E+05, 0.43615E+05, 0.43879E+05, 0.44143E+05, 0.44409E+05, 0.44677E+05, 0.44945E+05, 0.45215E+05 &
, 0.45487E+05, 0.45759E+05, 0.46033E+05, 0.46309E+05, 0.46585E+05, 0.46863E+05, 0.47143E+05, 0.47424E+05, 0.47706E+05, 0.47989E+05 &
, 0.48274E+05, 0.48560E+05, 0.48848E+05, 0.49137E+05, 0.49428E+05, 0.49720E+05, 0.50013E+05, 0.50308E+05, 0.50604E+05, 0.50901E+05 &
, 0.51200E+05, 0.51501E+05, 0.51803E+05, 0.52106E+05, 0.52411E+05, 0.52717E+05, 0.53025E+05, 0.53334E+05, 0.53644E+05, 0.53957E+05 &
, 0.54270E+05, 0.54585E+05, 0.54902E+05, 0.55220E+05, 0.55540E+05, 0.55861E+05, 0.56184E+05, 0.56508E+05, 0.56833E+05, 0.57161E+05 &
, 0.57490E+05, 0.57820E+05, 0.58152E+05, 0.58485E+05, 0.58820E+05, 0.59157E+05, 0.59495E+05, 0.59835E+05, 0.60176E+05, 0.60519E+05 &
, 0.60864E+05, 0.61210E+05, 0.61558E+05, 0.61907E+05, 0.62258E+05, 0.62611E+05, 0.62966E+05, 0.63321E+05, 0.63679E+05, 0.64038E+05 &
, 0.64399E+05, 0.64762E+05, 0.65126E+05, 0.65492E+05, 0.65860E+05, 0.66229E+05, 0.66600E+05, 0.66973E+05, 0.67347E+05, 0.67723E+05 &
, 0.68101E+05, 0.68481E+05, 0.68862E+05, 0.69245E+05, 0.69630E+05, 0.70017E+05, 0.70405E+05, 0.70795E+05, 0.71187E+05, 0.71581E+05 &
, 0.71976E+05, 0.72374E+05, 0.72773E+05, 0.73173E+05, 0.73576E+05, 0.73980E+05, 0.74387E+05, 0.74795E+05, 0.75205E+05, 0.75617E+05 &
, 0.76030E+05, 0.76446E+05, 0.76863E+05, 0.77282E+05, 0.77703E+05, 0.78126E+05, 0.78551E+05, 0.78978E+05, 0.79407E+05, 0.79837E+05 &
, 0.80270E+05, 0.80704E+05, 0.81140E+05, 0.81579E+05, 0.82019E+05, 0.82461E+05, 0.82905E+05, 0.83351E+05, 0.83799E+05, 0.84249E+05 /
 !            2          12
 DATA(QofT(          21 ,J),J=1,501)/   1374.90002441406       &
, 0.15112E+04, 0.16474E+04, 0.17837E+04, 0.19200E+04, 0.20563E+04, 0.21926E+04, 0.23289E+04, 0.24652E+04, 0.26015E+04, 0.27377E+04 &
, 0.28740E+04, 0.30103E+04, 0.31466E+04, 0.32829E+04, 0.34192E+04, 0.35555E+04, 0.36918E+04, 0.38281E+04, 0.39644E+04, 0.41007E+04 &
, 0.42370E+04, 0.43733E+04, 0.45096E+04, 0.46459E+04, 0.47822E+04, 0.49185E+04, 0.50548E+04, 0.51911E+04, 0.53275E+04, 0.54638E+04 &
, 0.56001E+04, 0.57365E+04, 0.58729E+04, 0.60092E+04, 0.61457E+04, 0.62821E+04, 0.64185E+04, 0.65550E+04, 0.66915E+04, 0.68281E+04 &
, 0.69647E+04, 0.71013E+04, 0.72380E+04, 0.73748E+04, 0.75117E+04, 0.76486E+04, 0.77856E+04, 0.79227E+04, 0.80599E+04, 0.81972E+04 &
, 0.83347E+04, 0.84722E+04, 0.86100E+04, 0.87478E+04, 0.88859E+04, 0.90241E+04, 0.91625E+04, 0.93010E+04, 0.94398E+04, 0.95788E+04 &
, 0.97181E+04, 0.98576E+04, 0.99973E+04, 0.10137E+05, 0.10278E+05, 0.10418E+05, 0.10559E+05, 0.10700E+05, 0.10842E+05, 0.10984E+05 &
, 0.11126E+05, 0.11269E+05, 0.11412E+05, 0.11555E+05, 0.11699E+05, 0.11843E+05, 0.11988E+05, 0.12133E+05, 0.12279E+05, 0.12425E+05 &
, 0.12571E+05, 0.12718E+05, 0.12866E+05, 0.13014E+05, 0.13163E+05, 0.13312E+05, 0.13462E+05, 0.13612E+05, 0.13763E+05, 0.13914E+05 &
, 0.14067E+05, 0.14219E+05, 0.14373E+05, 0.14527E+05, 0.14682E+05, 0.14837E+05, 0.14993E+05, 0.15150E+05, 0.15308E+05, 0.15466E+05 &
, 0.15625E+05, 0.15785E+05, 0.15945E+05, 0.16106E+05, 0.16268E+05, 0.16431E+05, 0.16595E+05, 0.16759E+05, 0.16925E+05, 0.17091E+05 &
, 0.17258E+05, 0.17425E+05, 0.17594E+05, 0.17764E+05, 0.17934E+05, 0.18105E+05, 0.18278E+05, 0.18451E+05, 0.18625E+05, 0.18800E+05 &
, 0.18976E+05, 0.19153E+05, 0.19330E+05, 0.19509E+05, 0.19689E+05, 0.19870E+05, 0.20052E+05, 0.20235E+05, 0.20418E+05, 0.20603E+05 &
, 0.20789E+05, 0.20976E+05, 0.21164E+05, 0.21353E+05, 0.21543E+05, 0.21735E+05, 0.21927E+05, 0.22120E+05, 0.22315E+05, 0.22511E+05 &
, 0.22708E+05, 0.22905E+05, 0.23105E+05, 0.23305E+05, 0.23506E+05, 0.23709E+05, 0.23913E+05, 0.24118E+05, 0.24324E+05, 0.24531E+05 &
, 0.24740E+05, 0.24950E+05, 0.25161E+05, 0.25373E+05, 0.25587E+05, 0.25801E+05, 0.26018E+05, 0.26235E+05, 0.26454E+05, 0.26674E+05 &
, 0.26895E+05, 0.27117E+05, 0.27341E+05, 0.27567E+05, 0.27793E+05, 0.28021E+05, 0.28250E+05, 0.28481E+05, 0.28713E+05, 0.28946E+05 &
, 0.29181E+05, 0.29417E+05, 0.29655E+05, 0.29894E+05, 0.30134E+05, 0.30376E+05, 0.30619E+05, 0.30864E+05, 0.31110E+05, 0.31358E+05 &
, 0.31607E+05, 0.31858E+05, 0.32110E+05, 0.32363E+05, 0.32619E+05, 0.32875E+05, 0.33133E+05, 0.33393E+05, 0.33654E+05, 0.33917E+05 &
, 0.34181E+05, 0.34447E+05, 0.34715E+05, 0.34984E+05, 0.35255E+05, 0.35527E+05, 0.35801E+05, 0.36076E+05, 0.36353E+05, 0.36632E+05 &
, 0.36913E+05, 0.37195E+05, 0.37478E+05, 0.37764E+05, 0.38051E+05, 0.38340E+05, 0.38630E+05, 0.38922E+05, 0.39216E+05, 0.39512E+05 &
, 0.39809E+05, 0.40109E+05, 0.40409E+05, 0.40712E+05, 0.41016E+05, 0.41323E+05, 0.41631E+05, 0.41940E+05, 0.42252E+05, 0.42565E+05 &
, 0.42881E+05, 0.43198E+05, 0.43517E+05, 0.43837E+05, 0.44160E+05, 0.44484E+05, 0.44811E+05, 0.45139E+05, 0.45469E+05, 0.45801E+05 &
, 0.46135E+05, 0.46471E+05, 0.46809E+05, 0.47149E+05, 0.47491E+05, 0.47834E+05, 0.48180E+05, 0.48528E+05, 0.48877E+05, 0.49229E+05 &
, 0.49583E+05, 0.49938E+05, 0.50296E+05, 0.50656E+05, 0.51018E+05, 0.51382E+05, 0.51748E+05, 0.52116E+05, 0.52486E+05, 0.52858E+05 &
, 0.53232E+05, 0.53609E+05, 0.53987E+05, 0.54368E+05, 0.54751E+05, 0.55136E+05, 0.55523E+05, 0.55913E+05, 0.56304E+05, 0.56698E+05 &
, 0.57094E+05, 0.57492E+05, 0.57893E+05, 0.58295E+05, 0.58700E+05, 0.59108E+05, 0.59517E+05, 0.59929E+05, 0.60343E+05, 0.60759E+05 &
, 0.61178E+05, 0.61599E+05, 0.62022E+05, 0.62448E+05, 0.62876E+05, 0.63306E+05, 0.63739E+05, 0.64174E+05, 0.64612E+05, 0.65052E+05 &
, 0.65494E+05, 0.65939E+05, 0.66386E+05, 0.66836E+05, 0.67288E+05, 0.67743E+05, 0.68200E+05, 0.68660E+05, 0.69122E+05, 0.69586E+05 &
, 0.70054E+05, 0.70523E+05, 0.70996E+05, 0.71470E+05, 0.71948E+05, 0.72428E+05, 0.72910E+05, 0.73396E+05, 0.73883E+05, 0.74374E+05 &
, 0.74867E+05, 0.75363E+05, 0.75861E+05, 0.76362E+05, 0.76866E+05, 0.77373E+05, 0.77882E+05, 0.78394E+05, 0.78908E+05, 0.79426E+05 &
, 0.79946E+05, 0.80469E+05, 0.80994E+05, 0.81523E+05, 0.82054E+05, 0.82588E+05, 0.83125E+05, 0.83665E+05, 0.84208E+05, 0.84753E+05 &
, 0.85302E+05, 0.85853E+05, 0.86407E+05, 0.86964E+05, 0.87524E+05, 0.88087E+05, 0.88653E+05, 0.89222E+05, 0.89794E+05, 0.90368E+05 &
, 0.90946E+05, 0.91527E+05, 0.92111E+05, 0.92698E+05, 0.93288E+05, 0.93881E+05, 0.94477E+05, 0.95076E+05, 0.95678E+05, 0.96284E+05 &
, 0.96892E+05, 0.97504E+05, 0.98119E+05, 0.98737E+05, 0.99358E+05, 0.99982E+05, 0.10061E+06, 0.10124E+06, 0.10187E+06, 0.10251E+06 &
, 0.10315E+06, 0.10380E+06, 0.10444E+06, 0.10509E+06, 0.10575E+06, 0.10640E+06, 0.10706E+06, 0.10773E+06, 0.10839E+06, 0.10907E+06 &
, 0.10974E+06, 0.11042E+06, 0.11110E+06, 0.11178E+06, 0.11247E+06, 0.11316E+06, 0.11385E+06, 0.11455E+06, 0.11525E+06, 0.11596E+06 &
, 0.11667E+06, 0.11738E+06, 0.11809E+06, 0.11881E+06, 0.11954E+06, 0.12026E+06, 0.12099E+06, 0.12172E+06, 0.12246E+06, 0.12320E+06 &
, 0.12395E+06, 0.12470E+06, 0.12545E+06, 0.12620E+06, 0.12696E+06, 0.12773E+06, 0.12849E+06, 0.12926E+06, 0.13004E+06, 0.13081E+06 &
, 0.13160E+06, 0.13238E+06, 0.13317E+06, 0.13397E+06, 0.13476E+06, 0.13556E+06, 0.13637E+06, 0.13718E+06, 0.13799E+06, 0.13881E+06 &
, 0.13963E+06, 0.14045E+06, 0.14128E+06, 0.14212E+06, 0.14295E+06, 0.14379E+06, 0.14464E+06, 0.14549E+06, 0.14634E+06, 0.14720E+06 &
, 0.14806E+06, 0.14892E+06, 0.14979E+06, 0.15067E+06, 0.15155E+06, 0.15243E+06, 0.15331E+06, 0.15420E+06, 0.15510E+06, 0.15600E+06 &
, 0.15690E+06, 0.15781E+06, 0.15872E+06, 0.15964E+06, 0.16056E+06, 0.16148E+06, 0.16241E+06, 0.16335E+06, 0.16428E+06, 0.16523E+06 &
, 0.16617E+06, 0.16712E+06, 0.16808E+06, 0.16904E+06, 0.17001E+06, 0.17097E+06, 0.17195E+06, 0.17293E+06, 0.17391E+06, 0.17490E+06 &
, 0.17589E+06, 0.17689E+06, 0.17789E+06, 0.17889E+06, 0.17990E+06, 0.18092E+06, 0.18194E+06, 0.18296E+06, 0.18399E+06, 0.18503E+06 &
, 0.18607E+06, 0.18711E+06, 0.18816E+06, 0.18921E+06, 0.19027E+06, 0.19133E+06, 0.19240E+06, 0.19348E+06, 0.19455E+06, 0.19564E+06 &
, 0.19672E+06, 0.19782E+06, 0.19891E+06, 0.20002E+06, 0.20112E+06, 0.20224E+06, 0.20335E+06, 0.20448E+06, 0.20560E+06, 0.20674E+06 &
, 0.20788E+06, 0.20902E+06, 0.21017E+06, 0.21132E+06, 0.21248E+06, 0.21364E+06, 0.21481E+06, 0.21599E+06, 0.21717E+06, 0.21835E+06 &
, 0.21954E+06, 0.22074E+06, 0.22194E+06, 0.22314E+06, 0.22435E+06, 0.22557E+06, 0.22679E+06, 0.22802E+06, 0.22925E+06, 0.23049E+06 &
, 0.23174E+06, 0.23299E+06, 0.23424E+06, 0.23550E+06, 0.23677E+06, 0.23804E+06, 0.23932E+06, 0.24060E+06, 0.24189E+06, 0.24319E+06 /
 !            2          13
 DATA(QofT(          22 ,J),J=1,501)/   17.9780006408691       &
, 0.19759E+02, 0.21540E+02, 0.23321E+02, 0.25102E+02, 0.26883E+02, 0.28664E+02, 0.30445E+02, 0.32226E+02, 0.34007E+02, 0.35788E+02 &
, 0.37570E+02, 0.39351E+02, 0.41132E+02, 0.42913E+02, 0.44694E+02, 0.46475E+02, 0.48257E+02, 0.50038E+02, 0.51819E+02, 0.53600E+02 &
, 0.55382E+02, 0.57163E+02, 0.58944E+02, 0.60725E+02, 0.62507E+02, 0.64288E+02, 0.66070E+02, 0.67851E+02, 0.69633E+02, 0.71415E+02 &
, 0.73197E+02, 0.74979E+02, 0.76761E+02, 0.78543E+02, 0.80326E+02, 0.82109E+02, 0.83893E+02, 0.85677E+02, 0.87461E+02, 0.89247E+02 &
, 0.91032E+02, 0.92819E+02, 0.94606E+02, 0.96395E+02, 0.98184E+02, 0.99974E+02, 0.10177E+03, 0.10356E+03, 0.10535E+03, 0.10715E+03 &
, 0.10895E+03, 0.11075E+03, 0.11255E+03, 0.11435E+03, 0.11616E+03, 0.11797E+03, 0.11978E+03, 0.12160E+03, 0.12341E+03, 0.12523E+03 &
, 0.12706E+03, 0.12888E+03, 0.13071E+03, 0.13255E+03, 0.13439E+03, 0.13623E+03, 0.13808E+03, 0.13993E+03, 0.14178E+03, 0.14365E+03 &
, 0.14551E+03, 0.14738E+03, 0.14926E+03, 0.15114E+03, 0.15303E+03, 0.15492E+03, 0.15682E+03, 0.15873E+03, 0.16064E+03, 0.16256E+03 &
, 0.16449E+03, 0.16642E+03, 0.16836E+03, 0.17030E+03, 0.17226E+03, 0.17422E+03, 0.17619E+03, 0.17817E+03, 0.18015E+03, 0.18215E+03 &
, 0.18415E+03, 0.18616E+03, 0.18818E+03, 0.19021E+03, 0.19225E+03, 0.19430E+03, 0.19635E+03, 0.19842E+03, 0.20050E+03, 0.20258E+03 &
, 0.20468E+03, 0.20678E+03, 0.20890E+03, 0.21102E+03, 0.21316E+03, 0.21531E+03, 0.21747E+03, 0.21964E+03, 0.22182E+03, 0.22401E+03 &
, 0.22621E+03, 0.22843E+03, 0.23065E+03, 0.23289E+03, 0.23514E+03, 0.23740E+03, 0.23968E+03, 0.24197E+03, 0.24426E+03, 0.24658E+03 &
, 0.24890E+03, 0.25124E+03, 0.25359E+03, 0.25595E+03, 0.25833E+03, 0.26072E+03, 0.26312E+03, 0.26554E+03, 0.26797E+03, 0.27041E+03 &
, 0.27287E+03, 0.27534E+03, 0.27783E+03, 0.28033E+03, 0.28285E+03, 0.28538E+03, 0.28792E+03, 0.29048E+03, 0.29305E+03, 0.29564E+03 &
, 0.29825E+03, 0.30087E+03, 0.30350E+03, 0.30615E+03, 0.30882E+03, 0.31150E+03, 0.31420E+03, 0.31691E+03, 0.31964E+03, 0.32238E+03 &
, 0.32514E+03, 0.32792E+03, 0.33072E+03, 0.33353E+03, 0.33635E+03, 0.33920E+03, 0.34206E+03, 0.34494E+03, 0.34783E+03, 0.35075E+03 &
, 0.35368E+03, 0.35662E+03, 0.35959E+03, 0.36257E+03, 0.36557E+03, 0.36859E+03, 0.37163E+03, 0.37468E+03, 0.37776E+03, 0.38085E+03 &
, 0.38396E+03, 0.38709E+03, 0.39023E+03, 0.39340E+03, 0.39658E+03, 0.39979E+03, 0.40301E+03, 0.40625E+03, 0.40951E+03, 0.41279E+03 &
, 0.41609E+03, 0.41941E+03, 0.42275E+03, 0.42611E+03, 0.42949E+03, 0.43289E+03, 0.43631E+03, 0.43975E+03, 0.44322E+03, 0.44670E+03 &
, 0.45020E+03, 0.45372E+03, 0.45727E+03, 0.46083E+03, 0.46442E+03, 0.46803E+03, 0.47165E+03, 0.47531E+03, 0.47898E+03, 0.48267E+03 &
, 0.48639E+03, 0.49012E+03, 0.49388E+03, 0.49767E+03, 0.50147E+03, 0.50530E+03, 0.50915E+03, 0.51302E+03, 0.51691E+03, 0.52083E+03 &
, 0.52477E+03, 0.52873E+03, 0.53272E+03, 0.53673E+03, 0.54076E+03, 0.54482E+03, 0.54890E+03, 0.55301E+03, 0.55713E+03, 0.56129E+03 &
, 0.56546E+03, 0.56966E+03, 0.57389E+03, 0.57814E+03, 0.58242E+03, 0.58671E+03, 0.59104E+03, 0.59539E+03, 0.59976E+03, 0.60416E+03 &
, 0.60859E+03, 0.61304E+03, 0.61752E+03, 0.62202E+03, 0.62655E+03, 0.63110E+03, 0.63568E+03, 0.64029E+03, 0.64492E+03, 0.64958E+03 &
, 0.65427E+03, 0.65898E+03, 0.66372E+03, 0.66849E+03, 0.67328E+03, 0.67811E+03, 0.68295E+03, 0.68783E+03, 0.69274E+03, 0.69767E+03 &
, 0.70263E+03, 0.70762E+03, 0.71263E+03, 0.71768E+03, 0.72275E+03, 0.72785E+03, 0.73298E+03, 0.73814E+03, 0.74333E+03, 0.74855E+03 &
, 0.75380E+03, 0.75907E+03, 0.76438E+03, 0.76972E+03, 0.77508E+03, 0.78048E+03, 0.78590E+03, 0.79136E+03, 0.79685E+03, 0.80236E+03 &
, 0.80791E+03, 0.81349E+03, 0.81910E+03, 0.82474E+03, 0.83041E+03, 0.83611E+03, 0.84185E+03, 0.84761E+03, 0.85341E+03, 0.85924E+03 &
, 0.86510E+03, 0.87100E+03, 0.87692E+03, 0.88288E+03, 0.88887E+03, 0.89490E+03, 0.90095E+03, 0.90705E+03, 0.91317E+03, 0.91933E+03 &
, 0.92552E+03, 0.93174E+03, 0.93800E+03, 0.94429E+03, 0.95062E+03, 0.95698E+03, 0.96337E+03, 0.96980E+03, 0.97627E+03, 0.98276E+03 &
, 0.98930E+03, 0.99587E+03, 0.10025E+04, 0.10091E+04, 0.10158E+04, 0.10225E+04, 0.10292E+04, 0.10360E+04, 0.10428E+04, 0.10497E+04 &
, 0.10566E+04, 0.10635E+04, 0.10705E+04, 0.10775E+04, 0.10845E+04, 0.10916E+04, 0.10987E+04, 0.11059E+04, 0.11131E+04, 0.11203E+04 &
, 0.11276E+04, 0.11349E+04, 0.11422E+04, 0.11496E+04, 0.11570E+04, 0.11645E+04, 0.11720E+04, 0.11795E+04, 0.11871E+04, 0.11947E+04 &
, 0.12024E+04, 0.12101E+04, 0.12178E+04, 0.12256E+04, 0.12334E+04, 0.12413E+04, 0.12492E+04, 0.12571E+04, 0.12651E+04, 0.12731E+04 &
, 0.12812E+04, 0.12893E+04, 0.12974E+04, 0.13056E+04, 0.13138E+04, 0.13221E+04, 0.13304E+04, 0.13388E+04, 0.13472E+04, 0.13556E+04 &
, 0.13641E+04, 0.13727E+04, 0.13812E+04, 0.13899E+04, 0.13985E+04, 0.14072E+04, 0.14160E+04, 0.14248E+04, 0.14336E+04, 0.14425E+04 &
, 0.14514E+04, 0.14604E+04, 0.14694E+04, 0.14785E+04, 0.14876E+04, 0.14968E+04, 0.15060E+04, 0.15152E+04, 0.15245E+04, 0.15338E+04 &
, 0.15432E+04, 0.15527E+04, 0.15622E+04, 0.15717E+04, 0.15813E+04, 0.15909E+04, 0.16006E+04, 0.16103E+04, 0.16201E+04, 0.16299E+04 &
, 0.16397E+04, 0.16497E+04, 0.16596E+04, 0.16696E+04, 0.16797E+04, 0.16898E+04, 0.17000E+04, 0.17102E+04, 0.17205E+04, 0.17308E+04 &
, 0.17411E+04, 0.17516E+04, 0.17620E+04, 0.17725E+04, 0.17831E+04, 0.17937E+04, 0.18044E+04, 0.18151E+04, 0.18259E+04, 0.18367E+04 &
, 0.18476E+04, 0.18585E+04, 0.18695E+04, 0.18806E+04, 0.18917E+04, 0.19028E+04, 0.19140E+04, 0.19253E+04, 0.19366E+04, 0.19479E+04 &
, 0.19594E+04, 0.19708E+04, 0.19824E+04, 0.19939E+04, 0.20056E+04, 0.20173E+04, 0.20290E+04, 0.20408E+04, 0.20527E+04, 0.20646E+04 &
, 0.20766E+04, 0.20886E+04, 0.21007E+04, 0.21129E+04, 0.21251E+04, 0.21373E+04, 0.21497E+04, 0.21620E+04, 0.21745E+04, 0.21870E+04 &
, 0.21995E+04, 0.22121E+04, 0.22248E+04, 0.22375E+04, 0.22503E+04, 0.22632E+04, 0.22761E+04, 0.22891E+04, 0.23021E+04, 0.23152E+04 &
, 0.23283E+04, 0.23416E+04, 0.23548E+04, 0.23682E+04, 0.23816E+04, 0.23950E+04, 0.24086E+04, 0.24222E+04, 0.24358E+04, 0.24495E+04 &
, 0.24633E+04, 0.24771E+04, 0.24910E+04, 0.25050E+04, 0.25190E+04, 0.25331E+04, 0.25473E+04, 0.25615E+04, 0.25758E+04, 0.25902E+04 &
, 0.26046E+04, 0.26191E+04, 0.26336E+04, 0.26483E+04, 0.26630E+04, 0.26777E+04, 0.26925E+04, 0.27074E+04, 0.27224E+04, 0.27374E+04 &
, 0.27525E+04, 0.27677E+04, 0.27829E+04, 0.27982E+04, 0.28135E+04, 0.28290E+04, 0.28445E+04, 0.28601E+04, 0.28757E+04, 0.28914E+04 &
, 0.29072E+04, 0.29231E+04, 0.29390E+04, 0.29550E+04, 0.29710E+04, 0.29872E+04, 0.30034E+04, 0.30197E+04, 0.30360E+04, 0.30525E+04 &
, 0.30690E+04, 0.30855E+04, 0.31022E+04, 0.31189E+04, 0.31357E+04, 0.31526E+04, 0.31695E+04, 0.31865E+04, 0.32036E+04, 0.32208E+04 /
 !            3           1
 DATA(QofT(          23 ,J),J=1,501)/   58.7010002136230       &
, 0.67667E+02, 0.77048E+02, 0.86827E+02, 0.96989E+02, 0.10752E+03, 0.11841E+03, 0.12964E+03, 0.14121E+03, 0.15310E+03, 0.16531E+03 &
, 0.17783E+03, 0.19065E+03, 0.20376E+03, 0.21717E+03, 0.23085E+03, 0.24482E+03, 0.25905E+03, 0.27355E+03, 0.28831E+03, 0.30333E+03 &
, 0.31861E+03, 0.33413E+03, 0.34990E+03, 0.36590E+03, 0.38215E+03, 0.39863E+03, 0.41534E+03, 0.43228E+03, 0.44945E+03, 0.46683E+03 &
, 0.48444E+03, 0.50227E+03, 0.52031E+03, 0.53856E+03, 0.55702E+03, 0.57569E+03, 0.59456E+03, 0.61364E+03, 0.63292E+03, 0.65240E+03 &
, 0.67208E+03, 0.69195E+03, 0.71202E+03, 0.73229E+03, 0.75274E+03, 0.77339E+03, 0.79422E+03, 0.81524E+03, 0.83645E+03, 0.85784E+03 &
, 0.87941E+03, 0.90118E+03, 0.92312E+03, 0.94525E+03, 0.96755E+03, 0.99005E+03, 0.10127E+04, 0.10356E+04, 0.10586E+04, 0.10818E+04 &
, 0.11052E+04, 0.11287E+04, 0.11525E+04, 0.11764E+04, 0.12005E+04, 0.12247E+04, 0.12492E+04, 0.12738E+04, 0.12986E+04, 0.13236E+04 &
, 0.13487E+04, 0.13741E+04, 0.13996E+04, 0.14252E+04, 0.14511E+04, 0.14772E+04, 0.15034E+04, 0.15298E+04, 0.15564E+04, 0.15831E+04 &
, 0.16101E+04, 0.16372E+04, 0.16645E+04, 0.16920E+04, 0.17197E+04, 0.17475E+04, 0.17756E+04, 0.18038E+04, 0.18323E+04, 0.18609E+04 &
, 0.18897E+04, 0.19187E+04, 0.19479E+04, 0.19773E+04, 0.20068E+04, 0.20366E+04, 0.20666E+04, 0.20968E+04, 0.21271E+04, 0.21577E+04 &
, 0.21885E+04, 0.22194E+04, 0.22506E+04, 0.22820E+04, 0.23136E+04, 0.23454E+04, 0.23774E+04, 0.24097E+04, 0.24421E+04, 0.24748E+04 &
, 0.25077E+04, 0.25408E+04, 0.25741E+04, 0.26076E+04, 0.26414E+04, 0.26754E+04, 0.27096E+04, 0.27440E+04, 0.27787E+04, 0.28136E+04 &
, 0.28487E+04, 0.28841E+04, 0.29197E+04, 0.29556E+04, 0.29917E+04, 0.30280E+04, 0.30646E+04, 0.31014E+04, 0.31385E+04, 0.31758E+04 &
, 0.32134E+04, 0.32513E+04, 0.32894E+04, 0.33277E+04, 0.33663E+04, 0.34052E+04, 0.34443E+04, 0.34837E+04, 0.35234E+04, 0.35634E+04 &
, 0.36035E+04, 0.36440E+04, 0.36848E+04, 0.37258E+04, 0.37672E+04, 0.38088E+04, 0.38506E+04, 0.38928E+04, 0.39353E+04, 0.39780E+04 &
, 0.40211E+04, 0.40645E+04, 0.41081E+04, 0.41520E+04, 0.41963E+04, 0.42408E+04, 0.42857E+04, 0.43308E+04, 0.43763E+04, 0.44220E+04 &
, 0.44682E+04, 0.45146E+04, 0.45613E+04, 0.46084E+04, 0.46557E+04, 0.47034E+04, 0.47514E+04, 0.47998E+04, 0.48485E+04, 0.48975E+04 &
, 0.49469E+04, 0.49966E+04, 0.50466E+04, 0.50970E+04, 0.51477E+04, 0.51988E+04, 0.52503E+04, 0.53020E+04, 0.53542E+04, 0.54067E+04 &
, 0.54595E+04, 0.55127E+04, 0.55663E+04, 0.56203E+04, 0.56746E+04, 0.57293E+04, 0.57843E+04, 0.58398E+04, 0.58956E+04, 0.59518E+04 &
, 0.60083E+04, 0.60653E+04, 0.61227E+04, 0.61805E+04, 0.62386E+04, 0.62971E+04, 0.63561E+04, 0.64154E+04, 0.64751E+04, 0.65353E+04 &
, 0.65958E+04, 0.66568E+04, 0.67182E+04, 0.67800E+04, 0.68422E+04, 0.69048E+04, 0.69678E+04, 0.70313E+04, 0.70952E+04, 0.71596E+04 &
, 0.72243E+04, 0.72895E+04, 0.73552E+04, 0.74213E+04, 0.74878E+04, 0.75548E+04, 0.76222E+04, 0.76900E+04, 0.77584E+04, 0.78272E+04 &
, 0.78964E+04, 0.79661E+04, 0.80363E+04, 0.81069E+04, 0.81780E+04, 0.82496E+04, 0.83216E+04, 0.83941E+04, 0.84671E+04, 0.85406E+04 &
, 0.86146E+04, 0.86891E+04, 0.87640E+04, 0.88394E+04, 0.89153E+04, 0.89918E+04, 0.90687E+04, 0.91462E+04, 0.92241E+04, 0.93026E+04 &
, 0.93815E+04, 0.94611E+04, 0.95410E+04, 0.96216E+04, 0.97026E+04, 0.97841E+04, 0.98662E+04, 0.99489E+04, 0.10032E+05, 0.10116E+05 &
, 0.10200E+05, 0.10285E+05, 0.10370E+05, 0.10456E+05, 0.10542E+05, 0.10629E+05, 0.10717E+05, 0.10805E+05, 0.10894E+05, 0.10983E+05 &
, 0.11073E+05, 0.11163E+05, 0.11254E+05, 0.11345E+05, 0.11438E+05, 0.11530E+05, 0.11624E+05, 0.11717E+05, 0.11812E+05, 0.11907E+05 &
, 0.12002E+05, 0.12099E+05, 0.12196E+05, 0.12293E+05, 0.12391E+05, 0.12490E+05, 0.12589E+05, 0.12689E+05, 0.12789E+05, 0.12891E+05 &
, 0.12992E+05, 0.13095E+05, 0.13198E+05, 0.13301E+05, 0.13406E+05, 0.13511E+05, 0.13616E+05, 0.13723E+05, 0.13829E+05, 0.13937E+05 &
, 0.14045E+05, 0.14154E+05, 0.14264E+05, 0.14374E+05, 0.14485E+05, 0.14596E+05, 0.14709E+05, 0.14822E+05, 0.14935E+05, 0.15049E+05 &
, 0.15164E+05, 0.15280E+05, 0.15396E+05, 0.15513E+05, 0.15631E+05, 0.15750E+05, 0.15869E+05, 0.15989E+05, 0.16109E+05, 0.16231E+05 &
, 0.16353E+05, 0.16476E+05, 0.16599E+05, 0.16723E+05, 0.16848E+05, 0.16974E+05, 0.17101E+05, 0.17228E+05, 0.17356E+05, 0.17484E+05 &
, 0.17614E+05, 0.17744E+05, 0.17875E+05, 0.18007E+05, 0.18139E+05, 0.18273E+05, 0.18407E+05, 0.18542E+05, 0.18677E+05, 0.18814E+05 &
, 0.18951E+05, 0.19089E+05, 0.19228E+05, 0.19367E+05, 0.19508E+05, 0.19649E+05, 0.19791E+05, 0.19934E+05, 0.20077E+05, 0.20222E+05 &
, 0.20367E+05, 0.20513E+05, 0.20660E+05, 0.20808E+05, 0.20957E+05, 0.21106E+05, 0.21257E+05, 0.21408E+05, 0.21560E+05, 0.21713E+05 &
, 0.21866E+05, 0.22021E+05, 0.22177E+05, 0.22333E+05, 0.22490E+05, 0.22648E+05, 0.22807E+05, 0.22967E+05, 0.23128E+05, 0.23290E+05 &
, 0.23452E+05, 0.23616E+05, 0.23780E+05, 0.23945E+05, 0.24112E+05, 0.24279E+05, 0.24447E+05, 0.24616E+05, 0.24786E+05, 0.24956E+05 &
, 0.25128E+05, 0.25301E+05, 0.25474E+05, 0.25649E+05, 0.25825E+05, 0.26001E+05, 0.26179E+05, 0.26357E+05, 0.26536E+05, 0.26717E+05 &
, 0.26898E+05, 0.27080E+05, 0.27264E+05, 0.27448E+05, 0.27633E+05, 0.27820E+05, 0.28007E+05, 0.28195E+05, 0.28384E+05, 0.28575E+05 &
, 0.28766E+05, 0.28958E+05, 0.29152E+05, 0.29346E+05, 0.29541E+05, 0.29738E+05, 0.29935E+05, 0.30134E+05, 0.30333E+05, 0.30534E+05 &
, 0.30736E+05, 0.30938E+05, 0.31142E+05, 0.31347E+05, 0.31553E+05, 0.31760E+05, 0.31968E+05, 0.32177E+05, 0.32387E+05, 0.32599E+05 &
, 0.32811E+05, 0.33025E+05, 0.33239E+05, 0.33455E+05, 0.33672E+05, 0.33890E+05, 0.34109E+05, 0.34330E+05, 0.34551E+05, 0.34773E+05 &
, 0.34997E+05, 0.35222E+05, 0.35448E+05, 0.35675E+05, 0.35903E+05, 0.36132E+05, 0.36363E+05, 0.36595E+05, 0.36828E+05, 0.37062E+05 &
, 0.37297E+05, 0.37534E+05, 0.37771E+05, 0.38010E+05, 0.38250E+05, 0.38491E+05, 0.38734E+05, 0.38978E+05, 0.39222E+05, 0.39469E+05 &
, 0.39716E+05, 0.39965E+05, 0.40214E+05, 0.40465E+05, 0.40718E+05, 0.40971E+05, 0.41226E+05, 0.41482E+05, 0.41740E+05, 0.41998E+05 &
, 0.42258E+05, 0.42519E+05, 0.42782E+05, 0.43046E+05, 0.43311E+05, 0.43577E+05, 0.43845E+05, 0.44113E+05, 0.44384E+05, 0.44655E+05 &
, 0.44928E+05, 0.45202E+05, 0.45478E+05, 0.45755E+05, 0.46033E+05, 0.46312E+05, 0.46593E+05, 0.46876E+05, 0.47159E+05, 0.47444E+05 &
, 0.47731E+05, 0.48018E+05, 0.48307E+05, 0.48598E+05, 0.48890E+05, 0.49183E+05, 0.49478E+05, 0.49774E+05, 0.50071E+05, 0.50370E+05 &
, 0.50670E+05, 0.50972E+05, 0.51275E+05, 0.51580E+05, 0.51886E+05, 0.52193E+05, 0.52502E+05, 0.52812E+05, 0.53124E+05, 0.53438E+05 &
, 0.53752E+05, 0.54068E+05, 0.54386E+05, 0.54705E+05, 0.55026E+05, 0.55348E+05, 0.55672E+05, 0.55997E+05, 0.56323E+05, 0.56652E+05 /
 !            3           2
 DATA(QofT(          24 ,J),J=1,501)/   125.290000915527       &
, 0.14443E+03, 0.16446E+03, 0.18534E+03, 0.20704E+03, 0.22952E+03, 0.25277E+03, 0.27675E+03, 0.30145E+03, 0.32684E+03, 0.35291E+03 &
, 0.37964E+03, 0.40702E+03, 0.43502E+03, 0.46364E+03, 0.49287E+03, 0.52268E+03, 0.55308E+03, 0.58404E+03, 0.61556E+03, 0.64763E+03 &
, 0.68024E+03, 0.71339E+03, 0.74705E+03, 0.78124E+03, 0.81593E+03, 0.85112E+03, 0.88680E+03, 0.92297E+03, 0.95962E+03, 0.99675E+03 &
, 0.10344E+04, 0.10724E+04, 0.11109E+04, 0.11499E+04, 0.11893E+04, 0.12292E+04, 0.12695E+04, 0.13102E+04, 0.13514E+04, 0.13930E+04 &
, 0.14350E+04, 0.14775E+04, 0.15203E+04, 0.15636E+04, 0.16073E+04, 0.16514E+04, 0.16959E+04, 0.17408E+04, 0.17861E+04, 0.18318E+04 &
, 0.18779E+04, 0.19244E+04, 0.19712E+04, 0.20185E+04, 0.20661E+04, 0.21142E+04, 0.21626E+04, 0.22114E+04, 0.22606E+04, 0.23102E+04 &
, 0.23602E+04, 0.24105E+04, 0.24613E+04, 0.25124E+04, 0.25639E+04, 0.26157E+04, 0.26680E+04, 0.27206E+04, 0.27736E+04, 0.28271E+04 &
, 0.28808E+04, 0.29350E+04, 0.29896E+04, 0.30445E+04, 0.30998E+04, 0.31556E+04, 0.32116E+04, 0.32681E+04, 0.33250E+04, 0.33823E+04 &
, 0.34400E+04, 0.34981E+04, 0.35565E+04, 0.36154E+04, 0.36747E+04, 0.37344E+04, 0.37944E+04, 0.38549E+04, 0.39158E+04, 0.39771E+04 &
, 0.40389E+04, 0.41010E+04, 0.41636E+04, 0.42266E+04, 0.42900E+04, 0.43538E+04, 0.44181E+04, 0.44828E+04, 0.45479E+04, 0.46135E+04 &
, 0.46795E+04, 0.47460E+04, 0.48129E+04, 0.48803E+04, 0.49481E+04, 0.50164E+04, 0.50852E+04, 0.51544E+04, 0.52240E+04, 0.52942E+04 &
, 0.53648E+04, 0.54359E+04, 0.55075E+04, 0.55796E+04, 0.56522E+04, 0.57252E+04, 0.57988E+04, 0.58728E+04, 0.59474E+04, 0.60225E+04 &
, 0.60981E+04, 0.61742E+04, 0.62509E+04, 0.63280E+04, 0.64057E+04, 0.64839E+04, 0.65627E+04, 0.66419E+04, 0.67218E+04, 0.68022E+04 &
, 0.68831E+04, 0.69647E+04, 0.70467E+04, 0.71294E+04, 0.72126E+04, 0.72963E+04, 0.73807E+04, 0.74657E+04, 0.75512E+04, 0.76374E+04 &
, 0.77241E+04, 0.78114E+04, 0.78994E+04, 0.79879E+04, 0.80771E+04, 0.81668E+04, 0.82572E+04, 0.83483E+04, 0.84400E+04, 0.85323E+04 &
, 0.86252E+04, 0.87189E+04, 0.88132E+04, 0.89081E+04, 0.90036E+04, 0.90998E+04, 0.91968E+04, 0.92944E+04, 0.93926E+04, 0.94916E+04 &
, 0.95913E+04, 0.96916E+04, 0.97927E+04, 0.98945E+04, 0.99969E+04, 0.10100E+05, 0.10204E+05, 0.10309E+05, 0.10414E+05, 0.10520E+05 &
, 0.10627E+05, 0.10735E+05, 0.10843E+05, 0.10952E+05, 0.11062E+05, 0.11172E+05, 0.11284E+05, 0.11396E+05, 0.11509E+05, 0.11623E+05 &
, 0.11737E+05, 0.11853E+05, 0.11969E+05, 0.12086E+05, 0.12203E+05, 0.12322E+05, 0.12441E+05, 0.12562E+05, 0.12682E+05, 0.12804E+05 &
, 0.12927E+05, 0.13051E+05, 0.13175E+05, 0.13300E+05, 0.13427E+05, 0.13554E+05, 0.13681E+05, 0.13810E+05, 0.13940E+05, 0.14071E+05 &
, 0.14202E+05, 0.14334E+05, 0.14468E+05, 0.14602E+05, 0.14737E+05, 0.14873E+05, 0.15010E+05, 0.15148E+05, 0.15287E+05, 0.15426E+05 &
, 0.15567E+05, 0.15709E+05, 0.15851E+05, 0.15995E+05, 0.16140E+05, 0.16285E+05, 0.16432E+05, 0.16579E+05, 0.16728E+05, 0.16878E+05 &
, 0.17028E+05, 0.17180E+05, 0.17332E+05, 0.17486E+05, 0.17641E+05, 0.17796E+05, 0.17953E+05, 0.18111E+05, 0.18270E+05, 0.18430E+05 &
, 0.18591E+05, 0.18753E+05, 0.18916E+05, 0.19080E+05, 0.19245E+05, 0.19412E+05, 0.19579E+05, 0.19748E+05, 0.19918E+05, 0.20089E+05 &
, 0.20260E+05, 0.20434E+05, 0.20608E+05, 0.20783E+05, 0.20960E+05, 0.21138E+05, 0.21316E+05, 0.21497E+05, 0.21678E+05, 0.21860E+05 &
, 0.22044E+05, 0.22229E+05, 0.22415E+05, 0.22602E+05, 0.22790E+05, 0.22980E+05, 0.23171E+05, 0.23363E+05, 0.23556E+05, 0.23751E+05 &
, 0.23947E+05, 0.24144E+05, 0.24342E+05, 0.24542E+05, 0.24743E+05, 0.24945E+05, 0.25148E+05, 0.25353E+05, 0.25559E+05, 0.25767E+05 &
, 0.25975E+05, 0.26185E+05, 0.26397E+05, 0.26609E+05, 0.26823E+05, 0.27039E+05, 0.27256E+05, 0.27474E+05, 0.27693E+05, 0.27914E+05 &
, 0.28136E+05, 0.28360E+05, 0.28585E+05, 0.28811E+05, 0.29039E+05, 0.29268E+05, 0.29499E+05, 0.29731E+05, 0.29965E+05, 0.30200E+05 &
, 0.30436E+05, 0.30674E+05, 0.30914E+05, 0.31154E+05, 0.31397E+05, 0.31641E+05, 0.31886E+05, 0.32133E+05, 0.32381E+05, 0.32631E+05 &
, 0.32882E+05, 0.33135E+05, 0.33389E+05, 0.33645E+05, 0.33903E+05, 0.34162E+05, 0.34423E+05, 0.34685E+05, 0.34949E+05, 0.35214E+05 &
, 0.35481E+05, 0.35750E+05, 0.36020E+05, 0.36291E+05, 0.36565E+05, 0.36840E+05, 0.37117E+05, 0.37395E+05, 0.37675E+05, 0.37957E+05 &
, 0.38240E+05, 0.38525E+05, 0.38812E+05, 0.39100E+05, 0.39390E+05, 0.39682E+05, 0.39975E+05, 0.40270E+05, 0.40567E+05, 0.40866E+05 &
, 0.41166E+05, 0.41469E+05, 0.41772E+05, 0.42078E+05, 0.42386E+05, 0.42695E+05, 0.43006E+05, 0.43319E+05, 0.43633E+05, 0.43950E+05 &
, 0.44268E+05, 0.44588E+05, 0.44910E+05, 0.45234E+05, 0.45559E+05, 0.45887E+05, 0.46216E+05, 0.46548E+05, 0.46881E+05, 0.47216E+05 &
, 0.47553E+05, 0.47892E+05, 0.48232E+05, 0.48575E+05, 0.48920E+05, 0.49266E+05, 0.49615E+05, 0.49965E+05, 0.50318E+05, 0.50672E+05 &
, 0.51028E+05, 0.51387E+05, 0.51747E+05, 0.52110E+05, 0.52474E+05, 0.52841E+05, 0.53209E+05, 0.53579E+05, 0.53952E+05, 0.54327E+05 &
, 0.54703E+05, 0.55082E+05, 0.55463E+05, 0.55846E+05, 0.56231E+05, 0.56618E+05, 0.57007E+05, 0.57399E+05, 0.57792E+05, 0.58188E+05 &
, 0.58586E+05, 0.58986E+05, 0.59388E+05, 0.59792E+05, 0.60199E+05, 0.60608E+05, 0.61019E+05, 0.61432E+05, 0.61847E+05, 0.62265E+05 &
, 0.62685E+05, 0.63107E+05, 0.63531E+05, 0.63958E+05, 0.64387E+05, 0.64818E+05, 0.65251E+05, 0.65687E+05, 0.66126E+05, 0.66566E+05 &
, 0.67009E+05, 0.67454E+05, 0.67901E+05, 0.68351E+05, 0.68803E+05, 0.69258E+05, 0.69715E+05, 0.70174E+05, 0.70636E+05, 0.71100E+05 &
, 0.71567E+05, 0.72036E+05, 0.72507E+05, 0.72981E+05, 0.73458E+05, 0.73937E+05, 0.74418E+05, 0.74902E+05, 0.75388E+05, 0.75877E+05 &
, 0.76368E+05, 0.76862E+05, 0.77359E+05, 0.77858E+05, 0.78359E+05, 0.78863E+05, 0.79370E+05, 0.79879E+05, 0.80391E+05, 0.80906E+05 &
, 0.81423E+05, 0.81943E+05, 0.82465E+05, 0.82990E+05, 0.83518E+05, 0.84048E+05, 0.84581E+05, 0.85117E+05, 0.85655E+05, 0.86196E+05 &
, 0.86740E+05, 0.87287E+05, 0.87836E+05, 0.88388E+05, 0.88943E+05, 0.89500E+05, 0.90061E+05, 0.90624E+05, 0.91190E+05, 0.91759E+05 &
, 0.92330E+05, 0.92905E+05, 0.93482E+05, 0.94062E+05, 0.94645E+05, 0.95231E+05, 0.95819E+05, 0.96411E+05, 0.97005E+05, 0.97603E+05 &
, 0.98203E+05, 0.98806E+05, 0.99413E+05, 0.10002E+06, 0.10063E+06, 0.10125E+06, 0.10187E+06, 0.10249E+06, 0.10311E+06, 0.10374E+06 &
, 0.10437E+06, 0.10500E+06, 0.10564E+06, 0.10628E+06, 0.10692E+06, 0.10757E+06, 0.10821E+06, 0.10887E+06, 0.10952E+06, 0.11018E+06 &
, 0.11084E+06, 0.11150E+06, 0.11217E+06, 0.11284E+06, 0.11351E+06, 0.11419E+06, 0.11487E+06, 0.11555E+06, 0.11624E+06, 0.11693E+06 &
, 0.11762E+06, 0.11832E+06, 0.11902E+06, 0.11972E+06, 0.12043E+06, 0.12114E+06, 0.12185E+06, 0.12257E+06, 0.12329E+06, 0.12401E+06 /
 !            3           3
 DATA(QofT(          25 ,J),J=1,501)/   61.2560005187988       &
, 0.70613E+02, 0.80403E+02, 0.90609E+02, 0.10121E+03, 0.11220E+03, 0.12357E+03, 0.13529E+03, 0.14736E+03, 0.15977E+03, 0.17251E+03 &
, 0.18558E+03, 0.19896E+03, 0.21264E+03, 0.22663E+03, 0.24092E+03, 0.25549E+03, 0.27034E+03, 0.28548E+03, 0.30088E+03, 0.31656E+03 &
, 0.33250E+03, 0.34869E+03, 0.36515E+03, 0.38185E+03, 0.39881E+03, 0.41601E+03, 0.43345E+03, 0.45113E+03, 0.46904E+03, 0.48719E+03 &
, 0.50557E+03, 0.52417E+03, 0.54299E+03, 0.56204E+03, 0.58131E+03, 0.60079E+03, 0.62048E+03, 0.64040E+03, 0.66052E+03, 0.68085E+03 &
, 0.70139E+03, 0.72213E+03, 0.74307E+03, 0.76422E+03, 0.78557E+03, 0.80712E+03, 0.82886E+03, 0.85080E+03, 0.87294E+03, 0.89526E+03 &
, 0.91779E+03, 0.94050E+03, 0.96341E+03, 0.98650E+03, 0.10098E+04, 0.10333E+04, 0.10569E+04, 0.10808E+04, 0.11048E+04, 0.11290E+04 &
, 0.11535E+04, 0.11781E+04, 0.12028E+04, 0.12278E+04, 0.12530E+04, 0.12783E+04, 0.13038E+04, 0.13295E+04, 0.13554E+04, 0.13815E+04 &
, 0.14078E+04, 0.14343E+04, 0.14609E+04, 0.14878E+04, 0.15148E+04, 0.15420E+04, 0.15694E+04, 0.15970E+04, 0.16248E+04, 0.16527E+04 &
, 0.16809E+04, 0.17093E+04, 0.17378E+04, 0.17666E+04, 0.17955E+04, 0.18247E+04, 0.18540E+04, 0.18835E+04, 0.19133E+04, 0.19432E+04 &
, 0.19733E+04, 0.20037E+04, 0.20342E+04, 0.20650E+04, 0.20960E+04, 0.21271E+04, 0.21585E+04, 0.21901E+04, 0.22219E+04, 0.22539E+04 &
, 0.22862E+04, 0.23186E+04, 0.23513E+04, 0.23842E+04, 0.24173E+04, 0.24507E+04, 0.24843E+04, 0.25180E+04, 0.25521E+04, 0.25863E+04 &
, 0.26208E+04, 0.26555E+04, 0.26905E+04, 0.27257E+04, 0.27611E+04, 0.27968E+04, 0.28328E+04, 0.28689E+04, 0.29054E+04, 0.29420E+04 &
, 0.29789E+04, 0.30161E+04, 0.30536E+04, 0.30912E+04, 0.31292E+04, 0.31674E+04, 0.32059E+04, 0.32446E+04, 0.32836E+04, 0.33229E+04 &
, 0.33624E+04, 0.34023E+04, 0.34423E+04, 0.34827E+04, 0.35234E+04, 0.35643E+04, 0.36056E+04, 0.36471E+04, 0.36889E+04, 0.37310E+04 &
, 0.37733E+04, 0.38161E+04, 0.38590E+04, 0.39023E+04, 0.39459E+04, 0.39898E+04, 0.40340E+04, 0.40785E+04, 0.41233E+04, 0.41684E+04 &
, 0.42139E+04, 0.42596E+04, 0.43057E+04, 0.43522E+04, 0.43989E+04, 0.44460E+04, 0.44934E+04, 0.45411E+04, 0.45892E+04, 0.46376E+04 &
, 0.46863E+04, 0.47354E+04, 0.47849E+04, 0.48346E+04, 0.48848E+04, 0.49352E+04, 0.49861E+04, 0.50373E+04, 0.50888E+04, 0.51408E+04 &
, 0.51930E+04, 0.52457E+04, 0.52987E+04, 0.53521E+04, 0.54059E+04, 0.54600E+04, 0.55145E+04, 0.55695E+04, 0.56247E+04, 0.56804E+04 &
, 0.57365E+04, 0.57929E+04, 0.58498E+04, 0.59071E+04, 0.59647E+04, 0.60228E+04, 0.60813E+04, 0.61401E+04, 0.61994E+04, 0.62591E+04 &
, 0.63192E+04, 0.63798E+04, 0.64407E+04, 0.65021E+04, 0.65639E+04, 0.66261E+04, 0.66887E+04, 0.67518E+04, 0.68154E+04, 0.68794E+04 &
, 0.69438E+04, 0.70086E+04, 0.70740E+04, 0.71397E+04, 0.72059E+04, 0.72726E+04, 0.73397E+04, 0.74073E+04, 0.74754E+04, 0.75439E+04 &
, 0.76129E+04, 0.76823E+04, 0.77523E+04, 0.78227E+04, 0.78936E+04, 0.79650E+04, 0.80368E+04, 0.81092E+04, 0.81821E+04, 0.82555E+04 &
, 0.83293E+04, 0.84036E+04, 0.84785E+04, 0.85538E+04, 0.86297E+04, 0.87061E+04, 0.87830E+04, 0.88604E+04, 0.89383E+04, 0.90168E+04 &
, 0.90958E+04, 0.91753E+04, 0.92553E+04, 0.93359E+04, 0.94171E+04, 0.94987E+04, 0.95809E+04, 0.96636E+04, 0.97470E+04, 0.98309E+04 &
, 0.99153E+04, 0.10000E+05, 0.10086E+05, 0.10172E+05, 0.10259E+05, 0.10346E+05, 0.10434E+05, 0.10522E+05, 0.10611E+05, 0.10701E+05 &
, 0.10791E+05, 0.10881E+05, 0.10973E+05, 0.11065E+05, 0.11157E+05, 0.11250E+05, 0.11344E+05, 0.11439E+05, 0.11533E+05, 0.11629E+05 &
, 0.11725E+05, 0.11822E+05, 0.11919E+05, 0.12017E+05, 0.12116E+05, 0.12216E+05, 0.12315E+05, 0.12416E+05, 0.12517E+05, 0.12619E+05 &
, 0.12722E+05, 0.12825E+05, 0.12929E+05, 0.13033E+05, 0.13139E+05, 0.13244E+05, 0.13351E+05, 0.13458E+05, 0.13566E+05, 0.13675E+05 &
, 0.13784E+05, 0.13894E+05, 0.14004E+05, 0.14116E+05, 0.14228E+05, 0.14340E+05, 0.14454E+05, 0.14568E+05, 0.14683E+05, 0.14798E+05 &
, 0.14915E+05, 0.15032E+05, 0.15149E+05, 0.15268E+05, 0.15387E+05, 0.15507E+05, 0.15627E+05, 0.15749E+05, 0.15871E+05, 0.15994E+05 &
, 0.16117E+05, 0.16242E+05, 0.16367E+05, 0.16493E+05, 0.16619E+05, 0.16747E+05, 0.16875E+05, 0.17004E+05, 0.17134E+05, 0.17264E+05 &
, 0.17396E+05, 0.17528E+05, 0.17661E+05, 0.17795E+05, 0.17929E+05, 0.18064E+05, 0.18201E+05, 0.18337E+05, 0.18475E+05, 0.18614E+05 &
, 0.18753E+05, 0.18894E+05, 0.19035E+05, 0.19177E+05, 0.19319E+05, 0.19463E+05, 0.19607E+05, 0.19753E+05, 0.19899E+05, 0.20046E+05 &
, 0.20194E+05, 0.20342E+05, 0.20492E+05, 0.20643E+05, 0.20794E+05, 0.20946E+05, 0.21099E+05, 0.21253E+05, 0.21408E+05, 0.21564E+05 &
, 0.21721E+05, 0.21878E+05, 0.22037E+05, 0.22196E+05, 0.22357E+05, 0.22518E+05, 0.22680E+05, 0.22843E+05, 0.23007E+05, 0.23172E+05 &
, 0.23338E+05, 0.23505E+05, 0.23673E+05, 0.23842E+05, 0.24011E+05, 0.24182E+05, 0.24354E+05, 0.24526E+05, 0.24700E+05, 0.24875E+05 &
, 0.25050E+05, 0.25227E+05, 0.25404E+05, 0.25583E+05, 0.25762E+05, 0.25943E+05, 0.26124E+05, 0.26307E+05, 0.26490E+05, 0.26675E+05 &
, 0.26861E+05, 0.27047E+05, 0.27235E+05, 0.27424E+05, 0.27613E+05, 0.27804E+05, 0.27996E+05, 0.28189E+05, 0.28383E+05, 0.28578E+05 &
, 0.28774E+05, 0.28971E+05, 0.29169E+05, 0.29369E+05, 0.29569E+05, 0.29770E+05, 0.29973E+05, 0.30177E+05, 0.30381E+05, 0.30587E+05 &
, 0.30794E+05, 0.31002E+05, 0.31211E+05, 0.31422E+05, 0.31633E+05, 0.31846E+05, 0.32060E+05, 0.32274E+05, 0.32490E+05, 0.32708E+05 &
, 0.32926E+05, 0.33145E+05, 0.33366E+05, 0.33588E+05, 0.33811E+05, 0.34035E+05, 0.34260E+05, 0.34487E+05, 0.34715E+05, 0.34944E+05 &
, 0.35174E+05, 0.35405E+05, 0.35637E+05, 0.35871E+05, 0.36106E+05, 0.36342E+05, 0.36580E+05, 0.36818E+05, 0.37058E+05, 0.37299E+05 &
, 0.37542E+05, 0.37785E+05, 0.38030E+05, 0.38276E+05, 0.38524E+05, 0.38772E+05, 0.39023E+05, 0.39274E+05, 0.39526E+05, 0.39780E+05 &
, 0.40035E+05, 0.40292E+05, 0.40549E+05, 0.40808E+05, 0.41069E+05, 0.41330E+05, 0.41593E+05, 0.41858E+05, 0.42123E+05, 0.42390E+05 &
, 0.42659E+05, 0.42928E+05, 0.43200E+05, 0.43472E+05, 0.43746E+05, 0.44021E+05, 0.44297E+05, 0.44575E+05, 0.44854E+05, 0.45135E+05 &
, 0.45417E+05, 0.45701E+05, 0.45986E+05, 0.46272E+05, 0.46560E+05, 0.46849E+05, 0.47139E+05, 0.47431E+05, 0.47725E+05, 0.48019E+05 &
, 0.48316E+05, 0.48614E+05, 0.48913E+05, 0.49213E+05, 0.49516E+05, 0.49819E+05, 0.50124E+05, 0.50431E+05, 0.50739E+05, 0.51049E+05 &
, 0.51360E+05, 0.51672E+05, 0.51986E+05, 0.52302E+05, 0.52619E+05, 0.52938E+05, 0.53258E+05, 0.53580E+05, 0.53903E+05, 0.54228E+05 &
, 0.54554E+05, 0.54882E+05, 0.55211E+05, 0.55543E+05, 0.55875E+05, 0.56209E+05, 0.56545E+05, 0.56883E+05, 0.57222E+05, 0.57562E+05 &
, 0.57904E+05, 0.58248E+05, 0.58594E+05, 0.58941E+05, 0.59289E+05, 0.59640E+05, 0.59992E+05, 0.60345E+05, 0.60701E+05, 0.61058E+05 /
 !            3           4
 DATA(QofT(          26 ,J),J=1,501)/   728.609985351562       &
, 0.83992E+03, 0.95638E+03, 0.10778E+04, 0.12040E+04, 0.13347E+04, 0.14699E+04, 0.16093E+04, 0.17529E+04, 0.19006E+04, 0.20521E+04 &
, 0.22076E+04, 0.23667E+04, 0.25296E+04, 0.26960E+04, 0.28659E+04, 0.30392E+04, 0.32160E+04, 0.33960E+04, 0.35793E+04, 0.37657E+04 &
, 0.39554E+04, 0.41481E+04, 0.43438E+04, 0.45426E+04, 0.47443E+04, 0.49489E+04, 0.51563E+04, 0.53666E+04, 0.55798E+04, 0.57956E+04 &
, 0.60143E+04, 0.62356E+04, 0.64595E+04, 0.66861E+04, 0.69153E+04, 0.71471E+04, 0.73815E+04, 0.76183E+04, 0.78577E+04, 0.80996E+04 &
, 0.83439E+04, 0.85906E+04, 0.88398E+04, 0.90914E+04, 0.93454E+04, 0.96017E+04, 0.98603E+04, 0.10121E+05, 0.10385E+05, 0.10650E+05 &
, 0.10918E+05, 0.11188E+05, 0.11461E+05, 0.11736E+05, 0.12013E+05, 0.12292E+05, 0.12574E+05, 0.12857E+05, 0.13143E+05, 0.13431E+05 &
, 0.13722E+05, 0.14014E+05, 0.14309E+05, 0.14606E+05, 0.14905E+05, 0.15207E+05, 0.15511E+05, 0.15817E+05, 0.16125E+05, 0.16435E+05 &
, 0.16747E+05, 0.17062E+05, 0.17379E+05, 0.17698E+05, 0.18020E+05, 0.18343E+05, 0.18669E+05, 0.18997E+05, 0.19328E+05, 0.19660E+05 &
, 0.19995E+05, 0.20332E+05, 0.20672E+05, 0.21014E+05, 0.21358E+05, 0.21704E+05, 0.22053E+05, 0.22404E+05, 0.22758E+05, 0.23114E+05 &
, 0.23472E+05, 0.23833E+05, 0.24196E+05, 0.24561E+05, 0.24929E+05, 0.25300E+05, 0.25673E+05, 0.26048E+05, 0.26426E+05, 0.26807E+05 &
, 0.27190E+05, 0.27575E+05, 0.27963E+05, 0.28354E+05, 0.28747E+05, 0.29143E+05, 0.29542E+05, 0.29943E+05, 0.30347E+05, 0.30754E+05 &
, 0.31163E+05, 0.31575E+05, 0.31991E+05, 0.32408E+05, 0.32829E+05, 0.33252E+05, 0.33679E+05, 0.34108E+05, 0.34540E+05, 0.34975E+05 &
, 0.35413E+05, 0.35854E+05, 0.36297E+05, 0.36744E+05, 0.37195E+05, 0.37647E+05, 0.38103E+05, 0.38563E+05, 0.39025E+05, 0.39490E+05 &
, 0.39959E+05, 0.40431E+05, 0.40906E+05, 0.41384E+05, 0.41866E+05, 0.42351E+05, 0.42839E+05, 0.43331E+05, 0.43826E+05, 0.44324E+05 &
, 0.44826E+05, 0.45331E+05, 0.45840E+05, 0.46352E+05, 0.46868E+05, 0.47387E+05, 0.47911E+05, 0.48437E+05, 0.48967E+05, 0.49501E+05 &
, 0.50039E+05, 0.50580E+05, 0.51125E+05, 0.51674E+05, 0.52226E+05, 0.52783E+05, 0.53343E+05, 0.53907E+05, 0.54475E+05, 0.55047E+05 &
, 0.55623E+05, 0.56203E+05, 0.56787E+05, 0.57375E+05, 0.57967E+05, 0.58563E+05, 0.59164E+05, 0.59768E+05, 0.60377E+05, 0.60990E+05 &
, 0.61607E+05, 0.62228E+05, 0.62854E+05, 0.63484E+05, 0.64118E+05, 0.64757E+05, 0.65400E+05, 0.66048E+05, 0.66700E+05, 0.67357E+05 &
, 0.68018E+05, 0.68683E+05, 0.69354E+05, 0.70029E+05, 0.70708E+05, 0.71393E+05, 0.72082E+05, 0.72776E+05, 0.73474E+05, 0.74177E+05 &
, 0.74886E+05, 0.75599E+05, 0.76316E+05, 0.77040E+05, 0.77768E+05, 0.78500E+05, 0.79238E+05, 0.79981E+05, 0.80729E+05, 0.81482E+05 &
, 0.82240E+05, 0.83003E+05, 0.83772E+05, 0.84546E+05, 0.85326E+05, 0.86110E+05, 0.86899E+05, 0.87695E+05, 0.88495E+05, 0.89301E+05 &
, 0.90112E+05, 0.90929E+05, 0.91752E+05, 0.92579E+05, 0.93414E+05, 0.94253E+05, 0.95097E+05, 0.95948E+05, 0.96804E+05, 0.97666E+05 &
, 0.98534E+05, 0.99408E+05, 0.10029E+06, 0.10117E+06, 0.10206E+06, 0.10296E+06, 0.10386E+06, 0.10477E+06, 0.10569E+06, 0.10661E+06 &
, 0.10754E+06, 0.10847E+06, 0.10941E+06, 0.11036E+06, 0.11131E+06, 0.11227E+06, 0.11323E+06, 0.11420E+06, 0.11518E+06, 0.11617E+06 &
, 0.11716E+06, 0.11815E+06, 0.11916E+06, 0.12017E+06, 0.12118E+06, 0.12221E+06, 0.12324E+06, 0.12427E+06, 0.12532E+06, 0.12637E+06 &
, 0.12742E+06, 0.12849E+06, 0.12956E+06, 0.13064E+06, 0.13172E+06, 0.13281E+06, 0.13391E+06, 0.13502E+06, 0.13613E+06, 0.13725E+06 &
, 0.13838E+06, 0.13951E+06, 0.14065E+06, 0.14180E+06, 0.14296E+06, 0.14412E+06, 0.14529E+06, 0.14647E+06, 0.14765E+06, 0.14885E+06 &
, 0.15005E+06, 0.15126E+06, 0.15247E+06, 0.15370E+06, 0.15493E+06, 0.15617E+06, 0.15741E+06, 0.15867E+06, 0.15993E+06, 0.16120E+06 &
, 0.16248E+06, 0.16377E+06, 0.16506E+06, 0.16636E+06, 0.16767E+06, 0.16899E+06, 0.17032E+06, 0.17165E+06, 0.17300E+06, 0.17435E+06 &
, 0.17571E+06, 0.17707E+06, 0.17845E+06, 0.17984E+06, 0.18123E+06, 0.18263E+06, 0.18404E+06, 0.18546E+06, 0.18689E+06, 0.18832E+06 &
, 0.18977E+06, 0.19122E+06, 0.19269E+06, 0.19416E+06, 0.19564E+06, 0.19713E+06, 0.19862E+06, 0.20013E+06, 0.20165E+06, 0.20317E+06 &
, 0.20471E+06, 0.20625E+06, 0.20780E+06, 0.20937E+06, 0.21094E+06, 0.21252E+06, 0.21411E+06, 0.21571E+06, 0.21732E+06, 0.21894E+06 &
, 0.22056E+06, 0.22220E+06, 0.22385E+06, 0.22550E+06, 0.22717E+06, 0.22885E+06, 0.23053E+06, 0.23223E+06, 0.23393E+06, 0.23565E+06 &
, 0.23738E+06, 0.23911E+06, 0.24086E+06, 0.24261E+06, 0.24438E+06, 0.24616E+06, 0.24794E+06, 0.24974E+06, 0.25155E+06, 0.25336E+06 &
, 0.25519E+06, 0.25703E+06, 0.25888E+06, 0.26074E+06, 0.26261E+06, 0.26449E+06, 0.26638E+06, 0.26828E+06, 0.27020E+06, 0.27212E+06 &
, 0.27406E+06, 0.27600E+06, 0.27796E+06, 0.27993E+06, 0.28190E+06, 0.28389E+06, 0.28590E+06, 0.28791E+06, 0.28993E+06, 0.29197E+06 &
, 0.29401E+06, 0.29607E+06, 0.29814E+06, 0.30022E+06, 0.30231E+06, 0.30441E+06, 0.30653E+06, 0.30866E+06, 0.31079E+06, 0.31295E+06 &
, 0.31511E+06, 0.31728E+06, 0.31947E+06, 0.32166E+06, 0.32387E+06, 0.32610E+06, 0.32833E+06, 0.33058E+06, 0.33284E+06, 0.33511E+06 &
, 0.33739E+06, 0.33969E+06, 0.34199E+06, 0.34431E+06, 0.34665E+06, 0.34899E+06, 0.35135E+06, 0.35372E+06, 0.35611E+06, 0.35850E+06 &
, 0.36091E+06, 0.36333E+06, 0.36577E+06, 0.36822E+06, 0.37068E+06, 0.37315E+06, 0.37564E+06, 0.37814E+06, 0.38065E+06, 0.38318E+06 &
, 0.38572E+06, 0.38827E+06, 0.39084E+06, 0.39342E+06, 0.39602E+06, 0.39862E+06, 0.40124E+06, 0.40388E+06, 0.40653E+06, 0.40919E+06 &
, 0.41187E+06, 0.41456E+06, 0.41726E+06, 0.41998E+06, 0.42271E+06, 0.42546E+06, 0.42822E+06, 0.43100E+06, 0.43378E+06, 0.43659E+06 &
, 0.43941E+06, 0.44224E+06, 0.44509E+06, 0.44795E+06, 0.45082E+06, 0.45371E+06, 0.45662E+06, 0.45954E+06, 0.46248E+06, 0.46543E+06 &
, 0.46839E+06, 0.47137E+06, 0.47437E+06, 0.47738E+06, 0.48040E+06, 0.48344E+06, 0.48650E+06, 0.48957E+06, 0.49266E+06, 0.49576E+06 &
, 0.49888E+06, 0.50201E+06, 0.50516E+06, 0.50833E+06, 0.51151E+06, 0.51470E+06, 0.51792E+06, 0.52115E+06, 0.52439E+06, 0.52765E+06 &
, 0.53093E+06, 0.53422E+06, 0.53753E+06, 0.54085E+06, 0.54419E+06, 0.54755E+06, 0.55092E+06, 0.55432E+06, 0.55772E+06, 0.56115E+06 &
, 0.56459E+06, 0.56805E+06, 0.57152E+06, 0.57501E+06, 0.57852E+06, 0.58204E+06, 0.58559E+06, 0.58915E+06, 0.59272E+06, 0.59632E+06 &
, 0.59993E+06, 0.60355E+06, 0.60720E+06, 0.61086E+06, 0.61454E+06, 0.61824E+06, 0.62196E+06, 0.62569E+06, 0.62945E+06, 0.63321E+06 &
, 0.63700E+06, 0.64081E+06, 0.64463E+06, 0.64847E+06, 0.65233E+06, 0.65621E+06, 0.66011E+06, 0.66402E+06, 0.66795E+06, 0.67191E+06 &
, 0.67588E+06, 0.67987E+06, 0.68387E+06, 0.68790E+06, 0.69194E+06, 0.69601E+06, 0.70009E+06, 0.70419E+06, 0.70831E+06, 0.71245E+06 /
 !            3           5
 DATA(QofT(          27 ,J),J=1,501)/   360.089996337891       &
, 0.41509E+03, 0.47264E+03, 0.53263E+03, 0.59497E+03, 0.65957E+03, 0.72636E+03, 0.79526E+03, 0.86622E+03, 0.93918E+03, 0.10141E+04 &
, 0.10909E+04, 0.11695E+04, 0.12500E+04, 0.13322E+04, 0.14162E+04, 0.15018E+04, 0.15891E+04, 0.16781E+04, 0.17687E+04, 0.18608E+04 &
, 0.19545E+04, 0.20497E+04, 0.21464E+04, 0.22446E+04, 0.23443E+04, 0.24454E+04, 0.25479E+04, 0.26518E+04, 0.27571E+04, 0.28638E+04 &
, 0.29718E+04, 0.30812E+04, 0.31918E+04, 0.33038E+04, 0.34170E+04, 0.35316E+04, 0.36474E+04, 0.37644E+04, 0.38827E+04, 0.40022E+04 &
, 0.41229E+04, 0.42448E+04, 0.43679E+04, 0.44922E+04, 0.46177E+04, 0.47444E+04, 0.48722E+04, 0.50011E+04, 0.51313E+04, 0.52625E+04 &
, 0.53949E+04, 0.55284E+04, 0.56630E+04, 0.57988E+04, 0.59356E+04, 0.60736E+04, 0.62127E+04, 0.63529E+04, 0.64941E+04, 0.66365E+04 &
, 0.67800E+04, 0.69245E+04, 0.70702E+04, 0.72170E+04, 0.73647E+04, 0.75137E+04, 0.76637E+04, 0.78148E+04, 0.79670E+04, 0.81202E+04 &
, 0.82746E+04, 0.84300E+04, 0.85866E+04, 0.87443E+04, 0.89030E+04, 0.90629E+04, 0.92239E+04, 0.93860E+04, 0.95491E+04, 0.97134E+04 &
, 0.98789E+04, 0.10045E+05, 0.10213E+05, 0.10382E+05, 0.10552E+05, 0.10723E+05, 0.10895E+05, 0.11069E+05, 0.11243E+05, 0.11419E+05 &
, 0.11596E+05, 0.11774E+05, 0.11954E+05, 0.12134E+05, 0.12316E+05, 0.12499E+05, 0.12683E+05, 0.12868E+05, 0.13055E+05, 0.13243E+05 &
, 0.13432E+05, 0.13622E+05, 0.13814E+05, 0.14007E+05, 0.14201E+05, 0.14397E+05, 0.14594E+05, 0.14792E+05, 0.14991E+05, 0.15192E+05 &
, 0.15395E+05, 0.15598E+05, 0.15803E+05, 0.16009E+05, 0.16217E+05, 0.16426E+05, 0.16637E+05, 0.16849E+05, 0.17062E+05, 0.17277E+05 &
, 0.17493E+05, 0.17711E+05, 0.17930E+05, 0.18151E+05, 0.18373E+05, 0.18597E+05, 0.18822E+05, 0.19049E+05, 0.19278E+05, 0.19508E+05 &
, 0.19739E+05, 0.19972E+05, 0.20207E+05, 0.20443E+05, 0.20681E+05, 0.20921E+05, 0.21162E+05, 0.21405E+05, 0.21649E+05, 0.21896E+05 &
, 0.22144E+05, 0.22393E+05, 0.22645E+05, 0.22898E+05, 0.23153E+05, 0.23409E+05, 0.23668E+05, 0.23928E+05, 0.24190E+05, 0.24454E+05 &
, 0.24719E+05, 0.24987E+05, 0.25256E+05, 0.25527E+05, 0.25801E+05, 0.26076E+05, 0.26353E+05, 0.26631E+05, 0.26912E+05, 0.27195E+05 &
, 0.27480E+05, 0.27766E+05, 0.28055E+05, 0.28346E+05, 0.28638E+05, 0.28933E+05, 0.29230E+05, 0.29529E+05, 0.29830E+05, 0.30133E+05 &
, 0.30438E+05, 0.30745E+05, 0.31054E+05, 0.31366E+05, 0.31679E+05, 0.31995E+05, 0.32314E+05, 0.32634E+05, 0.32956E+05, 0.33281E+05 &
, 0.33608E+05, 0.33937E+05, 0.34269E+05, 0.34602E+05, 0.34939E+05, 0.35277E+05, 0.35618E+05, 0.35961E+05, 0.36306E+05, 0.36654E+05 &
, 0.37005E+05, 0.37357E+05, 0.37713E+05, 0.38070E+05, 0.38430E+05, 0.38793E+05, 0.39158E+05, 0.39525E+05, 0.39895E+05, 0.40268E+05 &
, 0.40643E+05, 0.41021E+05, 0.41401E+05, 0.41784E+05, 0.42169E+05, 0.42557E+05, 0.42948E+05, 0.43342E+05, 0.43738E+05, 0.44137E+05 &
, 0.44538E+05, 0.44943E+05, 0.45349E+05, 0.45759E+05, 0.46172E+05, 0.46587E+05, 0.47005E+05, 0.47426E+05, 0.47850E+05, 0.48277E+05 &
, 0.48706E+05, 0.49139E+05, 0.49574E+05, 0.50012E+05, 0.50453E+05, 0.50898E+05, 0.51345E+05, 0.51795E+05, 0.52248E+05, 0.52704E+05 &
, 0.53163E+05, 0.53625E+05, 0.54091E+05, 0.54559E+05, 0.55031E+05, 0.55505E+05, 0.55983E+05, 0.56464E+05, 0.56948E+05, 0.57435E+05 &
, 0.57925E+05, 0.58419E+05, 0.58916E+05, 0.59416E+05, 0.59920E+05, 0.60427E+05, 0.60937E+05, 0.61450E+05, 0.61967E+05, 0.62487E+05 &
, 0.63010E+05, 0.63537E+05, 0.64068E+05, 0.64602E+05, 0.65139E+05, 0.65680E+05, 0.66224E+05, 0.66772E+05, 0.67323E+05, 0.67877E+05 &
, 0.68436E+05, 0.68998E+05, 0.69563E+05, 0.70132E+05, 0.70705E+05, 0.71281E+05, 0.71862E+05, 0.72445E+05, 0.73033E+05, 0.73624E+05 &
, 0.74219E+05, 0.74817E+05, 0.75420E+05, 0.76026E+05, 0.76637E+05, 0.77251E+05, 0.77868E+05, 0.78490E+05, 0.79116E+05, 0.79745E+05 &
, 0.80379E+05, 0.81016E+05, 0.81658E+05, 0.82303E+05, 0.82953E+05, 0.83606E+05, 0.84263E+05, 0.84925E+05, 0.85590E+05, 0.86260E+05 &
, 0.86934E+05, 0.87612E+05, 0.88294E+05, 0.88981E+05, 0.89671E+05, 0.90366E+05, 0.91066E+05, 0.91769E+05, 0.92476E+05, 0.93189E+05 &
, 0.93905E+05, 0.94625E+05, 0.95350E+05, 0.96080E+05, 0.96814E+05, 0.97552E+05, 0.98294E+05, 0.99041E+05, 0.99793E+05, 0.10055E+06 &
, 0.10131E+06, 0.10208E+06, 0.10285E+06, 0.10362E+06, 0.10440E+06, 0.10518E+06, 0.10597E+06, 0.10676E+06, 0.10756E+06, 0.10836E+06 &
, 0.10917E+06, 0.10998E+06, 0.11080E+06, 0.11162E+06, 0.11245E+06, 0.11328E+06, 0.11412E+06, 0.11496E+06, 0.11580E+06, 0.11665E+06 &
, 0.11751E+06, 0.11837E+06, 0.11924E+06, 0.12011E+06, 0.12098E+06, 0.12186E+06, 0.12275E+06, 0.12364E+06, 0.12454E+06, 0.12544E+06 &
, 0.12635E+06, 0.12726E+06, 0.12817E+06, 0.12910E+06, 0.13002E+06, 0.13096E+06, 0.13190E+06, 0.13284E+06, 0.13379E+06, 0.13474E+06 &
, 0.13570E+06, 0.13667E+06, 0.13764E+06, 0.13861E+06, 0.13960E+06, 0.14058E+06, 0.14158E+06, 0.14257E+06, 0.14358E+06, 0.14459E+06 &
, 0.14560E+06, 0.14662E+06, 0.14765E+06, 0.14868E+06, 0.14972E+06, 0.15076E+06, 0.15181E+06, 0.15287E+06, 0.15393E+06, 0.15500E+06 &
, 0.15607E+06, 0.15715E+06, 0.15823E+06, 0.15932E+06, 0.16042E+06, 0.16152E+06, 0.16263E+06, 0.16375E+06, 0.16487E+06, 0.16599E+06 &
, 0.16713E+06, 0.16827E+06, 0.16941E+06, 0.17056E+06, 0.17172E+06, 0.17288E+06, 0.17405E+06, 0.17523E+06, 0.17641E+06, 0.17760E+06 &
, 0.17880E+06, 0.18000E+06, 0.18121E+06, 0.18242E+06, 0.18365E+06, 0.18487E+06, 0.18611E+06, 0.18735E+06, 0.18860E+06, 0.18985E+06 &
, 0.19111E+06, 0.19238E+06, 0.19365E+06, 0.19494E+06, 0.19622E+06, 0.19752E+06, 0.19882E+06, 0.20013E+06, 0.20144E+06, 0.20276E+06 &
, 0.20409E+06, 0.20543E+06, 0.20677E+06, 0.20812E+06, 0.20948E+06, 0.21084E+06, 0.21221E+06, 0.21359E+06, 0.21497E+06, 0.21636E+06 &
, 0.21776E+06, 0.21917E+06, 0.22058E+06, 0.22200E+06, 0.22343E+06, 0.22487E+06, 0.22631E+06, 0.22776E+06, 0.22922E+06, 0.23068E+06 &
, 0.23215E+06, 0.23363E+06, 0.23512E+06, 0.23662E+06, 0.23812E+06, 0.23963E+06, 0.24115E+06, 0.24267E+06, 0.24420E+06, 0.24574E+06 &
, 0.24729E+06, 0.24885E+06, 0.25041E+06, 0.25198E+06, 0.25356E+06, 0.25515E+06, 0.25675E+06, 0.25835E+06, 0.25996E+06, 0.26158E+06 &
, 0.26321E+06, 0.26484E+06, 0.26648E+06, 0.26814E+06, 0.26979E+06, 0.27146E+06, 0.27314E+06, 0.27482E+06, 0.27651E+06, 0.27821E+06 &
, 0.27992E+06, 0.28164E+06, 0.28337E+06, 0.28510E+06, 0.28684E+06, 0.28859E+06, 0.29035E+06, 0.29212E+06, 0.29390E+06, 0.29568E+06 &
, 0.29748E+06, 0.29928E+06, 0.30109E+06, 0.30291E+06, 0.30474E+06, 0.30657E+06, 0.30842E+06, 0.31027E+06, 0.31214E+06, 0.31401E+06 &
, 0.31589E+06, 0.31778E+06, 0.31968E+06, 0.32159E+06, 0.32351E+06, 0.32543E+06, 0.32737E+06, 0.32931E+06, 0.33127E+06, 0.33323E+06 &
, 0.33520E+06, 0.33718E+06, 0.33917E+06, 0.34117E+06, 0.34318E+06, 0.34520E+06, 0.34723E+06, 0.34927E+06, 0.35132E+06, 0.35337E+06 /
 !            3           6
 DATA(QofT(          28 ,J),J=1,501)/   130.839996337891       &
, 0.15084E+03, 0.17176E+03, 0.19356E+03, 0.21623E+03, 0.23971E+03, 0.26399E+03, 0.28904E+03, 0.31483E+03, 0.34135E+03, 0.36858E+03 &
, 0.39650E+03, 0.42509E+03, 0.45434E+03, 0.48423E+03, 0.51475E+03, 0.54589E+03, 0.57764E+03, 0.60998E+03, 0.64290E+03, 0.67639E+03 &
, 0.71046E+03, 0.74507E+03, 0.78023E+03, 0.81593E+03, 0.85216E+03, 0.88892E+03, 0.92619E+03, 0.96397E+03, 0.10022E+04, 0.10410E+04 &
, 0.10803E+04, 0.11201E+04, 0.11603E+04, 0.12010E+04, 0.12422E+04, 0.12838E+04, 0.13259E+04, 0.13685E+04, 0.14115E+04, 0.14549E+04 &
, 0.14988E+04, 0.15431E+04, 0.15879E+04, 0.16331E+04, 0.16787E+04, 0.17248E+04, 0.17712E+04, 0.18181E+04, 0.18655E+04, 0.19132E+04 &
, 0.19613E+04, 0.20099E+04, 0.20589E+04, 0.21082E+04, 0.21580E+04, 0.22082E+04, 0.22588E+04, 0.23098E+04, 0.23612E+04, 0.24130E+04 &
, 0.24652E+04, 0.25179E+04, 0.25709E+04, 0.26243E+04, 0.26781E+04, 0.27323E+04, 0.27869E+04, 0.28420E+04, 0.28974E+04, 0.29532E+04 &
, 0.30094E+04, 0.30661E+04, 0.31231E+04, 0.31806E+04, 0.32384E+04, 0.32967E+04, 0.33554E+04, 0.34145E+04, 0.34740E+04, 0.35339E+04 &
, 0.35943E+04, 0.36550E+04, 0.37162E+04, 0.37778E+04, 0.38399E+04, 0.39023E+04, 0.39652E+04, 0.40285E+04, 0.40923E+04, 0.41565E+04 &
, 0.42211E+04, 0.42862E+04, 0.43518E+04, 0.44178E+04, 0.44842E+04, 0.45511E+04, 0.46185E+04, 0.46863E+04, 0.47546E+04, 0.48234E+04 &
, 0.48926E+04, 0.49623E+04, 0.50325E+04, 0.51032E+04, 0.51744E+04, 0.52460E+04, 0.53182E+04, 0.53908E+04, 0.54640E+04, 0.55377E+04 &
, 0.56118E+04, 0.56865E+04, 0.57617E+04, 0.58375E+04, 0.59137E+04, 0.59905E+04, 0.60678E+04, 0.61457E+04, 0.62241E+04, 0.63030E+04 &
, 0.63826E+04, 0.64626E+04, 0.65433E+04, 0.66245E+04, 0.67063E+04, 0.67886E+04, 0.68715E+04, 0.69550E+04, 0.70391E+04, 0.71238E+04 &
, 0.72090E+04, 0.72949E+04, 0.73815E+04, 0.74686E+04, 0.75563E+04, 0.76446E+04, 0.77337E+04, 0.78232E+04, 0.79135E+04, 0.80044E+04 &
, 0.80959E+04, 0.81881E+04, 0.82809E+04, 0.83744E+04, 0.84686E+04, 0.85634E+04, 0.86589E+04, 0.87552E+04, 0.88521E+04, 0.89496E+04 &
, 0.90478E+04, 0.91468E+04, 0.92465E+04, 0.93469E+04, 0.94480E+04, 0.95499E+04, 0.96525E+04, 0.97558E+04, 0.98598E+04, 0.99645E+04 &
, 0.10070E+05, 0.10176E+05, 0.10283E+05, 0.10391E+05, 0.10500E+05, 0.10609E+05, 0.10719E+05, 0.10830E+05, 0.10942E+05, 0.11054E+05 &
, 0.11168E+05, 0.11282E+05, 0.11397E+05, 0.11513E+05, 0.11629E+05, 0.11747E+05, 0.11865E+05, 0.11984E+05, 0.12104E+05, 0.12225E+05 &
, 0.12346E+05, 0.12469E+05, 0.12592E+05, 0.12716E+05, 0.12842E+05, 0.12968E+05, 0.13095E+05, 0.13223E+05, 0.13351E+05, 0.13481E+05 &
, 0.13611E+05, 0.13743E+05, 0.13875E+05, 0.14009E+05, 0.14143E+05, 0.14278E+05, 0.14414E+05, 0.14552E+05, 0.14690E+05, 0.14829E+05 &
, 0.14969E+05, 0.15110E+05, 0.15252E+05, 0.15395E+05, 0.15539E+05, 0.15684E+05, 0.15830E+05, 0.15977E+05, 0.16125E+05, 0.16274E+05 &
, 0.16424E+05, 0.16575E+05, 0.16728E+05, 0.16881E+05, 0.17035E+05, 0.17191E+05, 0.17347E+05, 0.17505E+05, 0.17663E+05, 0.17823E+05 &
, 0.17984E+05, 0.18146E+05, 0.18309E+05, 0.18473E+05, 0.18638E+05, 0.18805E+05, 0.18972E+05, 0.19141E+05, 0.19311E+05, 0.19482E+05 &
, 0.19654E+05, 0.19827E+05, 0.20002E+05, 0.20178E+05, 0.20355E+05, 0.20532E+05, 0.20712E+05, 0.20892E+05, 0.21074E+05, 0.21257E+05 &
, 0.21441E+05, 0.21626E+05, 0.21813E+05, 0.22001E+05, 0.22190E+05, 0.22380E+05, 0.22572E+05, 0.22765E+05, 0.22959E+05, 0.23155E+05 &
, 0.23351E+05, 0.23550E+05, 0.23749E+05, 0.23950E+05, 0.24152E+05, 0.24355E+05, 0.24560E+05, 0.24766E+05, 0.24973E+05, 0.25182E+05 &
, 0.25392E+05, 0.25604E+05, 0.25817E+05, 0.26031E+05, 0.26247E+05, 0.26464E+05, 0.26682E+05, 0.26902E+05, 0.27123E+05, 0.27346E+05 &
, 0.27570E+05, 0.27796E+05, 0.28023E+05, 0.28251E+05, 0.28481E+05, 0.28713E+05, 0.28946E+05, 0.29180E+05, 0.29416E+05, 0.29653E+05 &
, 0.29892E+05, 0.30133E+05, 0.30375E+05, 0.30618E+05, 0.30863E+05, 0.31110E+05, 0.31358E+05, 0.31608E+05, 0.31859E+05, 0.32112E+05 &
, 0.32366E+05, 0.32622E+05, 0.32880E+05, 0.33139E+05, 0.33400E+05, 0.33662E+05, 0.33927E+05, 0.34192E+05, 0.34459E+05, 0.34729E+05 &
, 0.34999E+05, 0.35272E+05, 0.35546E+05, 0.35821E+05, 0.36099E+05, 0.36378E+05, 0.36659E+05, 0.36941E+05, 0.37226E+05, 0.37512E+05 &
, 0.37799E+05, 0.38089E+05, 0.38380E+05, 0.38673E+05, 0.38968E+05, 0.39265E+05, 0.39563E+05, 0.39863E+05, 0.40165E+05, 0.40469E+05 &
, 0.40774E+05, 0.41082E+05, 0.41391E+05, 0.41702E+05, 0.42015E+05, 0.42330E+05, 0.42647E+05, 0.42965E+05, 0.43286E+05, 0.43608E+05 &
, 0.43933E+05, 0.44259E+05, 0.44587E+05, 0.44917E+05, 0.45249E+05, 0.45583E+05, 0.45919E+05, 0.46257E+05, 0.46596E+05, 0.46938E+05 &
, 0.47282E+05, 0.47628E+05, 0.47976E+05, 0.48326E+05, 0.48678E+05, 0.49032E+05, 0.49388E+05, 0.49746E+05, 0.50106E+05, 0.50468E+05 &
, 0.50832E+05, 0.51198E+05, 0.51567E+05, 0.51937E+05, 0.52310E+05, 0.52685E+05, 0.53062E+05, 0.53441E+05, 0.53822E+05, 0.54206E+05 &
, 0.54591E+05, 0.54979E+05, 0.55369E+05, 0.55761E+05, 0.56155E+05, 0.56551E+05, 0.56950E+05, 0.57351E+05, 0.57755E+05, 0.58160E+05 &
, 0.58568E+05, 0.58978E+05, 0.59390E+05, 0.59805E+05, 0.60222E+05, 0.60641E+05, 0.61063E+05, 0.61486E+05, 0.61913E+05, 0.62341E+05 &
, 0.62772E+05, 0.63206E+05, 0.63641E+05, 0.64079E+05, 0.64520E+05, 0.64963E+05, 0.65408E+05, 0.65856E+05, 0.66306E+05, 0.66758E+05 &
, 0.67213E+05, 0.67671E+05, 0.68131E+05, 0.68593E+05, 0.69058E+05, 0.69526E+05, 0.69996E+05, 0.70469E+05, 0.70944E+05, 0.71421E+05 &
, 0.71901E+05, 0.72384E+05, 0.72870E+05, 0.73357E+05, 0.73848E+05, 0.74341E+05, 0.74837E+05, 0.75335E+05, 0.75836E+05, 0.76340E+05 &
, 0.76846E+05, 0.77355E+05, 0.77867E+05, 0.78381E+05, 0.78898E+05, 0.79418E+05, 0.79940E+05, 0.80466E+05, 0.80994E+05, 0.81524E+05 &
, 0.82058E+05, 0.82594E+05, 0.83133E+05, 0.83675E+05, 0.84219E+05, 0.84767E+05, 0.85317E+05, 0.85871E+05, 0.86427E+05, 0.86985E+05 &
, 0.87547E+05, 0.88112E+05, 0.88679E+05, 0.89249E+05, 0.89823E+05, 0.90399E+05, 0.90978E+05, 0.91561E+05, 0.92145E+05, 0.92733E+05 &
, 0.93324E+05, 0.93919E+05, 0.94516E+05, 0.95116E+05, 0.95719E+05, 0.96325E+05, 0.96934E+05, 0.97546E+05, 0.98161E+05, 0.98780E+05 &
, 0.99401E+05, 0.10003E+06, 0.10065E+06, 0.10128E+06, 0.10192E+06, 0.10256E+06, 0.10320E+06, 0.10384E+06, 0.10449E+06, 0.10514E+06 &
, 0.10579E+06, 0.10644E+06, 0.10710E+06, 0.10777E+06, 0.10843E+06, 0.10910E+06, 0.10978E+06, 0.11045E+06, 0.11113E+06, 0.11181E+06 &
, 0.11250E+06, 0.11319E+06, 0.11388E+06, 0.11458E+06, 0.11527E+06, 0.11598E+06, 0.11668E+06, 0.11739E+06, 0.11811E+06, 0.11882E+06 &
, 0.11954E+06, 0.12026E+06, 0.12099E+06, 0.12172E+06, 0.12245E+06, 0.12319E+06, 0.12393E+06, 0.12468E+06, 0.12542E+06, 0.12618E+06 &
, 0.12693E+06, 0.12769E+06, 0.12845E+06, 0.12922E+06, 0.12999E+06, 0.13076E+06, 0.13154E+06, 0.13232E+06, 0.13310E+06, 0.13389E+06 /
 !            3           7
 DATA(QofT(          29 ,J),J=1,501)/   66.9420013427734       &
, 0.77173E+02, 0.87879E+02, 0.99040E+02, 0.11064E+03, 0.12266E+03, 0.13508E+03, 0.14790E+03, 0.16110E+03, 0.17467E+03, 0.18861E+03 &
, 0.20290E+03, 0.21753E+03, 0.23250E+03, 0.24780E+03, 0.26342E+03, 0.27936E+03, 0.29560E+03, 0.31215E+03, 0.32900E+03, 0.34615E+03 &
, 0.36358E+03, 0.38130E+03, 0.39929E+03, 0.41756E+03, 0.43611E+03, 0.45492E+03, 0.47399E+03, 0.49333E+03, 0.51292E+03, 0.53277E+03 &
, 0.55287E+03, 0.57321E+03, 0.59380E+03, 0.61464E+03, 0.63571E+03, 0.65702E+03, 0.67857E+03, 0.70035E+03, 0.72236E+03, 0.74460E+03 &
, 0.76706E+03, 0.78976E+03, 0.81267E+03, 0.83580E+03, 0.85915E+03, 0.88273E+03, 0.90652E+03, 0.93053E+03, 0.95475E+03, 0.97918E+03 &
, 0.10038E+04, 0.10287E+04, 0.10538E+04, 0.10790E+04, 0.11045E+04, 0.11302E+04, 0.11561E+04, 0.11822E+04, 0.12086E+04, 0.12351E+04 &
, 0.12618E+04, 0.12888E+04, 0.13159E+04, 0.13433E+04, 0.13708E+04, 0.13986E+04, 0.14265E+04, 0.14547E+04, 0.14831E+04, 0.15117E+04 &
, 0.15405E+04, 0.15695E+04, 0.15987E+04, 0.16281E+04, 0.16578E+04, 0.16876E+04, 0.17177E+04, 0.17480E+04, 0.17784E+04, 0.18091E+04 &
, 0.18400E+04, 0.18712E+04, 0.19025E+04, 0.19341E+04, 0.19659E+04, 0.19979E+04, 0.20301E+04, 0.20625E+04, 0.20952E+04, 0.21281E+04 &
, 0.21612E+04, 0.21946E+04, 0.22281E+04, 0.22620E+04, 0.22960E+04, 0.23303E+04, 0.23648E+04, 0.23996E+04, 0.24346E+04, 0.24698E+04 &
, 0.25053E+04, 0.25410E+04, 0.25769E+04, 0.26131E+04, 0.26496E+04, 0.26863E+04, 0.27233E+04, 0.27605E+04, 0.27980E+04, 0.28357E+04 &
, 0.28737E+04, 0.29120E+04, 0.29505E+04, 0.29893E+04, 0.30284E+04, 0.30677E+04, 0.31073E+04, 0.31472E+04, 0.31874E+04, 0.32278E+04 &
, 0.32685E+04, 0.33095E+04, 0.33508E+04, 0.33924E+04, 0.34343E+04, 0.34765E+04, 0.35189E+04, 0.35617E+04, 0.36048E+04, 0.36481E+04 &
, 0.36918E+04, 0.37358E+04, 0.37801E+04, 0.38247E+04, 0.38696E+04, 0.39148E+04, 0.39604E+04, 0.40062E+04, 0.40524E+04, 0.40990E+04 &
, 0.41458E+04, 0.41930E+04, 0.42405E+04, 0.42884E+04, 0.43366E+04, 0.43851E+04, 0.44340E+04, 0.44832E+04, 0.45328E+04, 0.45827E+04 &
, 0.46330E+04, 0.46836E+04, 0.47346E+04, 0.47860E+04, 0.48377E+04, 0.48898E+04, 0.49422E+04, 0.49951E+04, 0.50483E+04, 0.51019E+04 &
, 0.51558E+04, 0.52102E+04, 0.52649E+04, 0.53200E+04, 0.53756E+04, 0.54315E+04, 0.54878E+04, 0.55445E+04, 0.56016E+04, 0.56591E+04 &
, 0.57171E+04, 0.57754E+04, 0.58341E+04, 0.58933E+04, 0.59529E+04, 0.60129E+04, 0.60733E+04, 0.61342E+04, 0.61955E+04, 0.62572E+04 &
, 0.63193E+04, 0.63819E+04, 0.64450E+04, 0.65085E+04, 0.65724E+04, 0.66367E+04, 0.67016E+04, 0.67669E+04, 0.68326E+04, 0.68988E+04 &
, 0.69655E+04, 0.70326E+04, 0.71002E+04, 0.71683E+04, 0.72369E+04, 0.73060E+04, 0.73755E+04, 0.74455E+04, 0.75160E+04, 0.75870E+04 &
, 0.76585E+04, 0.77304E+04, 0.78029E+04, 0.78759E+04, 0.79494E+04, 0.80235E+04, 0.80980E+04, 0.81730E+04, 0.82485E+04, 0.83246E+04 &
, 0.84012E+04, 0.84783E+04, 0.85560E+04, 0.86342E+04, 0.87129E+04, 0.87922E+04, 0.88720E+04, 0.89523E+04, 0.90333E+04, 0.91148E+04 &
, 0.91968E+04, 0.92793E+04, 0.93625E+04, 0.94462E+04, 0.95305E+04, 0.96154E+04, 0.97008E+04, 0.97868E+04, 0.98734E+04, 0.99605E+04 &
, 0.10048E+05, 0.10137E+05, 0.10226E+05, 0.10315E+05, 0.10405E+05, 0.10496E+05, 0.10587E+05, 0.10679E+05, 0.10772E+05, 0.10865E+05 &
, 0.10959E+05, 0.11054E+05, 0.11149E+05, 0.11244E+05, 0.11341E+05, 0.11438E+05, 0.11535E+05, 0.11634E+05, 0.11733E+05, 0.11832E+05 &
, 0.11932E+05, 0.12033E+05, 0.12135E+05, 0.12237E+05, 0.12340E+05, 0.12443E+05, 0.12548E+05, 0.12653E+05, 0.12758E+05, 0.12865E+05 &
, 0.12972E+05, 0.13079E+05, 0.13188E+05, 0.13297E+05, 0.13407E+05, 0.13517E+05, 0.13628E+05, 0.13740E+05, 0.13853E+05, 0.13966E+05 &
, 0.14080E+05, 0.14195E+05, 0.14311E+05, 0.14427E+05, 0.14544E+05, 0.14662E+05, 0.14780E+05, 0.14900E+05, 0.15020E+05, 0.15140E+05 &
, 0.15262E+05, 0.15384E+05, 0.15507E+05, 0.15631E+05, 0.15756E+05, 0.15881E+05, 0.16008E+05, 0.16135E+05, 0.16263E+05, 0.16391E+05 &
, 0.16521E+05, 0.16651E+05, 0.16782E+05, 0.16914E+05, 0.17046E+05, 0.17180E+05, 0.17314E+05, 0.17449E+05, 0.17585E+05, 0.17722E+05 &
, 0.17860E+05, 0.17998E+05, 0.18137E+05, 0.18278E+05, 0.18419E+05, 0.18561E+05, 0.18703E+05, 0.18847E+05, 0.18992E+05, 0.19137E+05 &
, 0.19283E+05, 0.19430E+05, 0.19578E+05, 0.19727E+05, 0.19877E+05, 0.20028E+05, 0.20180E+05, 0.20332E+05, 0.20486E+05, 0.20640E+05 &
, 0.20795E+05, 0.20952E+05, 0.21109E+05, 0.21267E+05, 0.21426E+05, 0.21586E+05, 0.21747E+05, 0.21909E+05, 0.22071E+05, 0.22235E+05 &
, 0.22400E+05, 0.22566E+05, 0.22732E+05, 0.22900E+05, 0.23069E+05, 0.23238E+05, 0.23409E+05, 0.23581E+05, 0.23753E+05, 0.23927E+05 &
, 0.24101E+05, 0.24277E+05, 0.24454E+05, 0.24631E+05, 0.24810E+05, 0.24990E+05, 0.25171E+05, 0.25352E+05, 0.25535E+05, 0.25719E+05 &
, 0.25904E+05, 0.26090E+05, 0.26277E+05, 0.26465E+05, 0.26655E+05, 0.26845E+05, 0.27036E+05, 0.27229E+05, 0.27422E+05, 0.27617E+05 &
, 0.27813E+05, 0.28009E+05, 0.28207E+05, 0.28406E+05, 0.28607E+05, 0.28808E+05, 0.29010E+05, 0.29214E+05, 0.29418E+05, 0.29624E+05 &
, 0.29831E+05, 0.30039E+05, 0.30249E+05, 0.30459E+05, 0.30671E+05, 0.30883E+05, 0.31097E+05, 0.31312E+05, 0.31529E+05, 0.31746E+05 &
, 0.31965E+05, 0.32185E+05, 0.32406E+05, 0.32628E+05, 0.32851E+05, 0.33076E+05, 0.33302E+05, 0.33529E+05, 0.33758E+05, 0.33987E+05 &
, 0.34218E+05, 0.34450E+05, 0.34683E+05, 0.34918E+05, 0.35154E+05, 0.35391E+05, 0.35630E+05, 0.35869E+05, 0.36110E+05, 0.36352E+05 &
, 0.36596E+05, 0.36841E+05, 0.37087E+05, 0.37334E+05, 0.37583E+05, 0.37833E+05, 0.38084E+05, 0.38337E+05, 0.38591E+05, 0.38847E+05 &
, 0.39103E+05, 0.39362E+05, 0.39621E+05, 0.39882E+05, 0.40144E+05, 0.40407E+05, 0.40672E+05, 0.40939E+05, 0.41206E+05, 0.41475E+05 &
, 0.41746E+05, 0.42018E+05, 0.42291E+05, 0.42566E+05, 0.42842E+05, 0.43119E+05, 0.43398E+05, 0.43679E+05, 0.43960E+05, 0.44244E+05 &
, 0.44529E+05, 0.44815E+05, 0.45102E+05, 0.45391E+05, 0.45682E+05, 0.45974E+05, 0.46268E+05, 0.46563E+05, 0.46859E+05, 0.47157E+05 &
, 0.47457E+05, 0.47758E+05, 0.48060E+05, 0.48364E+05, 0.48670E+05, 0.48977E+05, 0.49286E+05, 0.49596E+05, 0.49908E+05, 0.50221E+05 &
, 0.50536E+05, 0.50852E+05, 0.51170E+05, 0.51490E+05, 0.51811E+05, 0.52134E+05, 0.52458E+05, 0.52784E+05, 0.53112E+05, 0.53441E+05 &
, 0.53772E+05, 0.54104E+05, 0.54439E+05, 0.54774E+05, 0.55112E+05, 0.55451E+05, 0.55791E+05, 0.56134E+05, 0.56478E+05, 0.56823E+05 &
, 0.57171E+05, 0.57520E+05, 0.57870E+05, 0.58223E+05, 0.58577E+05, 0.58933E+05, 0.59290E+05, 0.59650E+05, 0.60011E+05, 0.60373E+05 &
, 0.60738E+05, 0.61104E+05, 0.61472E+05, 0.61842E+05, 0.62213E+05, 0.62587E+05, 0.62962E+05, 0.63338E+05, 0.63717E+05, 0.64098E+05 &
, 0.64480E+05, 0.64864E+05, 0.65250E+05, 0.65637E+05, 0.66026E+05, 0.66418E+05, 0.66811E+05, 0.67206E+05, 0.67603E+05, 0.68002E+05 /
 !            3           8
 DATA(QofT(          30 ,J),J=1,501)/   768.869995117188       &
, 0.88635E+03, 0.10093E+04, 0.11374E+04, 0.12706E+04, 0.14086E+04, 0.15512E+04, 0.16984E+04, 0.18500E+04, 0.20058E+04, 0.21658E+04 &
, 0.23299E+04, 0.24979E+04, 0.26697E+04, 0.28454E+04, 0.30247E+04, 0.32077E+04, 0.33942E+04, 0.35843E+04, 0.37777E+04, 0.39745E+04 &
, 0.41747E+04, 0.43781E+04, 0.45847E+04, 0.47945E+04, 0.50074E+04, 0.52233E+04, 0.54423E+04, 0.56643E+04, 0.58893E+04, 0.61171E+04 &
, 0.63479E+04, 0.65815E+04, 0.68179E+04, 0.70570E+04, 0.72990E+04, 0.75436E+04, 0.77910E+04, 0.80411E+04, 0.82938E+04, 0.85491E+04 &
, 0.88070E+04, 0.90674E+04, 0.93305E+04, 0.95961E+04, 0.98642E+04, 0.10135E+05, 0.10408E+05, 0.10683E+05, 0.10961E+05, 0.11242E+05 &
, 0.11525E+05, 0.11810E+05, 0.12098E+05, 0.12388E+05, 0.12680E+05, 0.12975E+05, 0.13273E+05, 0.13572E+05, 0.13874E+05, 0.14178E+05 &
, 0.14485E+05, 0.14794E+05, 0.15106E+05, 0.15419E+05, 0.15736E+05, 0.16054E+05, 0.16375E+05, 0.16698E+05, 0.17024E+05, 0.17352E+05 &
, 0.17682E+05, 0.18014E+05, 0.18350E+05, 0.18687E+05, 0.19027E+05, 0.19369E+05, 0.19713E+05, 0.20060E+05, 0.20410E+05, 0.20761E+05 &
, 0.21116E+05, 0.21472E+05, 0.21832E+05, 0.22193E+05, 0.22558E+05, 0.22924E+05, 0.23293E+05, 0.23665E+05, 0.24039E+05, 0.24416E+05 &
, 0.24795E+05, 0.25177E+05, 0.25562E+05, 0.25949E+05, 0.26339E+05, 0.26731E+05, 0.27126E+05, 0.27524E+05, 0.27925E+05, 0.28328E+05 &
, 0.28734E+05, 0.29143E+05, 0.29555E+05, 0.29969E+05, 0.30386E+05, 0.30806E+05, 0.31229E+05, 0.31655E+05, 0.32084E+05, 0.32516E+05 &
, 0.32950E+05, 0.33388E+05, 0.33829E+05, 0.34273E+05, 0.34719E+05, 0.35169E+05, 0.35622E+05, 0.36078E+05, 0.36537E+05, 0.37000E+05 &
, 0.37466E+05, 0.37934E+05, 0.38406E+05, 0.38882E+05, 0.39361E+05, 0.39842E+05, 0.40328E+05, 0.40817E+05, 0.41309E+05, 0.41804E+05 &
, 0.42303E+05, 0.42806E+05, 0.43312E+05, 0.43821E+05, 0.44335E+05, 0.44852E+05, 0.45372E+05, 0.45896E+05, 0.46423E+05, 0.46955E+05 &
, 0.47490E+05, 0.48029E+05, 0.48572E+05, 0.49118E+05, 0.49669E+05, 0.50223E+05, 0.50781E+05, 0.51343E+05, 0.51909E+05, 0.52479E+05 &
, 0.53053E+05, 0.53632E+05, 0.54214E+05, 0.54800E+05, 0.55391E+05, 0.55985E+05, 0.56585E+05, 0.57188E+05, 0.57795E+05, 0.58407E+05 &
, 0.59023E+05, 0.59643E+05, 0.60268E+05, 0.60897E+05, 0.61530E+05, 0.62169E+05, 0.62811E+05, 0.63458E+05, 0.64110E+05, 0.64766E+05 &
, 0.65427E+05, 0.66093E+05, 0.66763E+05, 0.67438E+05, 0.68119E+05, 0.68803E+05, 0.69492E+05, 0.70187E+05, 0.70886E+05, 0.71590E+05 &
, 0.72299E+05, 0.73013E+05, 0.73732E+05, 0.74456E+05, 0.75186E+05, 0.75920E+05, 0.76660E+05, 0.77405E+05, 0.78154E+05, 0.78909E+05 &
, 0.79670E+05, 0.80436E+05, 0.81207E+05, 0.81983E+05, 0.82765E+05, 0.83553E+05, 0.84345E+05, 0.85144E+05, 0.85948E+05, 0.86758E+05 &
, 0.87573E+05, 0.88394E+05, 0.89220E+05, 0.90052E+05, 0.90891E+05, 0.91734E+05, 0.92584E+05, 0.93439E+05, 0.94301E+05, 0.95169E+05 &
, 0.96041E+05, 0.96921E+05, 0.97807E+05, 0.98698E+05, 0.99596E+05, 0.10050E+06, 0.10141E+06, 0.10233E+06, 0.10325E+06, 0.10418E+06 &
, 0.10511E+06, 0.10605E+06, 0.10700E+06, 0.10795E+06, 0.10892E+06, 0.10988E+06, 0.11086E+06, 0.11184E+06, 0.11282E+06, 0.11382E+06 &
, 0.11482E+06, 0.11583E+06, 0.11684E+06, 0.11786E+06, 0.11889E+06, 0.11992E+06, 0.12096E+06, 0.12201E+06, 0.12307E+06, 0.12413E+06 &
, 0.12520E+06, 0.12627E+06, 0.12736E+06, 0.12845E+06, 0.12955E+06, 0.13065E+06, 0.13176E+06, 0.13288E+06, 0.13401E+06, 0.13515E+06 &
, 0.13629E+06, 0.13744E+06, 0.13860E+06, 0.13976E+06, 0.14093E+06, 0.14211E+06, 0.14330E+06, 0.14450E+06, 0.14570E+06, 0.14691E+06 &
, 0.14813E+06, 0.14936E+06, 0.15059E+06, 0.15183E+06, 0.15308E+06, 0.15434E+06, 0.15561E+06, 0.15689E+06, 0.15817E+06, 0.15946E+06 &
, 0.16076E+06, 0.16207E+06, 0.16339E+06, 0.16471E+06, 0.16604E+06, 0.16739E+06, 0.16874E+06, 0.17009E+06, 0.17146E+06, 0.17284E+06 &
, 0.17422E+06, 0.17562E+06, 0.17702E+06, 0.17843E+06, 0.17985E+06, 0.18128E+06, 0.18272E+06, 0.18416E+06, 0.18562E+06, 0.18708E+06 &
, 0.18856E+06, 0.19004E+06, 0.19153E+06, 0.19304E+06, 0.19455E+06, 0.19607E+06, 0.19760E+06, 0.19914E+06, 0.20068E+06, 0.20224E+06 &
, 0.20381E+06, 0.20539E+06, 0.20697E+06, 0.20857E+06, 0.21018E+06, 0.21179E+06, 0.21342E+06, 0.21505E+06, 0.21670E+06, 0.21836E+06 &
, 0.22002E+06, 0.22170E+06, 0.22338E+06, 0.22508E+06, 0.22679E+06, 0.22850E+06, 0.23023E+06, 0.23197E+06, 0.23371E+06, 0.23547E+06 &
, 0.23724E+06, 0.23902E+06, 0.24081E+06, 0.24261E+06, 0.24442E+06, 0.24624E+06, 0.24808E+06, 0.24992E+06, 0.25177E+06, 0.25364E+06 &
, 0.25551E+06, 0.25740E+06, 0.25930E+06, 0.26121E+06, 0.26313E+06, 0.26506E+06, 0.26700E+06, 0.26896E+06, 0.27092E+06, 0.27290E+06 &
, 0.27489E+06, 0.27689E+06, 0.27890E+06, 0.28092E+06, 0.28296E+06, 0.28500E+06, 0.28706E+06, 0.28913E+06, 0.29121E+06, 0.29331E+06 &
, 0.29541E+06, 0.29753E+06, 0.29966E+06, 0.30180E+06, 0.30395E+06, 0.30612E+06, 0.30830E+06, 0.31049E+06, 0.31269E+06, 0.31491E+06 &
, 0.31714E+06, 0.31938E+06, 0.32163E+06, 0.32390E+06, 0.32617E+06, 0.32847E+06, 0.33077E+06, 0.33309E+06, 0.33542E+06, 0.33776E+06 &
, 0.34011E+06, 0.34248E+06, 0.34487E+06, 0.34726E+06, 0.34967E+06, 0.35209E+06, 0.35453E+06, 0.35697E+06, 0.35944E+06, 0.36191E+06 &
, 0.36440E+06, 0.36690E+06, 0.36942E+06, 0.37195E+06, 0.37449E+06, 0.37705E+06, 0.37962E+06, 0.38221E+06, 0.38481E+06, 0.38742E+06 &
, 0.39005E+06, 0.39269E+06, 0.39535E+06, 0.39801E+06, 0.40070E+06, 0.40340E+06, 0.40611E+06, 0.40884E+06, 0.41158E+06, 0.41434E+06 &
, 0.41711E+06, 0.41990E+06, 0.42270E+06, 0.42551E+06, 0.42835E+06, 0.43119E+06, 0.43405E+06, 0.43693E+06, 0.43982E+06, 0.44273E+06 &
, 0.44565E+06, 0.44859E+06, 0.45154E+06, 0.45451E+06, 0.45749E+06, 0.46049E+06, 0.46350E+06, 0.46654E+06, 0.46958E+06, 0.47264E+06 &
, 0.47572E+06, 0.47882E+06, 0.48192E+06, 0.48505E+06, 0.48819E+06, 0.49135E+06, 0.49452E+06, 0.49771E+06, 0.50092E+06, 0.50415E+06 &
, 0.50739E+06, 0.51064E+06, 0.51392E+06, 0.51721E+06, 0.52051E+06, 0.52384E+06, 0.52718E+06, 0.53053E+06, 0.53391E+06, 0.53730E+06 &
, 0.54071E+06, 0.54413E+06, 0.54757E+06, 0.55103E+06, 0.55451E+06, 0.55801E+06, 0.56152E+06, 0.56505E+06, 0.56860E+06, 0.57216E+06 &
, 0.57574E+06, 0.57934E+06, 0.58296E+06, 0.58660E+06, 0.59025E+06, 0.59393E+06, 0.59762E+06, 0.60133E+06, 0.60505E+06, 0.60880E+06 &
, 0.61256E+06, 0.61635E+06, 0.62015E+06, 0.62397E+06, 0.62780E+06, 0.63166E+06, 0.63554E+06, 0.63943E+06, 0.64335E+06, 0.64728E+06 &
, 0.65123E+06, 0.65520E+06, 0.65919E+06, 0.66320E+06, 0.66723E+06, 0.67128E+06, 0.67535E+06, 0.67943E+06, 0.68354E+06, 0.68767E+06 &
, 0.69182E+06, 0.69598E+06, 0.70017E+06, 0.70437E+06, 0.70860E+06, 0.71285E+06, 0.71711E+06, 0.72140E+06, 0.72571E+06, 0.73004E+06 &
, 0.73438E+06, 0.73875E+06, 0.74314E+06, 0.74755E+06, 0.75198E+06, 0.75643E+06, 0.76091E+06, 0.76540E+06, 0.76991E+06, 0.77445E+06 /
 !            3           9
 DATA(QofT(          31 ,J),J=1,501)/   778.090026855469       &
, 0.89699E+03, 0.10214E+04, 0.11511E+04, 0.12859E+04, 0.14256E+04, 0.15699E+04, 0.17189E+04, 0.18723E+04, 0.20301E+04, 0.21920E+04 &
, 0.23581E+04, 0.25281E+04, 0.27021E+04, 0.28799E+04, 0.30614E+04, 0.32466E+04, 0.34354E+04, 0.36277E+04, 0.38235E+04, 0.40228E+04 &
, 0.42253E+04, 0.44312E+04, 0.46404E+04, 0.48527E+04, 0.50682E+04, 0.52868E+04, 0.55084E+04, 0.57331E+04, 0.59608E+04, 0.61915E+04 &
, 0.64251E+04, 0.66615E+04, 0.69007E+04, 0.71429E+04, 0.73877E+04, 0.76354E+04, 0.78858E+04, 0.81389E+04, 0.83946E+04, 0.86530E+04 &
, 0.89141E+04, 0.91778E+04, 0.94441E+04, 0.97128E+04, 0.99842E+04, 0.10258E+05, 0.10535E+05, 0.10813E+05, 0.11095E+05, 0.11379E+05 &
, 0.11665E+05, 0.11954E+05, 0.12245E+05, 0.12539E+05, 0.12835E+05, 0.13134E+05, 0.13435E+05, 0.13738E+05, 0.14044E+05, 0.14352E+05 &
, 0.14662E+05, 0.14975E+05, 0.15290E+05, 0.15608E+05, 0.15928E+05, 0.16251E+05, 0.16575E+05, 0.16903E+05, 0.17232E+05, 0.17564E+05 &
, 0.17899E+05, 0.18235E+05, 0.18575E+05, 0.18916E+05, 0.19260E+05, 0.19607E+05, 0.19956E+05, 0.20307E+05, 0.20661E+05, 0.21017E+05 &
, 0.21376E+05, 0.21737E+05, 0.22101E+05, 0.22467E+05, 0.22836E+05, 0.23207E+05, 0.23581E+05, 0.23957E+05, 0.24336E+05, 0.24718E+05 &
, 0.25102E+05, 0.25489E+05, 0.25879E+05, 0.26271E+05, 0.26665E+05, 0.27063E+05, 0.27463E+05, 0.27866E+05, 0.28272E+05, 0.28680E+05 &
, 0.29091E+05, 0.29505E+05, 0.29922E+05, 0.30342E+05, 0.30764E+05, 0.31190E+05, 0.31618E+05, 0.32049E+05, 0.32484E+05, 0.32921E+05 &
, 0.33361E+05, 0.33804E+05, 0.34251E+05, 0.34700E+05, 0.35152E+05, 0.35608E+05, 0.36067E+05, 0.36528E+05, 0.36993E+05, 0.37461E+05 &
, 0.37933E+05, 0.38408E+05, 0.38885E+05, 0.39367E+05, 0.39851E+05, 0.40340E+05, 0.40831E+05, 0.41326E+05, 0.41824E+05, 0.42326E+05 &
, 0.42831E+05, 0.43340E+05, 0.43852E+05, 0.44368E+05, 0.44887E+05, 0.45410E+05, 0.45937E+05, 0.46468E+05, 0.47002E+05, 0.47540E+05 &
, 0.48081E+05, 0.48627E+05, 0.49176E+05, 0.49729E+05, 0.50286E+05, 0.50847E+05, 0.51412E+05, 0.51981E+05, 0.52554E+05, 0.53131E+05 &
, 0.53712E+05, 0.54297E+05, 0.54886E+05, 0.55480E+05, 0.56077E+05, 0.56679E+05, 0.57284E+05, 0.57895E+05, 0.58509E+05, 0.59128E+05 &
, 0.59751E+05, 0.60379E+05, 0.61011E+05, 0.61648E+05, 0.62289E+05, 0.62934E+05, 0.63584E+05, 0.64239E+05, 0.64898E+05, 0.65562E+05 &
, 0.66230E+05, 0.66903E+05, 0.67581E+05, 0.68264E+05, 0.68952E+05, 0.69644E+05, 0.70341E+05, 0.71043E+05, 0.71750E+05, 0.72462E+05 &
, 0.73179E+05, 0.73901E+05, 0.74629E+05, 0.75361E+05, 0.76098E+05, 0.76841E+05, 0.77588E+05, 0.78341E+05, 0.79100E+05, 0.79863E+05 &
, 0.80632E+05, 0.81405E+05, 0.82185E+05, 0.82970E+05, 0.83761E+05, 0.84556E+05, 0.85358E+05, 0.86165E+05, 0.86977E+05, 0.87796E+05 &
, 0.88619E+05, 0.89449E+05, 0.90285E+05, 0.91126E+05, 0.91973E+05, 0.92825E+05, 0.93684E+05, 0.94548E+05, 0.95419E+05, 0.96295E+05 &
, 0.97177E+05, 0.98066E+05, 0.98961E+05, 0.99861E+05, 0.10077E+06, 0.10168E+06, 0.10260E+06, 0.10353E+06, 0.10446E+06, 0.10540E+06 &
, 0.10634E+06, 0.10729E+06, 0.10825E+06, 0.10921E+06, 0.11018E+06, 0.11116E+06, 0.11214E+06, 0.11313E+06, 0.11413E+06, 0.11513E+06 &
, 0.11614E+06, 0.11716E+06, 0.11819E+06, 0.11922E+06, 0.12025E+06, 0.12130E+06, 0.12235E+06, 0.12341E+06, 0.12447E+06, 0.12555E+06 &
, 0.12663E+06, 0.12771E+06, 0.12881E+06, 0.12991E+06, 0.13102E+06, 0.13213E+06, 0.13326E+06, 0.13439E+06, 0.13552E+06, 0.13667E+06 &
, 0.13782E+06, 0.13898E+06, 0.14015E+06, 0.14133E+06, 0.14251E+06, 0.14370E+06, 0.14490E+06, 0.14611E+06, 0.14732E+06, 0.14854E+06 &
, 0.14978E+06, 0.15101E+06, 0.15226E+06, 0.15351E+06, 0.15478E+06, 0.15605E+06, 0.15733E+06, 0.15861E+06, 0.15991E+06, 0.16121E+06 &
, 0.16252E+06, 0.16384E+06, 0.16517E+06, 0.16651E+06, 0.16785E+06, 0.16921E+06, 0.17057E+06, 0.17194E+06, 0.17332E+06, 0.17471E+06 &
, 0.17611E+06, 0.17751E+06, 0.17893E+06, 0.18035E+06, 0.18178E+06, 0.18322E+06, 0.18467E+06, 0.18613E+06, 0.18760E+06, 0.18908E+06 &
, 0.19057E+06, 0.19206E+06, 0.19357E+06, 0.19508E+06, 0.19661E+06, 0.19814E+06, 0.19968E+06, 0.20124E+06, 0.20280E+06, 0.20437E+06 &
, 0.20595E+06, 0.20754E+06, 0.20914E+06, 0.21075E+06, 0.21237E+06, 0.21400E+06, 0.21564E+06, 0.21729E+06, 0.21895E+06, 0.22062E+06 &
, 0.22230E+06, 0.22399E+06, 0.22569E+06, 0.22740E+06, 0.22912E+06, 0.23085E+06, 0.23259E+06, 0.23434E+06, 0.23611E+06, 0.23788E+06 &
, 0.23966E+06, 0.24146E+06, 0.24326E+06, 0.24507E+06, 0.24690E+06, 0.24874E+06, 0.25058E+06, 0.25244E+06, 0.25431E+06, 0.25619E+06 &
, 0.25808E+06, 0.25998E+06, 0.26190E+06, 0.26382E+06, 0.26576E+06, 0.26771E+06, 0.26966E+06, 0.27163E+06, 0.27361E+06, 0.27561E+06 &
, 0.27761E+06, 0.27963E+06, 0.28165E+06, 0.28369E+06, 0.28574E+06, 0.28781E+06, 0.28988E+06, 0.29197E+06, 0.29406E+06, 0.29617E+06 &
, 0.29830E+06, 0.30043E+06, 0.30258E+06, 0.30473E+06, 0.30691E+06, 0.30909E+06, 0.31128E+06, 0.31349E+06, 0.31571E+06, 0.31794E+06 &
, 0.32019E+06, 0.32245E+06, 0.32472E+06, 0.32700E+06, 0.32930E+06, 0.33161E+06, 0.33393E+06, 0.33626E+06, 0.33861E+06, 0.34097E+06 &
, 0.34334E+06, 0.34573E+06, 0.34813E+06, 0.35054E+06, 0.35297E+06, 0.35541E+06, 0.35786E+06, 0.36033E+06, 0.36281E+06, 0.36530E+06 &
, 0.36781E+06, 0.37033E+06, 0.37286E+06, 0.37541E+06, 0.37797E+06, 0.38055E+06, 0.38314E+06, 0.38575E+06, 0.38836E+06, 0.39100E+06 &
, 0.39364E+06, 0.39631E+06, 0.39898E+06, 0.40167E+06, 0.40437E+06, 0.40709E+06, 0.40982E+06, 0.41257E+06, 0.41533E+06, 0.41811E+06 &
, 0.42090E+06, 0.42371E+06, 0.42653E+06, 0.42937E+06, 0.43222E+06, 0.43508E+06, 0.43797E+06, 0.44086E+06, 0.44377E+06, 0.44670E+06 &
, 0.44964E+06, 0.45260E+06, 0.45557E+06, 0.45856E+06, 0.46157E+06, 0.46459E+06, 0.46762E+06, 0.47067E+06, 0.47374E+06, 0.47683E+06 &
, 0.47992E+06, 0.48304E+06, 0.48617E+06, 0.48932E+06, 0.49248E+06, 0.49566E+06, 0.49886E+06, 0.50207E+06, 0.50530E+06, 0.50855E+06 &
, 0.51181E+06, 0.51509E+06, 0.51838E+06, 0.52169E+06, 0.52502E+06, 0.52837E+06, 0.53173E+06, 0.53511E+06, 0.53851E+06, 0.54192E+06 &
, 0.54535E+06, 0.54880E+06, 0.55227E+06, 0.55575E+06, 0.55925E+06, 0.56277E+06, 0.56630E+06, 0.56986E+06, 0.57343E+06, 0.57702E+06 &
, 0.58062E+06, 0.58425E+06, 0.58789E+06, 0.59155E+06, 0.59523E+06, 0.59893E+06, 0.60264E+06, 0.60638E+06, 0.61013E+06, 0.61390E+06 &
, 0.61768E+06, 0.62149E+06, 0.62532E+06, 0.62916E+06, 0.63302E+06, 0.63691E+06, 0.64081E+06, 0.64473E+06, 0.64867E+06, 0.65262E+06 &
, 0.65660E+06, 0.66060E+06, 0.66462E+06, 0.66865E+06, 0.67270E+06, 0.67678E+06, 0.68087E+06, 0.68499E+06, 0.68912E+06, 0.69327E+06 &
, 0.69744E+06, 0.70164E+06, 0.70585E+06, 0.71008E+06, 0.71434E+06, 0.71861E+06, 0.72290E+06, 0.72722E+06, 0.73155E+06, 0.73591E+06 &
, 0.74028E+06, 0.74468E+06, 0.74909E+06, 0.75353E+06, 0.75799E+06, 0.76247E+06, 0.76697E+06, 0.77149E+06, 0.77603E+06, 0.78059E+06 /
 !            3          10
 DATA(QofT(          32 ,J),J=1,501)/   760.640014648438       &
, 0.87685E+03, 0.99844E+03, 0.11252E+04, 0.12569E+04, 0.13934E+04, 0.15345E+04, 0.16801E+04, 0.18300E+04, 0.19842E+04, 0.21424E+04 &
, 0.23047E+04, 0.24709E+04, 0.26409E+04, 0.28146E+04, 0.29920E+04, 0.31730E+04, 0.33575E+04, 0.35455E+04, 0.37368E+04, 0.39315E+04 &
, 0.41295E+04, 0.43307E+04, 0.45350E+04, 0.47425E+04, 0.49531E+04, 0.51667E+04, 0.53833E+04, 0.56029E+04, 0.58254E+04, 0.60508E+04 &
, 0.62791E+04, 0.65101E+04, 0.67439E+04, 0.69805E+04, 0.72198E+04, 0.74618E+04, 0.77065E+04, 0.79538E+04, 0.82037E+04, 0.84562E+04 &
, 0.87113E+04, 0.89690E+04, 0.92291E+04, 0.94918E+04, 0.97569E+04, 0.10025E+05, 0.10295E+05, 0.10567E+05, 0.10842E+05, 0.11120E+05 &
, 0.11399E+05, 0.11682E+05, 0.11966E+05, 0.12253E+05, 0.12542E+05, 0.12834E+05, 0.13128E+05, 0.13424E+05, 0.13723E+05, 0.14024E+05 &
, 0.14327E+05, 0.14633E+05, 0.14941E+05, 0.15251E+05, 0.15564E+05, 0.15879E+05, 0.16196E+05, 0.16516E+05, 0.16838E+05, 0.17162E+05 &
, 0.17488E+05, 0.17817E+05, 0.18149E+05, 0.18482E+05, 0.18818E+05, 0.19157E+05, 0.19497E+05, 0.19840E+05, 0.20186E+05, 0.20534E+05 &
, 0.20884E+05, 0.21237E+05, 0.21592E+05, 0.21949E+05, 0.22309E+05, 0.22672E+05, 0.23037E+05, 0.23404E+05, 0.23774E+05, 0.24147E+05 &
, 0.24522E+05, 0.24900E+05, 0.25280E+05, 0.25663E+05, 0.26048E+05, 0.26436E+05, 0.26827E+05, 0.27220E+05, 0.27616E+05, 0.28015E+05 &
, 0.28416E+05, 0.28821E+05, 0.29228E+05, 0.29637E+05, 0.30050E+05, 0.30465E+05, 0.30884E+05, 0.31305E+05, 0.31729E+05, 0.32156E+05 &
, 0.32585E+05, 0.33018E+05, 0.33454E+05, 0.33892E+05, 0.34334E+05, 0.34779E+05, 0.35227E+05, 0.35678E+05, 0.36132E+05, 0.36590E+05 &
, 0.37050E+05, 0.37514E+05, 0.37981E+05, 0.38450E+05, 0.38924E+05, 0.39401E+05, 0.39881E+05, 0.40364E+05, 0.40851E+05, 0.41341E+05 &
, 0.41835E+05, 0.42331E+05, 0.42832E+05, 0.43336E+05, 0.43844E+05, 0.44355E+05, 0.44870E+05, 0.45388E+05, 0.45910E+05, 0.46435E+05 &
, 0.46965E+05, 0.47498E+05, 0.48035E+05, 0.48576E+05, 0.49120E+05, 0.49668E+05, 0.50221E+05, 0.50777E+05, 0.51337E+05, 0.51901E+05 &
, 0.52469E+05, 0.53041E+05, 0.53617E+05, 0.54197E+05, 0.54782E+05, 0.55370E+05, 0.55963E+05, 0.56560E+05, 0.57161E+05, 0.57766E+05 &
, 0.58375E+05, 0.58989E+05, 0.59608E+05, 0.60231E+05, 0.60857E+05, 0.61489E+05, 0.62125E+05, 0.62766E+05, 0.63411E+05, 0.64061E+05 &
, 0.64715E+05, 0.65374E+05, 0.66037E+05, 0.66706E+05, 0.67378E+05, 0.68057E+05, 0.68739E+05, 0.69426E+05, 0.70119E+05, 0.70816E+05 &
, 0.71518E+05, 0.72225E+05, 0.72937E+05, 0.73654E+05, 0.74376E+05, 0.75103E+05, 0.75836E+05, 0.76573E+05, 0.77316E+05, 0.78063E+05 &
, 0.78817E+05, 0.79575E+05, 0.80339E+05, 0.81108E+05, 0.81882E+05, 0.82662E+05, 0.83447E+05, 0.84238E+05, 0.85034E+05, 0.85837E+05 &
, 0.86644E+05, 0.87457E+05, 0.88276E+05, 0.89101E+05, 0.89930E+05, 0.90766E+05, 0.91608E+05, 0.92456E+05, 0.93310E+05, 0.94169E+05 &
, 0.95034E+05, 0.95905E+05, 0.96782E+05, 0.97666E+05, 0.98555E+05, 0.99451E+05, 0.10035E+06, 0.10126E+06, 0.10217E+06, 0.10309E+06 &
, 0.10402E+06, 0.10495E+06, 0.10589E+06, 0.10684E+06, 0.10779E+06, 0.10875E+06, 0.10971E+06, 0.11069E+06, 0.11166E+06, 0.11265E+06 &
, 0.11364E+06, 0.11464E+06, 0.11564E+06, 0.11666E+06, 0.11767E+06, 0.11870E+06, 0.11973E+06, 0.12077E+06, 0.12182E+06, 0.12287E+06 &
, 0.12393E+06, 0.12500E+06, 0.12607E+06, 0.12715E+06, 0.12824E+06, 0.12934E+06, 0.13044E+06, 0.13155E+06, 0.13267E+06, 0.13379E+06 &
, 0.13493E+06, 0.13607E+06, 0.13721E+06, 0.13837E+06, 0.13953E+06, 0.14070E+06, 0.14188E+06, 0.14307E+06, 0.14426E+06, 0.14546E+06 &
, 0.14667E+06, 0.14789E+06, 0.14911E+06, 0.15034E+06, 0.15158E+06, 0.15283E+06, 0.15409E+06, 0.15535E+06, 0.15663E+06, 0.15791E+06 &
, 0.15920E+06, 0.16049E+06, 0.16180E+06, 0.16311E+06, 0.16444E+06, 0.16577E+06, 0.16711E+06, 0.16846E+06, 0.16981E+06, 0.17118E+06 &
, 0.17255E+06, 0.17393E+06, 0.17532E+06, 0.17672E+06, 0.17813E+06, 0.17955E+06, 0.18098E+06, 0.18241E+06, 0.18386E+06, 0.18531E+06 &
, 0.18677E+06, 0.18825E+06, 0.18973E+06, 0.19122E+06, 0.19272E+06, 0.19422E+06, 0.19574E+06, 0.19727E+06, 0.19881E+06, 0.20035E+06 &
, 0.20191E+06, 0.20347E+06, 0.20505E+06, 0.20663E+06, 0.20823E+06, 0.20983E+06, 0.21144E+06, 0.21307E+06, 0.21470E+06, 0.21635E+06 &
, 0.21800E+06, 0.21966E+06, 0.22134E+06, 0.22302E+06, 0.22471E+06, 0.22642E+06, 0.22813E+06, 0.22985E+06, 0.23159E+06, 0.23333E+06 &
, 0.23509E+06, 0.23686E+06, 0.23863E+06, 0.24042E+06, 0.24222E+06, 0.24403E+06, 0.24584E+06, 0.24767E+06, 0.24951E+06, 0.25137E+06 &
, 0.25323E+06, 0.25510E+06, 0.25699E+06, 0.25888E+06, 0.26079E+06, 0.26271E+06, 0.26464E+06, 0.26657E+06, 0.26853E+06, 0.27049E+06 &
, 0.27246E+06, 0.27445E+06, 0.27645E+06, 0.27846E+06, 0.28047E+06, 0.28251E+06, 0.28455E+06, 0.28661E+06, 0.28867E+06, 0.29075E+06 &
, 0.29284E+06, 0.29495E+06, 0.29706E+06, 0.29919E+06, 0.30133E+06, 0.30348E+06, 0.30564E+06, 0.30782E+06, 0.31001E+06, 0.31221E+06 &
, 0.31442E+06, 0.31665E+06, 0.31888E+06, 0.32113E+06, 0.32340E+06, 0.32567E+06, 0.32796E+06, 0.33026E+06, 0.33258E+06, 0.33490E+06 &
, 0.33724E+06, 0.33960E+06, 0.34196E+06, 0.34434E+06, 0.34673E+06, 0.34914E+06, 0.35156E+06, 0.35399E+06, 0.35644E+06, 0.35890E+06 &
, 0.36137E+06, 0.36385E+06, 0.36635E+06, 0.36887E+06, 0.37139E+06, 0.37393E+06, 0.37649E+06, 0.37906E+06, 0.38164E+06, 0.38423E+06 &
, 0.38684E+06, 0.38947E+06, 0.39211E+06, 0.39476E+06, 0.39743E+06, 0.40011E+06, 0.40281E+06, 0.40552E+06, 0.40824E+06, 0.41098E+06 &
, 0.41373E+06, 0.41650E+06, 0.41929E+06, 0.42209E+06, 0.42490E+06, 0.42773E+06, 0.43057E+06, 0.43343E+06, 0.43630E+06, 0.43919E+06 &
, 0.44209E+06, 0.44501E+06, 0.44794E+06, 0.45090E+06, 0.45386E+06, 0.45684E+06, 0.45984E+06, 0.46285E+06, 0.46588E+06, 0.46892E+06 &
, 0.47198E+06, 0.47505E+06, 0.47814E+06, 0.48125E+06, 0.48437E+06, 0.48751E+06, 0.49066E+06, 0.49384E+06, 0.49702E+06, 0.50023E+06 &
, 0.50345E+06, 0.50668E+06, 0.50994E+06, 0.51320E+06, 0.51649E+06, 0.51979E+06, 0.52311E+06, 0.52645E+06, 0.52981E+06, 0.53318E+06 &
, 0.53657E+06, 0.53997E+06, 0.54339E+06, 0.54683E+06, 0.55029E+06, 0.55376E+06, 0.55725E+06, 0.56076E+06, 0.56429E+06, 0.56783E+06 &
, 0.57139E+06, 0.57497E+06, 0.57857E+06, 0.58219E+06, 0.58582E+06, 0.58947E+06, 0.59314E+06, 0.59683E+06, 0.60053E+06, 0.60425E+06 &
, 0.60800E+06, 0.61176E+06, 0.61554E+06, 0.61933E+06, 0.62315E+06, 0.62698E+06, 0.63084E+06, 0.63471E+06, 0.63860E+06, 0.64251E+06 &
, 0.64644E+06, 0.65039E+06, 0.65435E+06, 0.65834E+06, 0.66235E+06, 0.66637E+06, 0.67042E+06, 0.67448E+06, 0.67857E+06, 0.68267E+06 &
, 0.68679E+06, 0.69094E+06, 0.69510E+06, 0.69928E+06, 0.70348E+06, 0.70770E+06, 0.71195E+06, 0.71621E+06, 0.72049E+06, 0.72480E+06 &
, 0.72912E+06, 0.73346E+06, 0.73783E+06, 0.74222E+06, 0.74662E+06, 0.75105E+06, 0.75549E+06, 0.75996E+06, 0.76445E+06, 0.76896E+06 /
 !            3          11
 DATA(QofT(          33 ,J),J=1,501)/   4470.50000000000       &
, 0.51535E+04, 0.58681E+04, 0.66130E+04, 0.73871E+04, 0.81893E+04, 0.90186E+04, 0.98743E+04, 0.10755E+05, 0.11661E+05, 0.12591E+05 &
, 0.13545E+05, 0.14522E+05, 0.15521E+05, 0.16542E+05, 0.17584E+05, 0.18648E+05, 0.19733E+05, 0.20837E+05, 0.21962E+05, 0.23106E+05 &
, 0.24269E+05, 0.25452E+05, 0.26653E+05, 0.27872E+05, 0.29110E+05, 0.30365E+05, 0.31638E+05, 0.32929E+05, 0.34236E+05, 0.35561E+05 &
, 0.36903E+05, 0.38260E+05, 0.39635E+05, 0.41025E+05, 0.42432E+05, 0.43854E+05, 0.45292E+05, 0.46745E+05, 0.48214E+05, 0.49698E+05 &
, 0.51197E+05, 0.52711E+05, 0.54240E+05, 0.55784E+05, 0.57342E+05, 0.58915E+05, 0.60502E+05, 0.62104E+05, 0.63720E+05, 0.65350E+05 &
, 0.66994E+05, 0.68652E+05, 0.70325E+05, 0.72011E+05, 0.73711E+05, 0.75424E+05, 0.77152E+05, 0.78893E+05, 0.80648E+05, 0.82417E+05 &
, 0.84199E+05, 0.85995E+05, 0.87804E+05, 0.89628E+05, 0.91464E+05, 0.93314E+05, 0.95178E+05, 0.97056E+05, 0.98947E+05, 0.10085E+06 &
, 0.10277E+06, 0.10470E+06, 0.10665E+06, 0.10861E+06, 0.11058E+06, 0.11257E+06, 0.11457E+06, 0.11658E+06, 0.11861E+06, 0.12066E+06 &
, 0.12271E+06, 0.12478E+06, 0.12687E+06, 0.12897E+06, 0.13108E+06, 0.13321E+06, 0.13535E+06, 0.13751E+06, 0.13968E+06, 0.14187E+06 &
, 0.14407E+06, 0.14628E+06, 0.14852E+06, 0.15076E+06, 0.15302E+06, 0.15530E+06, 0.15759E+06, 0.15990E+06, 0.16222E+06, 0.16456E+06 &
, 0.16692E+06, 0.16929E+06, 0.17167E+06, 0.17408E+06, 0.17650E+06, 0.17893E+06, 0.18138E+06, 0.18385E+06, 0.18634E+06, 0.18884E+06 &
, 0.19136E+06, 0.19389E+06, 0.19645E+06, 0.19902E+06, 0.20161E+06, 0.20421E+06, 0.20684E+06, 0.20948E+06, 0.21214E+06, 0.21482E+06 &
, 0.21751E+06, 0.22023E+06, 0.22296E+06, 0.22571E+06, 0.22848E+06, 0.23128E+06, 0.23409E+06, 0.23691E+06, 0.23976E+06, 0.24263E+06 &
, 0.24552E+06, 0.24843E+06, 0.25136E+06, 0.25430E+06, 0.25727E+06, 0.26026E+06, 0.26327E+06, 0.26630E+06, 0.26936E+06, 0.27243E+06 &
, 0.27552E+06, 0.27864E+06, 0.28178E+06, 0.28494E+06, 0.28812E+06, 0.29133E+06, 0.29455E+06, 0.29780E+06, 0.30108E+06, 0.30437E+06 &
, 0.30769E+06, 0.31103E+06, 0.31440E+06, 0.31778E+06, 0.32120E+06, 0.32463E+06, 0.32809E+06, 0.33158E+06, 0.33509E+06, 0.33862E+06 &
, 0.34218E+06, 0.34576E+06, 0.34937E+06, 0.35300E+06, 0.35666E+06, 0.36035E+06, 0.36406E+06, 0.36779E+06, 0.37156E+06, 0.37535E+06 &
, 0.37916E+06, 0.38301E+06, 0.38688E+06, 0.39077E+06, 0.39470E+06, 0.39865E+06, 0.40262E+06, 0.40663E+06, 0.41067E+06, 0.41473E+06 &
, 0.41882E+06, 0.42294E+06, 0.42709E+06, 0.43127E+06, 0.43548E+06, 0.43971E+06, 0.44398E+06, 0.44828E+06, 0.45260E+06, 0.45696E+06 &
, 0.46134E+06, 0.46576E+06, 0.47021E+06, 0.47468E+06, 0.47919E+06, 0.48373E+06, 0.48830E+06, 0.49290E+06, 0.49754E+06, 0.50221E+06 &
, 0.50691E+06, 0.51164E+06, 0.51640E+06, 0.52120E+06, 0.52603E+06, 0.53089E+06, 0.53579E+06, 0.54072E+06, 0.54568E+06, 0.55068E+06 &
, 0.55571E+06, 0.56078E+06, 0.56588E+06, 0.57101E+06, 0.57619E+06, 0.58139E+06, 0.58663E+06, 0.59191E+06, 0.59723E+06, 0.60258E+06 &
, 0.60796E+06, 0.61338E+06, 0.61884E+06, 0.62433E+06, 0.62987E+06, 0.63543E+06, 0.64104E+06, 0.64669E+06, 0.65237E+06, 0.65809E+06 &
, 0.66385E+06, 0.66965E+06, 0.67548E+06, 0.68136E+06, 0.68727E+06, 0.69323E+06, 0.69922E+06, 0.70526E+06, 0.71133E+06, 0.71744E+06 &
, 0.72360E+06, 0.72979E+06, 0.73603E+06, 0.74231E+06, 0.74862E+06, 0.75499E+06, 0.76138E+06, 0.76783E+06, 0.77431E+06, 0.78084E+06 &
, 0.78741E+06, 0.79403E+06, 0.80069E+06, 0.80739E+06, 0.81413E+06, 0.82092E+06, 0.82775E+06, 0.83463E+06, 0.84155E+06, 0.84852E+06 &
, 0.85553E+06, 0.86258E+06, 0.86968E+06, 0.87683E+06, 0.88402E+06, 0.89127E+06, 0.89855E+06, 0.90588E+06, 0.91326E+06, 0.92069E+06 &
, 0.92816E+06, 0.93568E+06, 0.94325E+06, 0.95086E+06, 0.95853E+06, 0.96624E+06, 0.97401E+06, 0.98182E+06, 0.98968E+06, 0.99759E+06 &
, 0.10055E+07, 0.10136E+07, 0.10216E+07, 0.10297E+07, 0.10379E+07, 0.10461E+07, 0.10544E+07, 0.10627E+07, 0.10710E+07, 0.10795E+07 &
, 0.10879E+07, 0.10964E+07, 0.11050E+07, 0.11136E+07, 0.11223E+07, 0.11311E+07, 0.11399E+07, 0.11487E+07, 0.11576E+07, 0.11665E+07 &
, 0.11755E+07, 0.11846E+07, 0.11937E+07, 0.12029E+07, 0.12121E+07, 0.12214E+07, 0.12307E+07, 0.12401E+07, 0.12496E+07, 0.12591E+07 &
, 0.12687E+07, 0.12783E+07, 0.12880E+07, 0.12977E+07, 0.13075E+07, 0.13174E+07, 0.13273E+07, 0.13373E+07, 0.13473E+07, 0.13574E+07 &
, 0.13675E+07, 0.13778E+07, 0.13880E+07, 0.13984E+07, 0.14088E+07, 0.14192E+07, 0.14297E+07, 0.14403E+07, 0.14510E+07, 0.14617E+07 &
, 0.14724E+07, 0.14833E+07, 0.14942E+07, 0.15051E+07, 0.15161E+07, 0.15272E+07, 0.15384E+07, 0.15496E+07, 0.15609E+07, 0.15722E+07 &
, 0.15836E+07, 0.15951E+07, 0.16067E+07, 0.16183E+07, 0.16299E+07, 0.16417E+07, 0.16535E+07, 0.16654E+07, 0.16773E+07, 0.16893E+07 &
, 0.17014E+07, 0.17136E+07, 0.17258E+07, 0.17381E+07, 0.17504E+07, 0.17629E+07, 0.17753E+07, 0.17879E+07, 0.18006E+07, 0.18133E+07 &
, 0.18260E+07, 0.18389E+07, 0.18518E+07, 0.18648E+07, 0.18779E+07, 0.18910E+07, 0.19042E+07, 0.19175E+07, 0.19309E+07, 0.19443E+07 &
, 0.19579E+07, 0.19714E+07, 0.19851E+07, 0.19988E+07, 0.20126E+07, 0.20265E+07, 0.20405E+07, 0.20545E+07, 0.20687E+07, 0.20828E+07 &
, 0.20971E+07, 0.21115E+07, 0.21259E+07, 0.21404E+07, 0.21550E+07, 0.21697E+07, 0.21844E+07, 0.21992E+07, 0.22141E+07, 0.22291E+07 &
, 0.22442E+07, 0.22593E+07, 0.22745E+07, 0.22898E+07, 0.23052E+07, 0.23207E+07, 0.23363E+07, 0.23519E+07, 0.23676E+07, 0.23834E+07 &
, 0.23993E+07, 0.24153E+07, 0.24313E+07, 0.24475E+07, 0.24637E+07, 0.24800E+07, 0.24964E+07, 0.25129E+07, 0.25295E+07, 0.25461E+07 &
, 0.25629E+07, 0.25797E+07, 0.25966E+07, 0.26137E+07, 0.26307E+07, 0.26479E+07, 0.26652E+07, 0.26826E+07, 0.27000E+07, 0.27176E+07 &
, 0.27352E+07, 0.27529E+07, 0.27708E+07, 0.27887E+07, 0.28067E+07, 0.28248E+07, 0.28430E+07, 0.28612E+07, 0.28796E+07, 0.28981E+07 &
, 0.29166E+07, 0.29353E+07, 0.29541E+07, 0.29729E+07, 0.29918E+07, 0.30109E+07, 0.30300E+07, 0.30492E+07, 0.30686E+07, 0.30880E+07 &
, 0.31075E+07, 0.31271E+07, 0.31469E+07, 0.31667E+07, 0.31866E+07, 0.32066E+07, 0.32267E+07, 0.32469E+07, 0.32673E+07, 0.32877E+07 &
, 0.33082E+07, 0.33288E+07, 0.33496E+07, 0.33704E+07, 0.33913E+07, 0.34124E+07, 0.34335E+07, 0.34547E+07, 0.34761E+07, 0.34975E+07 &
, 0.35191E+07, 0.35407E+07, 0.35625E+07, 0.35844E+07, 0.36064E+07, 0.36284E+07, 0.36506E+07, 0.36729E+07, 0.36953E+07, 0.37179E+07 &
, 0.37405E+07, 0.37632E+07, 0.37861E+07, 0.38090E+07, 0.38321E+07, 0.38553E+07, 0.38786E+07, 0.39020E+07, 0.39255E+07, 0.39491E+07 &
, 0.39729E+07, 0.39967E+07, 0.40207E+07, 0.40447E+07, 0.40689E+07, 0.40932E+07, 0.41177E+07, 0.41422E+07, 0.41669E+07, 0.41916E+07 &
, 0.42165E+07, 0.42415E+07, 0.42667E+07, 0.42919E+07, 0.43173E+07, 0.43427E+07, 0.43683E+07, 0.43941E+07, 0.44199E+07, 0.44459E+07 /
 !            3          12
 DATA(QofT(          34 ,J),J=1,501)/   2261.80004882812       &
, 0.26074E+04, 0.29690E+04, 0.33459E+04, 0.37377E+04, 0.41436E+04, 0.45632E+04, 0.49962E+04, 0.54421E+04, 0.59005E+04, 0.63712E+04 &
, 0.68538E+04, 0.73480E+04, 0.78536E+04, 0.83702E+04, 0.88978E+04, 0.94361E+04, 0.99848E+04, 0.10544E+05, 0.11113E+05, 0.11692E+05 &
, 0.12281E+05, 0.12879E+05, 0.13487E+05, 0.14104E+05, 0.14730E+05, 0.15365E+05, 0.16010E+05, 0.16663E+05, 0.17324E+05, 0.17995E+05 &
, 0.18674E+05, 0.19361E+05, 0.20056E+05, 0.20760E+05, 0.21471E+05, 0.22191E+05, 0.22919E+05, 0.23654E+05, 0.24397E+05, 0.25148E+05 &
, 0.25907E+05, 0.26673E+05, 0.27447E+05, 0.28228E+05, 0.29017E+05, 0.29813E+05, 0.30616E+05, 0.31427E+05, 0.32244E+05, 0.33069E+05 &
, 0.33902E+05, 0.34741E+05, 0.35587E+05, 0.36440E+05, 0.37301E+05, 0.38168E+05, 0.39042E+05, 0.39924E+05, 0.40812E+05, 0.41707E+05 &
, 0.42609E+05, 0.43518E+05, 0.44434E+05, 0.45357E+05, 0.46286E+05, 0.47223E+05, 0.48167E+05, 0.49117E+05, 0.50074E+05, 0.51038E+05 &
, 0.52009E+05, 0.52987E+05, 0.53972E+05, 0.54964E+05, 0.55963E+05, 0.56969E+05, 0.57982E+05, 0.59002E+05, 0.60029E+05, 0.61063E+05 &
, 0.62105E+05, 0.63153E+05, 0.64209E+05, 0.65272E+05, 0.66342E+05, 0.67419E+05, 0.68504E+05, 0.69596E+05, 0.70696E+05, 0.71802E+05 &
, 0.72917E+05, 0.74039E+05, 0.75168E+05, 0.76306E+05, 0.77450E+05, 0.78603E+05, 0.79764E+05, 0.80932E+05, 0.82109E+05, 0.83292E+05 &
, 0.84484E+05, 0.85685E+05, 0.86893E+05, 0.88109E+05, 0.89334E+05, 0.90567E+05, 0.91809E+05, 0.93058E+05, 0.94316E+05, 0.95583E+05 &
, 0.96858E+05, 0.98143E+05, 0.99435E+05, 0.10074E+06, 0.10205E+06, 0.10337E+06, 0.10469E+06, 0.10603E+06, 0.10738E+06, 0.10873E+06 &
, 0.11010E+06, 0.11147E+06, 0.11286E+06, 0.11425E+06, 0.11565E+06, 0.11707E+06, 0.11849E+06, 0.11992E+06, 0.12136E+06, 0.12281E+06 &
, 0.12428E+06, 0.12575E+06, 0.12723E+06, 0.12872E+06, 0.13022E+06, 0.13174E+06, 0.13326E+06, 0.13479E+06, 0.13634E+06, 0.13789E+06 &
, 0.13946E+06, 0.14104E+06, 0.14263E+06, 0.14423E+06, 0.14584E+06, 0.14746E+06, 0.14909E+06, 0.15073E+06, 0.15239E+06, 0.15406E+06 &
, 0.15574E+06, 0.15743E+06, 0.15913E+06, 0.16084E+06, 0.16257E+06, 0.16431E+06, 0.16606E+06, 0.16782E+06, 0.16959E+06, 0.17138E+06 &
, 0.17318E+06, 0.17499E+06, 0.17682E+06, 0.17866E+06, 0.18051E+06, 0.18237E+06, 0.18425E+06, 0.18614E+06, 0.18804E+06, 0.18996E+06 &
, 0.19189E+06, 0.19383E+06, 0.19579E+06, 0.19776E+06, 0.19974E+06, 0.20174E+06, 0.20375E+06, 0.20578E+06, 0.20782E+06, 0.20987E+06 &
, 0.21194E+06, 0.21402E+06, 0.21612E+06, 0.21823E+06, 0.22036E+06, 0.22250E+06, 0.22466E+06, 0.22683E+06, 0.22901E+06, 0.23121E+06 &
, 0.23343E+06, 0.23566E+06, 0.23791E+06, 0.24017E+06, 0.24245E+06, 0.24475E+06, 0.24706E+06, 0.24938E+06, 0.25172E+06, 0.25408E+06 &
, 0.25646E+06, 0.25885E+06, 0.26125E+06, 0.26368E+06, 0.26612E+06, 0.26858E+06, 0.27105E+06, 0.27354E+06, 0.27605E+06, 0.27857E+06 &
, 0.28111E+06, 0.28367E+06, 0.28625E+06, 0.28885E+06, 0.29146E+06, 0.29409E+06, 0.29673E+06, 0.29940E+06, 0.30208E+06, 0.30478E+06 &
, 0.30750E+06, 0.31024E+06, 0.31300E+06, 0.31578E+06, 0.31857E+06, 0.32138E+06, 0.32421E+06, 0.32706E+06, 0.32993E+06, 0.33282E+06 &
, 0.33573E+06, 0.33866E+06, 0.34160E+06, 0.34457E+06, 0.34756E+06, 0.35056E+06, 0.35359E+06, 0.35663E+06, 0.35970E+06, 0.36279E+06 &
, 0.36589E+06, 0.36902E+06, 0.37217E+06, 0.37534E+06, 0.37853E+06, 0.38174E+06, 0.38497E+06, 0.38822E+06, 0.39150E+06, 0.39479E+06 &
, 0.39811E+06, 0.40145E+06, 0.40481E+06, 0.40819E+06, 0.41159E+06, 0.41502E+06, 0.41847E+06, 0.42194E+06, 0.42543E+06, 0.42894E+06 &
, 0.43248E+06, 0.43604E+06, 0.43963E+06, 0.44323E+06, 0.44686E+06, 0.45052E+06, 0.45419E+06, 0.45789E+06, 0.46161E+06, 0.46536E+06 &
, 0.46913E+06, 0.47293E+06, 0.47674E+06, 0.48059E+06, 0.48445E+06, 0.48835E+06, 0.49226E+06, 0.49620E+06, 0.50017E+06, 0.50416E+06 &
, 0.50817E+06, 0.51221E+06, 0.51628E+06, 0.52037E+06, 0.52448E+06, 0.52863E+06, 0.53279E+06, 0.53699E+06, 0.54121E+06, 0.54545E+06 &
, 0.54972E+06, 0.55402E+06, 0.55835E+06, 0.56270E+06, 0.56708E+06, 0.57148E+06, 0.57591E+06, 0.58037E+06, 0.58486E+06, 0.58937E+06 &
, 0.59391E+06, 0.59848E+06, 0.60308E+06, 0.60770E+06, 0.61236E+06, 0.61704E+06, 0.62175E+06, 0.62648E+06, 0.63125E+06, 0.63604E+06 &
, 0.64086E+06, 0.64572E+06, 0.65060E+06, 0.65551E+06, 0.66045E+06, 0.66542E+06, 0.67041E+06, 0.67544E+06, 0.68050E+06, 0.68559E+06 &
, 0.69071E+06, 0.69586E+06, 0.70104E+06, 0.70625E+06, 0.71149E+06, 0.71676E+06, 0.72206E+06, 0.72740E+06, 0.73276E+06, 0.73816E+06 &
, 0.74359E+06, 0.74904E+06, 0.75454E+06, 0.76006E+06, 0.76561E+06, 0.77120E+06, 0.77682E+06, 0.78247E+06, 0.78816E+06, 0.79387E+06 &
, 0.79963E+06, 0.80541E+06, 0.81123E+06, 0.81708E+06, 0.82296E+06, 0.82888E+06, 0.83483E+06, 0.84081E+06, 0.84683E+06, 0.85289E+06 &
, 0.85897E+06, 0.86509E+06, 0.87125E+06, 0.87745E+06, 0.88367E+06, 0.88993E+06, 0.89623E+06, 0.90256E+06, 0.90893E+06, 0.91533E+06 &
, 0.92177E+06, 0.92825E+06, 0.93476E+06, 0.94131E+06, 0.94789E+06, 0.95452E+06, 0.96117E+06, 0.96786E+06, 0.97460E+06, 0.98137E+06 &
, 0.98817E+06, 0.99502E+06, 0.10019E+07, 0.10088E+07, 0.10158E+07, 0.10228E+07, 0.10298E+07, 0.10369E+07, 0.10440E+07, 0.10511E+07 &
, 0.10583E+07, 0.10656E+07, 0.10728E+07, 0.10801E+07, 0.10875E+07, 0.10949E+07, 0.11023E+07, 0.11098E+07, 0.11173E+07, 0.11248E+07 &
, 0.11324E+07, 0.11400E+07, 0.11477E+07, 0.11554E+07, 0.11631E+07, 0.11709E+07, 0.11788E+07, 0.11866E+07, 0.11946E+07, 0.12025E+07 &
, 0.12105E+07, 0.12186E+07, 0.12266E+07, 0.12348E+07, 0.12429E+07, 0.12512E+07, 0.12594E+07, 0.12677E+07, 0.12761E+07, 0.12845E+07 &
, 0.12929E+07, 0.13014E+07, 0.13099E+07, 0.13184E+07, 0.13271E+07, 0.13357E+07, 0.13444E+07, 0.13531E+07, 0.13619E+07, 0.13708E+07 &
, 0.13796E+07, 0.13886E+07, 0.13975E+07, 0.14066E+07, 0.14156E+07, 0.14247E+07, 0.14339E+07, 0.14431E+07, 0.14523E+07, 0.14616E+07 &
, 0.14710E+07, 0.14804E+07, 0.14898E+07, 0.14993E+07, 0.15088E+07, 0.15184E+07, 0.15281E+07, 0.15377E+07, 0.15475E+07, 0.15572E+07 &
, 0.15671E+07, 0.15769E+07, 0.15869E+07, 0.15968E+07, 0.16069E+07, 0.16169E+07, 0.16271E+07, 0.16372E+07, 0.16475E+07, 0.16578E+07 &
, 0.16681E+07, 0.16785E+07, 0.16889E+07, 0.16994E+07, 0.17099E+07, 0.17205E+07, 0.17311E+07, 0.17418E+07, 0.17526E+07, 0.17634E+07 &
, 0.17742E+07, 0.17851E+07, 0.17961E+07, 0.18071E+07, 0.18181E+07, 0.18292E+07, 0.18404E+07, 0.18516E+07, 0.18629E+07, 0.18742E+07 &
, 0.18856E+07, 0.18971E+07, 0.19086E+07, 0.19201E+07, 0.19317E+07, 0.19434E+07, 0.19551E+07, 0.19669E+07, 0.19787E+07, 0.19906E+07 &
, 0.20025E+07, 0.20145E+07, 0.20266E+07, 0.20387E+07, 0.20509E+07, 0.20631E+07, 0.20754E+07, 0.20878E+07, 0.21002E+07, 0.21126E+07 &
, 0.21252E+07, 0.21377E+07, 0.21504E+07, 0.21631E+07, 0.21758E+07, 0.21887E+07, 0.22015E+07, 0.22145E+07, 0.22275E+07, 0.22405E+07 /
 !            3          13
 DATA(QofT(          35 ,J),J=1,501)/   69.9550018310547       &
, 0.80648E+02, 0.91836E+02, 0.10350E+03, 0.11562E+03, 0.12818E+03, 0.14116E+03, 0.15456E+03, 0.16836E+03, 0.18254E+03, 0.19711E+03 &
, 0.21204E+03, 0.22733E+03, 0.24298E+03, 0.25896E+03, 0.27529E+03, 0.29194E+03, 0.30892E+03, 0.32622E+03, 0.34383E+03, 0.36175E+03 &
, 0.37996E+03, 0.39848E+03, 0.41729E+03, 0.43638E+03, 0.45576E+03, 0.47542E+03, 0.49535E+03, 0.51556E+03, 0.53604E+03, 0.55678E+03 &
, 0.57779E+03, 0.59905E+03, 0.62057E+03, 0.64234E+03, 0.66437E+03, 0.68664E+03, 0.70916E+03, 0.73192E+03, 0.75492E+03, 0.77816E+03 &
, 0.80164E+03, 0.82536E+03, 0.84931E+03, 0.87349E+03, 0.89790E+03, 0.92254E+03, 0.94740E+03, 0.97249E+03, 0.99781E+03, 0.10234E+04 &
, 0.10491E+04, 0.10751E+04, 0.11013E+04, 0.11277E+04, 0.11544E+04, 0.11812E+04, 0.12083E+04, 0.12356E+04, 0.12631E+04, 0.12909E+04 &
, 0.13188E+04, 0.13470E+04, 0.13754E+04, 0.14040E+04, 0.14328E+04, 0.14619E+04, 0.14911E+04, 0.15206E+04, 0.15503E+04, 0.15802E+04 &
, 0.16103E+04, 0.16407E+04, 0.16712E+04, 0.17020E+04, 0.17330E+04, 0.17643E+04, 0.17957E+04, 0.18274E+04, 0.18593E+04, 0.18915E+04 &
, 0.19238E+04, 0.19564E+04, 0.19892E+04, 0.20223E+04, 0.20556E+04, 0.20891E+04, 0.21229E+04, 0.21569E+04, 0.21911E+04, 0.22256E+04 &
, 0.22603E+04, 0.22953E+04, 0.23305E+04, 0.23659E+04, 0.24016E+04, 0.24376E+04, 0.24738E+04, 0.25103E+04, 0.25470E+04, 0.25839E+04 &
, 0.26212E+04, 0.26587E+04, 0.26964E+04, 0.27344E+04, 0.27727E+04, 0.28113E+04, 0.28501E+04, 0.28892E+04, 0.29286E+04, 0.29683E+04 &
, 0.30082E+04, 0.30484E+04, 0.30890E+04, 0.31298E+04, 0.31708E+04, 0.32122E+04, 0.32539E+04, 0.32959E+04, 0.33382E+04, 0.33807E+04 &
, 0.34236E+04, 0.34668E+04, 0.35103E+04, 0.35541E+04, 0.35982E+04, 0.36426E+04, 0.36874E+04, 0.37325E+04, 0.37779E+04, 0.38236E+04 &
, 0.38697E+04, 0.39161E+04, 0.39628E+04, 0.40098E+04, 0.40573E+04, 0.41050E+04, 0.41531E+04, 0.42015E+04, 0.42503E+04, 0.42995E+04 &
, 0.43489E+04, 0.43988E+04, 0.44490E+04, 0.44996E+04, 0.45506E+04, 0.46019E+04, 0.46536E+04, 0.47057E+04, 0.47581E+04, 0.48109E+04 &
, 0.48642E+04, 0.49178E+04, 0.49718E+04, 0.50261E+04, 0.50809E+04, 0.51361E+04, 0.51917E+04, 0.52477E+04, 0.53041E+04, 0.53609E+04 &
, 0.54181E+04, 0.54757E+04, 0.55338E+04, 0.55922E+04, 0.56511E+04, 0.57104E+04, 0.57702E+04, 0.58304E+04, 0.58910E+04, 0.59521E+04 &
, 0.60136E+04, 0.60756E+04, 0.61380E+04, 0.62008E+04, 0.62641E+04, 0.63279E+04, 0.63921E+04, 0.64568E+04, 0.65220E+04, 0.65877E+04 &
, 0.66538E+04, 0.67203E+04, 0.67874E+04, 0.68550E+04, 0.69230E+04, 0.69915E+04, 0.70605E+04, 0.71301E+04, 0.72001E+04, 0.72706E+04 &
, 0.73416E+04, 0.74132E+04, 0.74852E+04, 0.75578E+04, 0.76308E+04, 0.77044E+04, 0.77785E+04, 0.78533E+04, 0.79284E+04, 0.80042E+04 &
, 0.80804E+04, 0.81572E+04, 0.82345E+04, 0.83124E+04, 0.83909E+04, 0.84699E+04, 0.85494E+04, 0.86296E+04, 0.87102E+04, 0.87915E+04 &
, 0.88733E+04, 0.89558E+04, 0.90387E+04, 0.91223E+04, 0.92064E+04, 0.92912E+04, 0.93765E+04, 0.94624E+04, 0.95490E+04, 0.96361E+04 &
, 0.97238E+04, 0.98121E+04, 0.99011E+04, 0.99907E+04, 0.10081E+05, 0.10172E+05, 0.10263E+05, 0.10355E+05, 0.10448E+05, 0.10541E+05 &
, 0.10635E+05, 0.10730E+05, 0.10825E+05, 0.10921E+05, 0.11018E+05, 0.11115E+05, 0.11213E+05, 0.11311E+05, 0.11411E+05, 0.11511E+05 &
, 0.11611E+05, 0.11713E+05, 0.11815E+05, 0.11917E+05, 0.12021E+05, 0.12125E+05, 0.12229E+05, 0.12335E+05, 0.12441E+05, 0.12548E+05 &
, 0.12656E+05, 0.12764E+05, 0.12873E+05, 0.12983E+05, 0.13093E+05, 0.13204E+05, 0.13316E+05, 0.13429E+05, 0.13543E+05, 0.13657E+05 &
, 0.13772E+05, 0.13887E+05, 0.14004E+05, 0.14121E+05, 0.14239E+05, 0.14358E+05, 0.14478E+05, 0.14598E+05, 0.14719E+05, 0.14841E+05 &
, 0.14964E+05, 0.15087E+05, 0.15212E+05, 0.15337E+05, 0.15463E+05, 0.15589E+05, 0.15717E+05, 0.15845E+05, 0.15975E+05, 0.16105E+05 &
, 0.16236E+05, 0.16367E+05, 0.16500E+05, 0.16633E+05, 0.16768E+05, 0.16903E+05, 0.17039E+05, 0.17176E+05, 0.17313E+05, 0.17452E+05 &
, 0.17591E+05, 0.17732E+05, 0.17873E+05, 0.18015E+05, 0.18158E+05, 0.18302E+05, 0.18447E+05, 0.18592E+05, 0.18739E+05, 0.18887E+05 &
, 0.19035E+05, 0.19185E+05, 0.19335E+05, 0.19486E+05, 0.19638E+05, 0.19791E+05, 0.19945E+05, 0.20101E+05, 0.20256E+05, 0.20414E+05 &
, 0.20571E+05, 0.20730E+05, 0.20890E+05, 0.21051E+05, 0.21213E+05, 0.21375E+05, 0.21539E+05, 0.21704E+05, 0.21870E+05, 0.22037E+05 &
, 0.22205E+05, 0.22373E+05, 0.22543E+05, 0.22714E+05, 0.22886E+05, 0.23059E+05, 0.23233E+05, 0.23408E+05, 0.23584E+05, 0.23761E+05 &
, 0.23939E+05, 0.24118E+05, 0.24298E+05, 0.24480E+05, 0.24662E+05, 0.24846E+05, 0.25030E+05, 0.25216E+05, 0.25403E+05, 0.25591E+05 &
, 0.25780E+05, 0.25970E+05, 0.26161E+05, 0.26353E+05, 0.26547E+05, 0.26741E+05, 0.26937E+05, 0.27134E+05, 0.27332E+05, 0.27531E+05 &
, 0.27731E+05, 0.27932E+05, 0.28135E+05, 0.28339E+05, 0.28544E+05, 0.28750E+05, 0.28957E+05, 0.29166E+05, 0.29375E+05, 0.29586E+05 &
, 0.29798E+05, 0.30011E+05, 0.30226E+05, 0.30442E+05, 0.30659E+05, 0.30877E+05, 0.31096E+05, 0.31317E+05, 0.31539E+05, 0.31762E+05 &
, 0.31986E+05, 0.32212E+05, 0.32439E+05, 0.32667E+05, 0.32897E+05, 0.33127E+05, 0.33359E+05, 0.33593E+05, 0.33827E+05, 0.34063E+05 &
, 0.34301E+05, 0.34539E+05, 0.34779E+05, 0.35020E+05, 0.35263E+05, 0.35507E+05, 0.35752E+05, 0.35998E+05, 0.36246E+05, 0.36496E+05 &
, 0.36746E+05, 0.36998E+05, 0.37251E+05, 0.37506E+05, 0.37762E+05, 0.38020E+05, 0.38279E+05, 0.38539E+05, 0.38801E+05, 0.39064E+05 &
, 0.39329E+05, 0.39595E+05, 0.39862E+05, 0.40131E+05, 0.40401E+05, 0.40673E+05, 0.40946E+05, 0.41221E+05, 0.41497E+05, 0.41775E+05 &
, 0.42054E+05, 0.42334E+05, 0.42616E+05, 0.42900E+05, 0.43185E+05, 0.43471E+05, 0.43759E+05, 0.44049E+05, 0.44340E+05, 0.44632E+05 &
, 0.44927E+05, 0.45222E+05, 0.45519E+05, 0.45818E+05, 0.46119E+05, 0.46421E+05, 0.46724E+05, 0.47029E+05, 0.47336E+05, 0.47644E+05 &
, 0.47954E+05, 0.48265E+05, 0.48578E+05, 0.48893E+05, 0.49209E+05, 0.49527E+05, 0.49847E+05, 0.50168E+05, 0.50490E+05, 0.50815E+05 &
, 0.51141E+05, 0.51469E+05, 0.51798E+05, 0.52129E+05, 0.52462E+05, 0.52796E+05, 0.53133E+05, 0.53470E+05, 0.53810E+05, 0.54151E+05 &
, 0.54494E+05, 0.54839E+05, 0.55185E+05, 0.55534E+05, 0.55883E+05, 0.56235E+05, 0.56589E+05, 0.56944E+05, 0.57301E+05, 0.57660E+05 &
, 0.58020E+05, 0.58382E+05, 0.58747E+05, 0.59112E+05, 0.59480E+05, 0.59850E+05, 0.60221E+05, 0.60594E+05, 0.60969E+05, 0.61346E+05 &
, 0.61725E+05, 0.62105E+05, 0.62488E+05, 0.62872E+05, 0.63258E+05, 0.63646E+05, 0.64036E+05, 0.64428E+05, 0.64822E+05, 0.65218E+05 &
, 0.65615E+05, 0.66015E+05, 0.66416E+05, 0.66819E+05, 0.67225E+05, 0.67632E+05, 0.68041E+05, 0.68452E+05, 0.68865E+05, 0.69281E+05 &
, 0.69698E+05, 0.70117E+05, 0.70538E+05, 0.70961E+05, 0.71386E+05, 0.71813E+05, 0.72242E+05, 0.72674E+05, 0.73107E+05, 0.73542E+05 /
 !            3          14
 DATA(QofT(          36 ,J),J=1,501)/   812.340026855469       &
, 0.93648E+03, 0.10664E+04, 0.12018E+04, 0.13425E+04, 0.14883E+04, 0.16391E+04, 0.17946E+04, 0.19548E+04, 0.21195E+04, 0.22886E+04 &
, 0.24620E+04, 0.26395E+04, 0.28211E+04, 0.30067E+04, 0.31963E+04, 0.33896E+04, 0.35868E+04, 0.37876E+04, 0.39920E+04, 0.42000E+04 &
, 0.44115E+04, 0.46265E+04, 0.48448E+04, 0.50665E+04, 0.52915E+04, 0.55198E+04, 0.57512E+04, 0.59858E+04, 0.62235E+04, 0.64644E+04 &
, 0.67082E+04, 0.69550E+04, 0.72049E+04, 0.74577E+04, 0.77133E+04, 0.79719E+04, 0.82334E+04, 0.84976E+04, 0.87647E+04, 0.90345E+04 &
, 0.93070E+04, 0.95824E+04, 0.98604E+04, 0.10141E+05, 0.10424E+05, 0.10710E+05, 0.10999E+05, 0.11290E+05, 0.11584E+05, 0.11881E+05 &
, 0.12180E+05, 0.12481E+05, 0.12786E+05, 0.13092E+05, 0.13402E+05, 0.13713E+05, 0.14028E+05, 0.14344E+05, 0.14664E+05, 0.14986E+05 &
, 0.15310E+05, 0.15637E+05, 0.15966E+05, 0.16298E+05, 0.16633E+05, 0.16970E+05, 0.17309E+05, 0.17651E+05, 0.17995E+05, 0.18342E+05 &
, 0.18692E+05, 0.19044E+05, 0.19398E+05, 0.19756E+05, 0.20115E+05, 0.20478E+05, 0.20842E+05, 0.21210E+05, 0.21580E+05, 0.21952E+05 &
, 0.22328E+05, 0.22706E+05, 0.23086E+05, 0.23469E+05, 0.23855E+05, 0.24244E+05, 0.24635E+05, 0.25029E+05, 0.25426E+05, 0.25826E+05 &
, 0.26228E+05, 0.26633E+05, 0.27041E+05, 0.27452E+05, 0.27866E+05, 0.28282E+05, 0.28701E+05, 0.29124E+05, 0.29549E+05, 0.29977E+05 &
, 0.30408E+05, 0.30842E+05, 0.31280E+05, 0.31720E+05, 0.32163E+05, 0.32610E+05, 0.33059E+05, 0.33512E+05, 0.33968E+05, 0.34427E+05 &
, 0.34889E+05, 0.35355E+05, 0.35824E+05, 0.36296E+05, 0.36771E+05, 0.37250E+05, 0.37732E+05, 0.38218E+05, 0.38707E+05, 0.39199E+05 &
, 0.39695E+05, 0.40195E+05, 0.40698E+05, 0.41204E+05, 0.41715E+05, 0.42228E+05, 0.42746E+05, 0.43267E+05, 0.43792E+05, 0.44320E+05 &
, 0.44853E+05, 0.45389E+05, 0.45929E+05, 0.46473E+05, 0.47021E+05, 0.47573E+05, 0.48128E+05, 0.48688E+05, 0.49252E+05, 0.49819E+05 &
, 0.50392E+05, 0.50967E+05, 0.51547E+05, 0.52132E+05, 0.52720E+05, 0.53313E+05, 0.53910E+05, 0.54511E+05, 0.55117E+05, 0.55727E+05 &
, 0.56341E+05, 0.56960E+05, 0.57583E+05, 0.58211E+05, 0.58843E+05, 0.59481E+05, 0.60122E+05, 0.60768E+05, 0.61419E+05, 0.62074E+05 &
, 0.62735E+05, 0.63400E+05, 0.64070E+05, 0.64744E+05, 0.65423E+05, 0.66108E+05, 0.66798E+05, 0.67492E+05, 0.68191E+05, 0.68895E+05 &
, 0.69605E+05, 0.70319E+05, 0.71039E+05, 0.71764E+05, 0.72494E+05, 0.73229E+05, 0.73970E+05, 0.74716E+05, 0.75467E+05, 0.76224E+05 &
, 0.76986E+05, 0.77754E+05, 0.78527E+05, 0.79305E+05, 0.80090E+05, 0.80880E+05, 0.81675E+05, 0.82476E+05, 0.83283E+05, 0.84095E+05 &
, 0.84914E+05, 0.85739E+05, 0.86569E+05, 0.87405E+05, 0.88246E+05, 0.89094E+05, 0.89948E+05, 0.90808E+05, 0.91675E+05, 0.92546E+05 &
, 0.93425E+05, 0.94309E+05, 0.95200E+05, 0.96097E+05, 0.97000E+05, 0.97911E+05, 0.98827E+05, 0.99749E+05, 0.10068E+06, 0.10161E+06 &
, 0.10256E+06, 0.10350E+06, 0.10446E+06, 0.10542E+06, 0.10639E+06, 0.10737E+06, 0.10835E+06, 0.10934E+06, 0.11033E+06, 0.11134E+06 &
, 0.11235E+06, 0.11336E+06, 0.11439E+06, 0.11542E+06, 0.11646E+06, 0.11750E+06, 0.11855E+06, 0.11961E+06, 0.12068E+06, 0.12175E+06 &
, 0.12283E+06, 0.12392E+06, 0.12502E+06, 0.12612E+06, 0.12724E+06, 0.12835E+06, 0.12948E+06, 0.13061E+06, 0.13176E+06, 0.13290E+06 &
, 0.13406E+06, 0.13523E+06, 0.13640E+06, 0.13758E+06, 0.13877E+06, 0.13997E+06, 0.14117E+06, 0.14238E+06, 0.14360E+06, 0.14483E+06 &
, 0.14607E+06, 0.14732E+06, 0.14857E+06, 0.14983E+06, 0.15110E+06, 0.15238E+06, 0.15367E+06, 0.15496E+06, 0.15627E+06, 0.15758E+06 &
, 0.15890E+06, 0.16023E+06, 0.16157E+06, 0.16292E+06, 0.16428E+06, 0.16564E+06, 0.16702E+06, 0.16840E+06, 0.16979E+06, 0.17119E+06 &
, 0.17260E+06, 0.17402E+06, 0.17545E+06, 0.17689E+06, 0.17834E+06, 0.17979E+06, 0.18126E+06, 0.18273E+06, 0.18422E+06, 0.18571E+06 &
, 0.18722E+06, 0.18873E+06, 0.19025E+06, 0.19179E+06, 0.19333E+06, 0.19488E+06, 0.19644E+06, 0.19801E+06, 0.19959E+06, 0.20119E+06 &
, 0.20279E+06, 0.20440E+06, 0.20602E+06, 0.20766E+06, 0.20930E+06, 0.21095E+06, 0.21261E+06, 0.21429E+06, 0.21597E+06, 0.21766E+06 &
, 0.21937E+06, 0.22108E+06, 0.22281E+06, 0.22455E+06, 0.22630E+06, 0.22805E+06, 0.22982E+06, 0.23160E+06, 0.23339E+06, 0.23520E+06 &
, 0.23701E+06, 0.23883E+06, 0.24067E+06, 0.24251E+06, 0.24437E+06, 0.24624E+06, 0.24812E+06, 0.25001E+06, 0.25191E+06, 0.25383E+06 &
, 0.25575E+06, 0.25769E+06, 0.25964E+06, 0.26160E+06, 0.26357E+06, 0.26556E+06, 0.26756E+06, 0.26956E+06, 0.27158E+06, 0.27361E+06 &
, 0.27566E+06, 0.27772E+06, 0.27978E+06, 0.28187E+06, 0.28396E+06, 0.28606E+06, 0.28818E+06, 0.29031E+06, 0.29245E+06, 0.29461E+06 &
, 0.29678E+06, 0.29896E+06, 0.30115E+06, 0.30336E+06, 0.30558E+06, 0.30781E+06, 0.31005E+06, 0.31231E+06, 0.31458E+06, 0.31687E+06 &
, 0.31916E+06, 0.32147E+06, 0.32380E+06, 0.32614E+06, 0.32849E+06, 0.33085E+06, 0.33323E+06, 0.33562E+06, 0.33802E+06, 0.34044E+06 &
, 0.34288E+06, 0.34532E+06, 0.34778E+06, 0.35025E+06, 0.35274E+06, 0.35525E+06, 0.35776E+06, 0.36029E+06, 0.36284E+06, 0.36539E+06 &
, 0.36797E+06, 0.37055E+06, 0.37316E+06, 0.37577E+06, 0.37840E+06, 0.38105E+06, 0.38371E+06, 0.38639E+06, 0.38908E+06, 0.39178E+06 &
, 0.39450E+06, 0.39724E+06, 0.39998E+06, 0.40275E+06, 0.40553E+06, 0.40833E+06, 0.41114E+06, 0.41396E+06, 0.41681E+06, 0.41966E+06 &
, 0.42253E+06, 0.42542E+06, 0.42833E+06, 0.43125E+06, 0.43418E+06, 0.43713E+06, 0.44010E+06, 0.44309E+06, 0.44609E+06, 0.44910E+06 &
, 0.45214E+06, 0.45518E+06, 0.45825E+06, 0.46133E+06, 0.46443E+06, 0.46754E+06, 0.47067E+06, 0.47382E+06, 0.47698E+06, 0.48016E+06 &
, 0.48336E+06, 0.48658E+06, 0.48981E+06, 0.49306E+06, 0.49632E+06, 0.49960E+06, 0.50291E+06, 0.50622E+06, 0.50956E+06, 0.51291E+06 &
, 0.51628E+06, 0.51967E+06, 0.52307E+06, 0.52650E+06, 0.52994E+06, 0.53340E+06, 0.53687E+06, 0.54037E+06, 0.54388E+06, 0.54741E+06 &
, 0.55096E+06, 0.55453E+06, 0.55811E+06, 0.56172E+06, 0.56534E+06, 0.56898E+06, 0.57264E+06, 0.57632E+06, 0.58002E+06, 0.58373E+06 &
, 0.58747E+06, 0.59122E+06, 0.59500E+06, 0.59879E+06, 0.60260E+06, 0.60643E+06, 0.61028E+06, 0.61415E+06, 0.61804E+06, 0.62195E+06 &
, 0.62587E+06, 0.62982E+06, 0.63379E+06, 0.63778E+06, 0.64179E+06, 0.64581E+06, 0.64986E+06, 0.65393E+06, 0.65802E+06, 0.66212E+06 &
, 0.66625E+06, 0.67040E+06, 0.67457E+06, 0.67876E+06, 0.68297E+06, 0.68720E+06, 0.69145E+06, 0.69573E+06, 0.70002E+06, 0.70434E+06 &
, 0.70867E+06, 0.71303E+06, 0.71741E+06, 0.72181E+06, 0.72623E+06, 0.73067E+06, 0.73514E+06, 0.73962E+06, 0.74413E+06, 0.74866E+06 &
, 0.75321E+06, 0.75778E+06, 0.76238E+06, 0.76700E+06, 0.77164E+06, 0.77630E+06, 0.78098E+06, 0.78569E+06, 0.79042E+06, 0.79517E+06 &
, 0.79994E+06, 0.80474E+06, 0.80956E+06, 0.81440E+06, 0.81927E+06, 0.82416E+06, 0.82907E+06, 0.83401E+06, 0.83896E+06, 0.84395E+06 /
 !            3          15
 DATA(QofT(          37 ,J),J=1,501)/   410.940002441406       &
, 0.47375E+03, 0.53948E+03, 0.60799E+03, 0.67919E+03, 0.75297E+03, 0.82925E+03, 0.90795E+03, 0.98899E+03, 0.10723E+04, 0.11579E+04 &
, 0.12456E+04, 0.13354E+04, 0.14273E+04, 0.15212E+04, 0.16171E+04, 0.17150E+04, 0.18147E+04, 0.19163E+04, 0.20198E+04, 0.21250E+04 &
, 0.22320E+04, 0.23408E+04, 0.24513E+04, 0.25634E+04, 0.26773E+04, 0.27927E+04, 0.29098E+04, 0.30285E+04, 0.31488E+04, 0.32707E+04 &
, 0.33941E+04, 0.35190E+04, 0.36454E+04, 0.37733E+04, 0.39026E+04, 0.40335E+04, 0.41658E+04, 0.42995E+04, 0.44346E+04, 0.45711E+04 &
, 0.47091E+04, 0.48483E+04, 0.49890E+04, 0.51311E+04, 0.52744E+04, 0.54192E+04, 0.55652E+04, 0.57126E+04, 0.58613E+04, 0.60113E+04 &
, 0.61627E+04, 0.63153E+04, 0.64692E+04, 0.66244E+04, 0.67809E+04, 0.69387E+04, 0.70977E+04, 0.72581E+04, 0.74197E+04, 0.75826E+04 &
, 0.77468E+04, 0.79121E+04, 0.80789E+04, 0.82468E+04, 0.84161E+04, 0.85866E+04, 0.87584E+04, 0.89315E+04, 0.91058E+04, 0.92814E+04 &
, 0.94583E+04, 0.96365E+04, 0.98160E+04, 0.99967E+04, 0.10179E+05, 0.10362E+05, 0.10547E+05, 0.10733E+05, 0.10920E+05, 0.11109E+05 &
, 0.11299E+05, 0.11490E+05, 0.11683E+05, 0.11877E+05, 0.12072E+05, 0.12269E+05, 0.12467E+05, 0.12666E+05, 0.12867E+05, 0.13069E+05 &
, 0.13273E+05, 0.13478E+05, 0.13685E+05, 0.13892E+05, 0.14102E+05, 0.14313E+05, 0.14525E+05, 0.14739E+05, 0.14954E+05, 0.15171E+05 &
, 0.15389E+05, 0.15609E+05, 0.15830E+05, 0.16053E+05, 0.16278E+05, 0.16504E+05, 0.16731E+05, 0.16960E+05, 0.17191E+05, 0.17423E+05 &
, 0.17657E+05, 0.17893E+05, 0.18130E+05, 0.18369E+05, 0.18610E+05, 0.18852E+05, 0.19096E+05, 0.19342E+05, 0.19590E+05, 0.19839E+05 &
, 0.20090E+05, 0.20343E+05, 0.20597E+05, 0.20854E+05, 0.21112E+05, 0.21372E+05, 0.21634E+05, 0.21898E+05, 0.22163E+05, 0.22431E+05 &
, 0.22700E+05, 0.22971E+05, 0.23245E+05, 0.23520E+05, 0.23797E+05, 0.24076E+05, 0.24358E+05, 0.24641E+05, 0.24926E+05, 0.25213E+05 &
, 0.25502E+05, 0.25794E+05, 0.26087E+05, 0.26383E+05, 0.26680E+05, 0.26980E+05, 0.27282E+05, 0.27586E+05, 0.27893E+05, 0.28201E+05 &
, 0.28512E+05, 0.28825E+05, 0.29140E+05, 0.29458E+05, 0.29778E+05, 0.30100E+05, 0.30424E+05, 0.30751E+05, 0.31080E+05, 0.31412E+05 &
, 0.31746E+05, 0.32082E+05, 0.32421E+05, 0.32762E+05, 0.33106E+05, 0.33452E+05, 0.33800E+05, 0.34151E+05, 0.34505E+05, 0.34861E+05 &
, 0.35220E+05, 0.35581E+05, 0.35945E+05, 0.36311E+05, 0.36680E+05, 0.37052E+05, 0.37426E+05, 0.37804E+05, 0.38183E+05, 0.38566E+05 &
, 0.38951E+05, 0.39339E+05, 0.39730E+05, 0.40123E+05, 0.40520E+05, 0.40919E+05, 0.41321E+05, 0.41726E+05, 0.42133E+05, 0.42544E+05 &
, 0.42958E+05, 0.43374E+05, 0.43793E+05, 0.44216E+05, 0.44642E+05, 0.45070E+05, 0.45502E+05, 0.45936E+05, 0.46373E+05, 0.46814E+05 &
, 0.47258E+05, 0.47705E+05, 0.48155E+05, 0.48608E+05, 0.49064E+05, 0.49524E+05, 0.49987E+05, 0.50453E+05, 0.50922E+05, 0.51395E+05 &
, 0.51871E+05, 0.52350E+05, 0.52833E+05, 0.53318E+05, 0.53807E+05, 0.54300E+05, 0.54796E+05, 0.55296E+05, 0.55799E+05, 0.56305E+05 &
, 0.56815E+05, 0.57328E+05, 0.57845E+05, 0.58366E+05, 0.58890E+05, 0.59418E+05, 0.59949E+05, 0.60484E+05, 0.61023E+05, 0.61565E+05 &
, 0.62111E+05, 0.62661E+05, 0.63214E+05, 0.63772E+05, 0.64332E+05, 0.64897E+05, 0.65466E+05, 0.66038E+05, 0.66615E+05, 0.67195E+05 &
, 0.67779E+05, 0.68367E+05, 0.68959E+05, 0.69555E+05, 0.70155E+05, 0.70759E+05, 0.71367E+05, 0.71979E+05, 0.72595E+05, 0.73216E+05 &
, 0.73840E+05, 0.74469E+05, 0.75102E+05, 0.75738E+05, 0.76379E+05, 0.77024E+05, 0.77674E+05, 0.78328E+05, 0.78986E+05, 0.79649E+05 &
, 0.80315E+05, 0.80987E+05, 0.81662E+05, 0.82342E+05, 0.83027E+05, 0.83715E+05, 0.84409E+05, 0.85107E+05, 0.85809E+05, 0.86516E+05 &
, 0.87228E+05, 0.87944E+05, 0.88664E+05, 0.89390E+05, 0.90120E+05, 0.90855E+05, 0.91594E+05, 0.92339E+05, 0.93087E+05, 0.93841E+05 &
, 0.94600E+05, 0.95363E+05, 0.96131E+05, 0.96904E+05, 0.97682E+05, 0.98465E+05, 0.99253E+05, 0.10005E+06, 0.10084E+06, 0.10165E+06 &
, 0.10245E+06, 0.10327E+06, 0.10408E+06, 0.10491E+06, 0.10574E+06, 0.10657E+06, 0.10741E+06, 0.10825E+06, 0.10910E+06, 0.10996E+06 &
, 0.11082E+06, 0.11168E+06, 0.11255E+06, 0.11343E+06, 0.11431E+06, 0.11519E+06, 0.11609E+06, 0.11698E+06, 0.11789E+06, 0.11879E+06 &
, 0.11971E+06, 0.12063E+06, 0.12155E+06, 0.12248E+06, 0.12342E+06, 0.12436E+06, 0.12531E+06, 0.12626E+06, 0.12722E+06, 0.12819E+06 &
, 0.12916E+06, 0.13014E+06, 0.13112E+06, 0.13211E+06, 0.13310E+06, 0.13410E+06, 0.13511E+06, 0.13612E+06, 0.13714E+06, 0.13816E+06 &
, 0.13919E+06, 0.14023E+06, 0.14127E+06, 0.14232E+06, 0.14338E+06, 0.14444E+06, 0.14550E+06, 0.14658E+06, 0.14766E+06, 0.14874E+06 &
, 0.14984E+06, 0.15094E+06, 0.15204E+06, 0.15315E+06, 0.15427E+06, 0.15540E+06, 0.15653E+06, 0.15767E+06, 0.15881E+06, 0.15996E+06 &
, 0.16112E+06, 0.16228E+06, 0.16345E+06, 0.16463E+06, 0.16582E+06, 0.16701E+06, 0.16820E+06, 0.16941E+06, 0.17062E+06, 0.17184E+06 &
, 0.17307E+06, 0.17430E+06, 0.17554E+06, 0.17678E+06, 0.17804E+06, 0.17930E+06, 0.18057E+06, 0.18184E+06, 0.18312E+06, 0.18441E+06 &
, 0.18571E+06, 0.18701E+06, 0.18832E+06, 0.18964E+06, 0.19097E+06, 0.19230E+06, 0.19364E+06, 0.19499E+06, 0.19634E+06, 0.19770E+06 &
, 0.19907E+06, 0.20045E+06, 0.20184E+06, 0.20323E+06, 0.20463E+06, 0.20604E+06, 0.20745E+06, 0.20888E+06, 0.21031E+06, 0.21175E+06 &
, 0.21319E+06, 0.21465E+06, 0.21611E+06, 0.21758E+06, 0.21906E+06, 0.22055E+06, 0.22204E+06, 0.22355E+06, 0.22506E+06, 0.22658E+06 &
, 0.22810E+06, 0.22964E+06, 0.23118E+06, 0.23273E+06, 0.23429E+06, 0.23586E+06, 0.23744E+06, 0.23902E+06, 0.24062E+06, 0.24222E+06 &
, 0.24383E+06, 0.24545E+06, 0.24707E+06, 0.24871E+06, 0.25035E+06, 0.25201E+06, 0.25367E+06, 0.25534E+06, 0.25702E+06, 0.25871E+06 &
, 0.26040E+06, 0.26211E+06, 0.26382E+06, 0.26555E+06, 0.26728E+06, 0.26902E+06, 0.27077E+06, 0.27253E+06, 0.27430E+06, 0.27608E+06 &
, 0.27786E+06, 0.27966E+06, 0.28147E+06, 0.28328E+06, 0.28510E+06, 0.28694E+06, 0.28878E+06, 0.29063E+06, 0.29249E+06, 0.29436E+06 &
, 0.29625E+06, 0.29813E+06, 0.30003E+06, 0.30194E+06, 0.30386E+06, 0.30579E+06, 0.30773E+06, 0.30968E+06, 0.31163E+06, 0.31360E+06 &
, 0.31558E+06, 0.31757E+06, 0.31956E+06, 0.32157E+06, 0.32359E+06, 0.32561E+06, 0.32765E+06, 0.32970E+06, 0.33176E+06, 0.33382E+06 &
, 0.33590E+06, 0.33799E+06, 0.34009E+06, 0.34220E+06, 0.34432E+06, 0.34645E+06, 0.34859E+06, 0.35074E+06, 0.35290E+06, 0.35507E+06 &
, 0.35725E+06, 0.35945E+06, 0.36165E+06, 0.36386E+06, 0.36609E+06, 0.36832E+06, 0.37057E+06, 0.37283E+06, 0.37510E+06, 0.37738E+06 &
, 0.37967E+06, 0.38197E+06, 0.38428E+06, 0.38660E+06, 0.38894E+06, 0.39129E+06, 0.39364E+06, 0.39601E+06, 0.39839E+06, 0.40078E+06 &
, 0.40318E+06, 0.40560E+06, 0.40802E+06, 0.41046E+06, 0.41291E+06, 0.41537E+06, 0.41784E+06, 0.42032E+06, 0.42282E+06, 0.42533E+06 /
 !            3          16
 DATA(QofT(          38 ,J),J=1,501)/   4775.79980468750       &
, 0.55056E+04, 0.62693E+04, 0.70654E+04, 0.78927E+04, 0.87500E+04, 0.96363E+04, 0.10551E+05, 0.11492E+05, 0.12461E+05, 0.13455E+05 &
, 0.14474E+05, 0.15518E+05, 0.16585E+05, 0.17677E+05, 0.18791E+05, 0.19928E+05, 0.21086E+05, 0.22267E+05, 0.23469E+05, 0.24692E+05 &
, 0.25935E+05, 0.27199E+05, 0.28483E+05, 0.29786E+05, 0.31109E+05, 0.32450E+05, 0.33811E+05, 0.35190E+05, 0.36588E+05, 0.38004E+05 &
, 0.39437E+05, 0.40888E+05, 0.42357E+05, 0.43843E+05, 0.45346E+05, 0.46867E+05, 0.48403E+05, 0.49957E+05, 0.51527E+05, 0.53113E+05 &
, 0.54716E+05, 0.56334E+05, 0.57968E+05, 0.59618E+05, 0.61284E+05, 0.62965E+05, 0.64662E+05, 0.66374E+05, 0.68102E+05, 0.69845E+05 &
, 0.71603E+05, 0.73376E+05, 0.75164E+05, 0.76967E+05, 0.78785E+05, 0.80617E+05, 0.82465E+05, 0.84327E+05, 0.86204E+05, 0.88096E+05 &
, 0.90002E+05, 0.91924E+05, 0.93860E+05, 0.95810E+05, 0.97775E+05, 0.99755E+05, 0.10175E+06, 0.10376E+06, 0.10578E+06, 0.10782E+06 &
, 0.10988E+06, 0.11194E+06, 0.11403E+06, 0.11613E+06, 0.11824E+06, 0.12037E+06, 0.12251E+06, 0.12467E+06, 0.12684E+06, 0.12903E+06 &
, 0.13123E+06, 0.13345E+06, 0.13569E+06, 0.13794E+06, 0.14021E+06, 0.14249E+06, 0.14479E+06, 0.14710E+06, 0.14943E+06, 0.15177E+06 &
, 0.15414E+06, 0.15651E+06, 0.15891E+06, 0.16132E+06, 0.16375E+06, 0.16619E+06, 0.16865E+06, 0.17113E+06, 0.17362E+06, 0.17614E+06 &
, 0.17867E+06, 0.18121E+06, 0.18378E+06, 0.18636E+06, 0.18896E+06, 0.19158E+06, 0.19421E+06, 0.19687E+06, 0.19954E+06, 0.20223E+06 &
, 0.20494E+06, 0.20767E+06, 0.21042E+06, 0.21319E+06, 0.21597E+06, 0.21878E+06, 0.22160E+06, 0.22445E+06, 0.22731E+06, 0.23020E+06 &
, 0.23310E+06, 0.23603E+06, 0.23897E+06, 0.24194E+06, 0.24493E+06, 0.24793E+06, 0.25096E+06, 0.25401E+06, 0.25709E+06, 0.26018E+06 &
, 0.26330E+06, 0.26643E+06, 0.26959E+06, 0.27278E+06, 0.27598E+06, 0.27921E+06, 0.28246E+06, 0.28573E+06, 0.28903E+06, 0.29235E+06 &
, 0.29569E+06, 0.29906E+06, 0.30245E+06, 0.30587E+06, 0.30930E+06, 0.31277E+06, 0.31626E+06, 0.31977E+06, 0.32331E+06, 0.32688E+06 &
, 0.33046E+06, 0.33408E+06, 0.33772E+06, 0.34139E+06, 0.34508E+06, 0.34880E+06, 0.35255E+06, 0.35632E+06, 0.36012E+06, 0.36394E+06 &
, 0.36780E+06, 0.37168E+06, 0.37559E+06, 0.37953E+06, 0.38349E+06, 0.38748E+06, 0.39150E+06, 0.39555E+06, 0.39964E+06, 0.40375E+06 &
, 0.40788E+06, 0.41205E+06, 0.41625E+06, 0.42047E+06, 0.42473E+06, 0.42902E+06, 0.43333E+06, 0.43768E+06, 0.44206E+06, 0.44647E+06 &
, 0.45091E+06, 0.45539E+06, 0.45989E+06, 0.46443E+06, 0.46900E+06, 0.47360E+06, 0.47823E+06, 0.48290E+06, 0.48760E+06, 0.49233E+06 &
, 0.49710E+06, 0.50190E+06, 0.50673E+06, 0.51160E+06, 0.51650E+06, 0.52144E+06, 0.52641E+06, 0.53141E+06, 0.53645E+06, 0.54153E+06 &
, 0.54664E+06, 0.55179E+06, 0.55697E+06, 0.56219E+06, 0.56744E+06, 0.57274E+06, 0.57807E+06, 0.58343E+06, 0.58884E+06, 0.59428E+06 &
, 0.59976E+06, 0.60527E+06, 0.61083E+06, 0.61643E+06, 0.62206E+06, 0.62772E+06, 0.63344E+06, 0.63919E+06, 0.64497E+06, 0.65080E+06 &
, 0.65667E+06, 0.66258E+06, 0.66853E+06, 0.67452E+06, 0.68055E+06, 0.68662E+06, 0.69273E+06, 0.69889E+06, 0.70509E+06, 0.71133E+06 &
, 0.71761E+06, 0.72393E+06, 0.73030E+06, 0.73671E+06, 0.74316E+06, 0.74966E+06, 0.75620E+06, 0.76278E+06, 0.76941E+06, 0.77608E+06 &
, 0.78280E+06, 0.78956E+06, 0.79637E+06, 0.80322E+06, 0.81012E+06, 0.81706E+06, 0.82406E+06, 0.83109E+06, 0.83818E+06, 0.84531E+06 &
, 0.85249E+06, 0.85971E+06, 0.86698E+06, 0.87430E+06, 0.88167E+06, 0.88909E+06, 0.89656E+06, 0.90407E+06, 0.91163E+06, 0.91925E+06 &
, 0.92691E+06, 0.93462E+06, 0.94239E+06, 0.95020E+06, 0.95807E+06, 0.96598E+06, 0.97394E+06, 0.98196E+06, 0.99003E+06, 0.99815E+06 &
, 0.10063E+07, 0.10146E+07, 0.10228E+07, 0.10312E+07, 0.10396E+07, 0.10480E+07, 0.10565E+07, 0.10650E+07, 0.10736E+07, 0.10823E+07 &
, 0.10910E+07, 0.10998E+07, 0.11086E+07, 0.11175E+07, 0.11264E+07, 0.11354E+07, 0.11444E+07, 0.11535E+07, 0.11627E+07, 0.11719E+07 &
, 0.11812E+07, 0.11905E+07, 0.11999E+07, 0.12094E+07, 0.12189E+07, 0.12285E+07, 0.12381E+07, 0.12478E+07, 0.12575E+07, 0.12673E+07 &
, 0.12772E+07, 0.12871E+07, 0.12971E+07, 0.13072E+07, 0.13173E+07, 0.13274E+07, 0.13377E+07, 0.13480E+07, 0.13583E+07, 0.13688E+07 &
, 0.13793E+07, 0.13898E+07, 0.14004E+07, 0.14111E+07, 0.14218E+07, 0.14327E+07, 0.14435E+07, 0.14545E+07, 0.14655E+07, 0.14765E+07 &
, 0.14877E+07, 0.14989E+07, 0.15102E+07, 0.15215E+07, 0.15329E+07, 0.15444E+07, 0.15559E+07, 0.15675E+07, 0.15792E+07, 0.15910E+07 &
, 0.16028E+07, 0.16147E+07, 0.16266E+07, 0.16386E+07, 0.16507E+07, 0.16629E+07, 0.16752E+07, 0.16875E+07, 0.16999E+07, 0.17123E+07 &
, 0.17248E+07, 0.17374E+07, 0.17501E+07, 0.17629E+07, 0.17757E+07, 0.17886E+07, 0.18015E+07, 0.18146E+07, 0.18277E+07, 0.18409E+07 &
, 0.18542E+07, 0.18675E+07, 0.18810E+07, 0.18945E+07, 0.19080E+07, 0.19217E+07, 0.19354E+07, 0.19492E+07, 0.19631E+07, 0.19771E+07 &
, 0.19911E+07, 0.20053E+07, 0.20195E+07, 0.20338E+07, 0.20481E+07, 0.20626E+07, 0.20771E+07, 0.20917E+07, 0.21064E+07, 0.21212E+07 &
, 0.21360E+07, 0.21510E+07, 0.21660E+07, 0.21811E+07, 0.21963E+07, 0.22115E+07, 0.22269E+07, 0.22423E+07, 0.22579E+07, 0.22735E+07 &
, 0.22892E+07, 0.23050E+07, 0.23208E+07, 0.23368E+07, 0.23528E+07, 0.23690E+07, 0.23852E+07, 0.24015E+07, 0.24179E+07, 0.24344E+07 &
, 0.24509E+07, 0.24676E+07, 0.24844E+07, 0.25012E+07, 0.25181E+07, 0.25352E+07, 0.25523E+07, 0.25695E+07, 0.25868E+07, 0.26042E+07 &
, 0.26217E+07, 0.26393E+07, 0.26569E+07, 0.26747E+07, 0.26926E+07, 0.27105E+07, 0.27286E+07, 0.27467E+07, 0.27650E+07, 0.27833E+07 &
, 0.28018E+07, 0.28203E+07, 0.28389E+07, 0.28577E+07, 0.28765E+07, 0.28954E+07, 0.29145E+07, 0.29336E+07, 0.29528E+07, 0.29721E+07 &
, 0.29916E+07, 0.30111E+07, 0.30307E+07, 0.30505E+07, 0.30703E+07, 0.30902E+07, 0.31103E+07, 0.31304E+07, 0.31507E+07, 0.31710E+07 &
, 0.31915E+07, 0.32120E+07, 0.32327E+07, 0.32535E+07, 0.32743E+07, 0.32953E+07, 0.33164E+07, 0.33376E+07, 0.33589E+07, 0.33803E+07 &
, 0.34018E+07, 0.34235E+07, 0.34452E+07, 0.34671E+07, 0.34890E+07, 0.35111E+07, 0.35333E+07, 0.35556E+07, 0.35780E+07, 0.36005E+07 &
, 0.36231E+07, 0.36458E+07, 0.36687E+07, 0.36917E+07, 0.37147E+07, 0.37379E+07, 0.37612E+07, 0.37847E+07, 0.38082E+07, 0.38319E+07 &
, 0.38557E+07, 0.38795E+07, 0.39036E+07, 0.39277E+07, 0.39519E+07, 0.39763E+07, 0.40008E+07, 0.40254E+07, 0.40501E+07, 0.40750E+07 &
, 0.40999E+07, 0.41250E+07, 0.41502E+07, 0.41756E+07, 0.42010E+07, 0.42266E+07, 0.42523E+07, 0.42781E+07, 0.43040E+07, 0.43301E+07 &
, 0.43563E+07, 0.43827E+07, 0.44091E+07, 0.44357E+07, 0.44624E+07, 0.44892E+07, 0.45162E+07, 0.45433E+07, 0.45705E+07, 0.45978E+07 &
, 0.46253E+07, 0.46529E+07, 0.46806E+07, 0.47085E+07, 0.47365E+07, 0.47647E+07, 0.47929E+07, 0.48213E+07, 0.48498E+07, 0.48785E+07 /
 !            3          17
 DATA(QofT(          39 ,J),J=1,501)/   2362.00000000000       &
, 0.27229E+04, 0.31006E+04, 0.34943E+04, 0.39033E+04, 0.43273E+04, 0.47656E+04, 0.52177E+04, 0.56834E+04, 0.61622E+04, 0.66537E+04 &
, 0.71577E+04, 0.76739E+04, 0.82019E+04, 0.87415E+04, 0.92925E+04, 0.98546E+04, 0.10428E+05, 0.11011E+05, 0.11606E+05, 0.12211E+05 &
, 0.12825E+05, 0.13450E+05, 0.14085E+05, 0.14730E+05, 0.15384E+05, 0.16047E+05, 0.16720E+05, 0.17402E+05, 0.18093E+05, 0.18793E+05 &
, 0.19502E+05, 0.20220E+05, 0.20946E+05, 0.21681E+05, 0.22424E+05, 0.23176E+05, 0.23936E+05, 0.24704E+05, 0.25480E+05, 0.26264E+05 &
, 0.27057E+05, 0.27857E+05, 0.28665E+05, 0.29481E+05, 0.30305E+05, 0.31136E+05, 0.31975E+05, 0.32822E+05, 0.33676E+05, 0.34538E+05 &
, 0.35407E+05, 0.36284E+05, 0.37168E+05, 0.38059E+05, 0.38958E+05, 0.39864E+05, 0.40778E+05, 0.41698E+05, 0.42626E+05, 0.43561E+05 &
, 0.44504E+05, 0.45454E+05, 0.46411E+05, 0.47375E+05, 0.48347E+05, 0.49326E+05, 0.50312E+05, 0.51305E+05, 0.52306E+05, 0.53314E+05 &
, 0.54329E+05, 0.55351E+05, 0.56382E+05, 0.57419E+05, 0.58463E+05, 0.59515E+05, 0.60575E+05, 0.61642E+05, 0.62716E+05, 0.63798E+05 &
, 0.64888E+05, 0.65985E+05, 0.67089E+05, 0.68201E+05, 0.69322E+05, 0.70450E+05, 0.71585E+05, 0.72729E+05, 0.73880E+05, 0.75039E+05 &
, 0.76207E+05, 0.77382E+05, 0.78565E+05, 0.79757E+05, 0.80957E+05, 0.82165E+05, 0.83381E+05, 0.84607E+05, 0.85840E+05, 0.87081E+05 &
, 0.88332E+05, 0.89590E+05, 0.90858E+05, 0.92135E+05, 0.93420E+05, 0.94714E+05, 0.96017E+05, 0.97329E+05, 0.98651E+05, 0.99980E+05 &
, 0.10132E+06, 0.10267E+06, 0.10403E+06, 0.10539E+06, 0.10677E+06, 0.10816E+06, 0.10956E+06, 0.11096E+06, 0.11238E+06, 0.11380E+06 &
, 0.11524E+06, 0.11669E+06, 0.11814E+06, 0.11961E+06, 0.12109E+06, 0.12257E+06, 0.12407E+06, 0.12558E+06, 0.12710E+06, 0.12863E+06 &
, 0.13017E+06, 0.13172E+06, 0.13328E+06, 0.13486E+06, 0.13644E+06, 0.13804E+06, 0.13964E+06, 0.14126E+06, 0.14289E+06, 0.14454E+06 &
, 0.14619E+06, 0.14785E+06, 0.14953E+06, 0.15122E+06, 0.15292E+06, 0.15464E+06, 0.15636E+06, 0.15810E+06, 0.15985E+06, 0.16161E+06 &
, 0.16339E+06, 0.16518E+06, 0.16698E+06, 0.16879E+06, 0.17062E+06, 0.17246E+06, 0.17431E+06, 0.17618E+06, 0.17806E+06, 0.17995E+06 &
, 0.18186E+06, 0.18378E+06, 0.18571E+06, 0.18766E+06, 0.18962E+06, 0.19160E+06, 0.19359E+06, 0.19559E+06, 0.19761E+06, 0.19965E+06 &
, 0.20169E+06, 0.20376E+06, 0.20583E+06, 0.20793E+06, 0.21003E+06, 0.21215E+06, 0.21429E+06, 0.21644E+06, 0.21861E+06, 0.22079E+06 &
, 0.22299E+06, 0.22521E+06, 0.22744E+06, 0.22968E+06, 0.23194E+06, 0.23422E+06, 0.23652E+06, 0.23883E+06, 0.24115E+06, 0.24350E+06 &
, 0.24586E+06, 0.24823E+06, 0.25063E+06, 0.25304E+06, 0.25546E+06, 0.25791E+06, 0.26037E+06, 0.26285E+06, 0.26534E+06, 0.26786E+06 &
, 0.27039E+06, 0.27294E+06, 0.27550E+06, 0.27809E+06, 0.28069E+06, 0.28331E+06, 0.28595E+06, 0.28861E+06, 0.29129E+06, 0.29399E+06 &
, 0.29670E+06, 0.29943E+06, 0.30218E+06, 0.30496E+06, 0.30775E+06, 0.31056E+06, 0.31338E+06, 0.31623E+06, 0.31910E+06, 0.32199E+06 &
, 0.32490E+06, 0.32782E+06, 0.33077E+06, 0.33374E+06, 0.33673E+06, 0.33974E+06, 0.34276E+06, 0.34581E+06, 0.34888E+06, 0.35198E+06 &
, 0.35509E+06, 0.35822E+06, 0.36138E+06, 0.36455E+06, 0.36775E+06, 0.37097E+06, 0.37422E+06, 0.37748E+06, 0.38076E+06, 0.38407E+06 &
, 0.38740E+06, 0.39075E+06, 0.39413E+06, 0.39752E+06, 0.40094E+06, 0.40439E+06, 0.40785E+06, 0.41134E+06, 0.41485E+06, 0.41839E+06 &
, 0.42194E+06, 0.42552E+06, 0.42913E+06, 0.43276E+06, 0.43641E+06, 0.44009E+06, 0.44379E+06, 0.44752E+06, 0.45127E+06, 0.45504E+06 &
, 0.45884E+06, 0.46267E+06, 0.46652E+06, 0.47039E+06, 0.47429E+06, 0.47822E+06, 0.48217E+06, 0.48614E+06, 0.49015E+06, 0.49417E+06 &
, 0.49823E+06, 0.50231E+06, 0.50642E+06, 0.51055E+06, 0.51471E+06, 0.51889E+06, 0.52310E+06, 0.52734E+06, 0.53161E+06, 0.53591E+06 &
, 0.54022E+06, 0.54457E+06, 0.54895E+06, 0.55335E+06, 0.55778E+06, 0.56225E+06, 0.56673E+06, 0.57125E+06, 0.57579E+06, 0.58037E+06 &
, 0.58497E+06, 0.58960E+06, 0.59426E+06, 0.59895E+06, 0.60366E+06, 0.60841E+06, 0.61319E+06, 0.61799E+06, 0.62283E+06, 0.62769E+06 &
, 0.63259E+06, 0.63752E+06, 0.64247E+06, 0.64746E+06, 0.65248E+06, 0.65753E+06, 0.66261E+06, 0.66772E+06, 0.67286E+06, 0.67803E+06 &
, 0.68324E+06, 0.68848E+06, 0.69374E+06, 0.69904E+06, 0.70438E+06, 0.70974E+06, 0.71514E+06, 0.72057E+06, 0.72603E+06, 0.73152E+06 &
, 0.73705E+06, 0.74261E+06, 0.74821E+06, 0.75383E+06, 0.75950E+06, 0.76519E+06, 0.77092E+06, 0.77668E+06, 0.78248E+06, 0.78831E+06 &
, 0.79418E+06, 0.80008E+06, 0.80602E+06, 0.81199E+06, 0.81799E+06, 0.82403E+06, 0.83011E+06, 0.83622E+06, 0.84237E+06, 0.84856E+06 &
, 0.85478E+06, 0.86103E+06, 0.86732E+06, 0.87365E+06, 0.88002E+06, 0.88642E+06, 0.89286E+06, 0.89934E+06, 0.90586E+06, 0.91241E+06 &
, 0.91899E+06, 0.92562E+06, 0.93229E+06, 0.93899E+06, 0.94573E+06, 0.95252E+06, 0.95933E+06, 0.96619E+06, 0.97309E+06, 0.98002E+06 &
, 0.98700E+06, 0.99401E+06, 0.10011E+07, 0.10082E+07, 0.10153E+07, 0.10225E+07, 0.10297E+07, 0.10369E+07, 0.10442E+07, 0.10516E+07 &
, 0.10590E+07, 0.10664E+07, 0.10738E+07, 0.10813E+07, 0.10889E+07, 0.10965E+07, 0.11041E+07, 0.11118E+07, 0.11195E+07, 0.11272E+07 &
, 0.11350E+07, 0.11429E+07, 0.11507E+07, 0.11587E+07, 0.11666E+07, 0.11746E+07, 0.11827E+07, 0.11908E+07, 0.11990E+07, 0.12071E+07 &
, 0.12154E+07, 0.12237E+07, 0.12320E+07, 0.12403E+07, 0.12488E+07, 0.12572E+07, 0.12657E+07, 0.12743E+07, 0.12829E+07, 0.12915E+07 &
, 0.13002E+07, 0.13089E+07, 0.13177E+07, 0.13265E+07, 0.13354E+07, 0.13443E+07, 0.13533E+07, 0.13623E+07, 0.13714E+07, 0.13805E+07 &
, 0.13897E+07, 0.13989E+07, 0.14081E+07, 0.14174E+07, 0.14268E+07, 0.14362E+07, 0.14457E+07, 0.14552E+07, 0.14647E+07, 0.14743E+07 &
, 0.14840E+07, 0.14937E+07, 0.15034E+07, 0.15132E+07, 0.15231E+07, 0.15330E+07, 0.15430E+07, 0.15530E+07, 0.15630E+07, 0.15731E+07 &
, 0.15833E+07, 0.15935E+07, 0.16038E+07, 0.16141E+07, 0.16245E+07, 0.16349E+07, 0.16454E+07, 0.16559E+07, 0.16665E+07, 0.16772E+07 &
, 0.16879E+07, 0.16986E+07, 0.17094E+07, 0.17203E+07, 0.17312E+07, 0.17422E+07, 0.17532E+07, 0.17643E+07, 0.17754E+07, 0.17866E+07 &
, 0.17978E+07, 0.18091E+07, 0.18205E+07, 0.18319E+07, 0.18434E+07, 0.18549E+07, 0.18665E+07, 0.18782E+07, 0.18898E+07, 0.19016E+07 &
, 0.19134E+07, 0.19253E+07, 0.19372E+07, 0.19492E+07, 0.19613E+07, 0.19734E+07, 0.19856E+07, 0.19978E+07, 0.20101E+07, 0.20225E+07 &
, 0.20349E+07, 0.20473E+07, 0.20599E+07, 0.20725E+07, 0.20851E+07, 0.20978E+07, 0.21106E+07, 0.21234E+07, 0.21363E+07, 0.21493E+07 &
, 0.21623E+07, 0.21754E+07, 0.21886E+07, 0.22018E+07, 0.22151E+07, 0.22284E+07, 0.22418E+07, 0.22553E+07, 0.22688E+07, 0.22824E+07 &
, 0.22961E+07, 0.23098E+07, 0.23236E+07, 0.23374E+07, 0.23514E+07, 0.23654E+07, 0.23794E+07, 0.23935E+07, 0.24077E+07, 0.24220E+07 /
 !            3          18
 DATA(QofT(          40 ,J),J=1,501)/   13880.0000000000       &
, 0.16000E+05, 0.18220E+05, 0.20533E+05, 0.22937E+05, 0.25428E+05, 0.28003E+05, 0.30660E+05, 0.33396E+05, 0.36210E+05, 0.39098E+05 &
, 0.42059E+05, 0.45092E+05, 0.48195E+05, 0.51366E+05, 0.54603E+05, 0.57906E+05, 0.61274E+05, 0.64704E+05, 0.68197E+05, 0.71750E+05 &
, 0.75363E+05, 0.79035E+05, 0.82765E+05, 0.86552E+05, 0.90395E+05, 0.94294E+05, 0.98247E+05, 0.10225E+06, 0.10632E+06, 0.11043E+06 &
, 0.11459E+06, 0.11881E+06, 0.12308E+06, 0.12740E+06, 0.13176E+06, 0.13618E+06, 0.14065E+06, 0.14516E+06, 0.14972E+06, 0.15433E+06 &
, 0.15899E+06, 0.16369E+06, 0.16844E+06, 0.17323E+06, 0.17807E+06, 0.18296E+06, 0.18789E+06, 0.19286E+06, 0.19788E+06, 0.20294E+06 &
, 0.20805E+06, 0.21320E+06, 0.21839E+06, 0.22363E+06, 0.22891E+06, 0.23424E+06, 0.23960E+06, 0.24501E+06, 0.25046E+06, 0.25596E+06 &
, 0.26149E+06, 0.26707E+06, 0.27270E+06, 0.27836E+06, 0.28407E+06, 0.28982E+06, 0.29561E+06, 0.30144E+06, 0.30732E+06, 0.31324E+06 &
, 0.31920E+06, 0.32521E+06, 0.33126E+06, 0.33735E+06, 0.34348E+06, 0.34966E+06, 0.35588E+06, 0.36214E+06, 0.36845E+06, 0.37480E+06 &
, 0.38120E+06, 0.38764E+06, 0.39412E+06, 0.40065E+06, 0.40723E+06, 0.41385E+06, 0.42051E+06, 0.42722E+06, 0.43398E+06, 0.44078E+06 &
, 0.44763E+06, 0.45453E+06, 0.46147E+06, 0.46846E+06, 0.47550E+06, 0.48259E+06, 0.48972E+06, 0.49690E+06, 0.50414E+06, 0.51142E+06 &
, 0.51875E+06, 0.52613E+06, 0.53356E+06, 0.54105E+06, 0.54858E+06, 0.55617E+06, 0.56380E+06, 0.57149E+06, 0.57923E+06, 0.58703E+06 &
, 0.59488E+06, 0.60278E+06, 0.61074E+06, 0.61875E+06, 0.62682E+06, 0.63494E+06, 0.64312E+06, 0.65136E+06, 0.65965E+06, 0.66800E+06 &
, 0.67641E+06, 0.68488E+06, 0.69340E+06, 0.70198E+06, 0.71063E+06, 0.71933E+06, 0.72809E+06, 0.73692E+06, 0.74580E+06, 0.75476E+06 &
, 0.76377E+06, 0.77284E+06, 0.78198E+06, 0.79118E+06, 0.80045E+06, 0.80978E+06, 0.81918E+06, 0.82865E+06, 0.83817E+06, 0.84777E+06 &
, 0.85743E+06, 0.86717E+06, 0.87697E+06, 0.88684E+06, 0.89677E+06, 0.90679E+06, 0.91687E+06, 0.92702E+06, 0.93725E+06, 0.94754E+06 &
, 0.95791E+06, 0.96835E+06, 0.97887E+06, 0.98945E+06, 0.10001E+07, 0.10109E+07, 0.10217E+07, 0.10326E+07, 0.10435E+07, 0.10546E+07 &
, 0.10657E+07, 0.10769E+07, 0.10882E+07, 0.10996E+07, 0.11110E+07, 0.11225E+07, 0.11341E+07, 0.11458E+07, 0.11576E+07, 0.11694E+07 &
, 0.11814E+07, 0.11934E+07, 0.12055E+07, 0.12177E+07, 0.12300E+07, 0.12424E+07, 0.12548E+07, 0.12673E+07, 0.12800E+07, 0.12927E+07 &
, 0.13055E+07, 0.13184E+07, 0.13314E+07, 0.13445E+07, 0.13576E+07, 0.13709E+07, 0.13843E+07, 0.13977E+07, 0.14113E+07, 0.14249E+07 &
, 0.14386E+07, 0.14525E+07, 0.14664E+07, 0.14804E+07, 0.14945E+07, 0.15088E+07, 0.15231E+07, 0.15375E+07, 0.15520E+07, 0.15667E+07 &
, 0.15814E+07, 0.15962E+07, 0.16111E+07, 0.16262E+07, 0.16413E+07, 0.16566E+07, 0.16719E+07, 0.16874E+07, 0.17029E+07, 0.17186E+07 &
, 0.17344E+07, 0.17503E+07, 0.17663E+07, 0.17824E+07, 0.17986E+07, 0.18149E+07, 0.18313E+07, 0.18479E+07, 0.18645E+07, 0.18813E+07 &
, 0.18982E+07, 0.19152E+07, 0.19323E+07, 0.19496E+07, 0.19669E+07, 0.19844E+07, 0.20020E+07, 0.20197E+07, 0.20375E+07, 0.20555E+07 &
, 0.20736E+07, 0.20917E+07, 0.21101E+07, 0.21285E+07, 0.21471E+07, 0.21657E+07, 0.21846E+07, 0.22035E+07, 0.22226E+07, 0.22418E+07 &
, 0.22611E+07, 0.22805E+07, 0.23001E+07, 0.23198E+07, 0.23396E+07, 0.23596E+07, 0.23797E+07, 0.23999E+07, 0.24203E+07, 0.24408E+07 &
, 0.24614E+07, 0.24822E+07, 0.25031E+07, 0.25242E+07, 0.25454E+07, 0.25667E+07, 0.25881E+07, 0.26097E+07, 0.26315E+07, 0.26534E+07 &
, 0.26754E+07, 0.26975E+07, 0.27199E+07, 0.27423E+07, 0.27649E+07, 0.27877E+07, 0.28106E+07, 0.28336E+07, 0.28568E+07, 0.28801E+07 &
, 0.29036E+07, 0.29272E+07, 0.29510E+07, 0.29749E+07, 0.29990E+07, 0.30233E+07, 0.30477E+07, 0.30722E+07, 0.30969E+07, 0.31218E+07 &
, 0.31468E+07, 0.31720E+07, 0.31973E+07, 0.32228E+07, 0.32485E+07, 0.32743E+07, 0.33003E+07, 0.33264E+07, 0.33527E+07, 0.33792E+07 &
, 0.34058E+07, 0.34326E+07, 0.34596E+07, 0.34867E+07, 0.35140E+07, 0.35415E+07, 0.35691E+07, 0.35969E+07, 0.36249E+07, 0.36531E+07 &
, 0.36814E+07, 0.37099E+07, 0.37386E+07, 0.37674E+07, 0.37965E+07, 0.38256E+07, 0.38550E+07, 0.38846E+07, 0.39143E+07, 0.39442E+07 &
, 0.39743E+07, 0.40046E+07, 0.40351E+07, 0.40657E+07, 0.40966E+07, 0.41276E+07, 0.41588E+07, 0.41902E+07, 0.42218E+07, 0.42535E+07 &
, 0.42855E+07, 0.43176E+07, 0.43500E+07, 0.43825E+07, 0.44152E+07, 0.44481E+07, 0.44812E+07, 0.45145E+07, 0.45481E+07, 0.45818E+07 &
, 0.46157E+07, 0.46497E+07, 0.46840E+07, 0.47185E+07, 0.47532E+07, 0.47881E+07, 0.48233E+07, 0.48586E+07, 0.48941E+07, 0.49298E+07 &
, 0.49657E+07, 0.50019E+07, 0.50382E+07, 0.50748E+07, 0.51115E+07, 0.51485E+07, 0.51857E+07, 0.52231E+07, 0.52607E+07, 0.52986E+07 &
, 0.53366E+07, 0.53749E+07, 0.54133E+07, 0.54521E+07, 0.54910E+07, 0.55301E+07, 0.55695E+07, 0.56091E+07, 0.56489E+07, 0.56889E+07 &
, 0.57292E+07, 0.57697E+07, 0.58104E+07, 0.58513E+07, 0.58925E+07, 0.59339E+07, 0.59756E+07, 0.60174E+07, 0.60596E+07, 0.61019E+07 &
, 0.61444E+07, 0.61873E+07, 0.62303E+07, 0.62736E+07, 0.63171E+07, 0.63609E+07, 0.64049E+07, 0.64491E+07, 0.64936E+07, 0.65383E+07 &
, 0.65833E+07, 0.66285E+07, 0.66740E+07, 0.67197E+07, 0.67657E+07, 0.68119E+07, 0.68584E+07, 0.69051E+07, 0.69521E+07, 0.69993E+07 &
, 0.70468E+07, 0.70945E+07, 0.71425E+07, 0.71908E+07, 0.72393E+07, 0.72881E+07, 0.73371E+07, 0.73864E+07, 0.74360E+07, 0.74858E+07 &
, 0.75359E+07, 0.75862E+07, 0.76369E+07, 0.76878E+07, 0.77389E+07, 0.77904E+07, 0.78421E+07, 0.78941E+07, 0.79463E+07, 0.79988E+07 &
, 0.80517E+07, 0.81047E+07, 0.81581E+07, 0.82117E+07, 0.82656E+07, 0.83198E+07, 0.83743E+07, 0.84291E+07, 0.84841E+07, 0.85395E+07 &
, 0.85951E+07, 0.86510E+07, 0.87072E+07, 0.87637E+07, 0.88205E+07, 0.88776E+07, 0.89349E+07, 0.89926E+07, 0.90505E+07, 0.91088E+07 &
, 0.91674E+07, 0.92262E+07, 0.92854E+07, 0.93448E+07, 0.94046E+07, 0.94647E+07, 0.95250E+07, 0.95857E+07, 0.96467E+07, 0.97080E+07 &
, 0.97696E+07, 0.98315E+07, 0.98937E+07, 0.99562E+07, 0.10019E+08, 0.10082E+08, 0.10146E+08, 0.10209E+08, 0.10274E+08, 0.10338E+08 &
, 0.10403E+08, 0.10468E+08, 0.10533E+08, 0.10599E+08, 0.10665E+08, 0.10731E+08, 0.10798E+08, 0.10865E+08, 0.10933E+08, 0.11000E+08 &
, 0.11068E+08, 0.11137E+08, 0.11205E+08, 0.11274E+08, 0.11344E+08, 0.11413E+08, 0.11483E+08, 0.11554E+08, 0.11625E+08, 0.11696E+08 &
, 0.11767E+08, 0.11839E+08, 0.11911E+08, 0.11983E+08, 0.12056E+08, 0.12129E+08, 0.12203E+08, 0.12277E+08, 0.12351E+08, 0.12426E+08 &
, 0.12501E+08, 0.12576E+08, 0.12652E+08, 0.12728E+08, 0.12804E+08, 0.12881E+08, 0.12958E+08, 0.13035E+08, 0.13113E+08, 0.13191E+08 &
, 0.13270E+08, 0.13349E+08, 0.13428E+08, 0.13508E+08, 0.13588E+08, 0.13669E+08, 0.13749E+08, 0.13831E+08, 0.13912E+08, 0.13994E+08 /
 !            4           1
 DATA(QofT(          41 ,J),J=1,501)/   301.600006103516       &
, 0.33146E+03, 0.36132E+03, 0.39118E+03, 0.42103E+03, 0.45089E+03, 0.48075E+03, 0.51061E+03, 0.54047E+03, 0.57033E+03, 0.60019E+03 &
, 0.63005E+03, 0.65991E+03, 0.68977E+03, 0.71964E+03, 0.74950E+03, 0.77936E+03, 0.80922E+03, 0.83908E+03, 0.86894E+03, 0.89881E+03 &
, 0.92867E+03, 0.95853E+03, 0.98840E+03, 0.10183E+04, 0.10481E+04, 0.10780E+04, 0.11079E+04, 0.11377E+04, 0.11676E+04, 0.11975E+04 &
, 0.12274E+04, 0.12573E+04, 0.12872E+04, 0.13171E+04, 0.13470E+04, 0.13769E+04, 0.14068E+04, 0.14368E+04, 0.14667E+04, 0.14967E+04 &
, 0.15267E+04, 0.15567E+04, 0.15868E+04, 0.16168E+04, 0.16469E+04, 0.16770E+04, 0.17072E+04, 0.17373E+04, 0.17676E+04, 0.17978E+04 &
, 0.18281E+04, 0.18585E+04, 0.18889E+04, 0.19193E+04, 0.19498E+04, 0.19804E+04, 0.20110E+04, 0.20417E+04, 0.20725E+04, 0.21033E+04 &
, 0.21342E+04, 0.21652E+04, 0.21962E+04, 0.22274E+04, 0.22586E+04, 0.22899E+04, 0.23214E+04, 0.23529E+04, 0.23845E+04, 0.24162E+04 &
, 0.24481E+04, 0.24800E+04, 0.25121E+04, 0.25443E+04, 0.25766E+04, 0.26090E+04, 0.26416E+04, 0.26743E+04, 0.27071E+04, 0.27400E+04 &
, 0.27732E+04, 0.28064E+04, 0.28398E+04, 0.28733E+04, 0.29071E+04, 0.29409E+04, 0.29749E+04, 0.30091E+04, 0.30435E+04, 0.30780E+04 &
, 0.31127E+04, 0.31475E+04, 0.31826E+04, 0.32178E+04, 0.32532E+04, 0.32888E+04, 0.33246E+04, 0.33605E+04, 0.33967E+04, 0.34331E+04 &
, 0.34696E+04, 0.35064E+04, 0.35434E+04, 0.35805E+04, 0.36179E+04, 0.36555E+04, 0.36933E+04, 0.37314E+04, 0.37696E+04, 0.38081E+04 &
, 0.38468E+04, 0.38857E+04, 0.39248E+04, 0.39642E+04, 0.40038E+04, 0.40437E+04, 0.40838E+04, 0.41241E+04, 0.41647E+04, 0.42055E+04 &
, 0.42466E+04, 0.42879E+04, 0.43295E+04, 0.43713E+04, 0.44134E+04, 0.44557E+04, 0.44983E+04, 0.45412E+04, 0.45844E+04, 0.46278E+04 &
, 0.46714E+04, 0.47154E+04, 0.47596E+04, 0.48041E+04, 0.48489E+04, 0.48939E+04, 0.49393E+04, 0.49849E+04, 0.50308E+04, 0.50770E+04 &
, 0.51235E+04, 0.51703E+04, 0.52174E+04, 0.52648E+04, 0.53124E+04, 0.53604E+04, 0.54087E+04, 0.54573E+04, 0.55062E+04, 0.55554E+04 &
, 0.56049E+04, 0.56548E+04, 0.57049E+04, 0.57554E+04, 0.58062E+04, 0.58573E+04, 0.59088E+04, 0.59605E+04, 0.60126E+04, 0.60650E+04 &
, 0.61178E+04, 0.61709E+04, 0.62243E+04, 0.62781E+04, 0.63322E+04, 0.63867E+04, 0.64415E+04, 0.64966E+04, 0.65521E+04, 0.66080E+04 &
, 0.66642E+04, 0.67208E+04, 0.67777E+04, 0.68350E+04, 0.68926E+04, 0.69506E+04, 0.70090E+04, 0.70677E+04, 0.71268E+04, 0.71863E+04 &
, 0.72462E+04, 0.73064E+04, 0.73670E+04, 0.74280E+04, 0.74894E+04, 0.75512E+04, 0.76133E+04, 0.76758E+04, 0.77388E+04, 0.78021E+04 &
, 0.78658E+04, 0.79299E+04, 0.79945E+04, 0.80594E+04, 0.81247E+04, 0.81905E+04, 0.82566E+04, 0.83232E+04, 0.83901E+04, 0.84575E+04 &
, 0.85253E+04, 0.85935E+04, 0.86622E+04, 0.87313E+04, 0.88007E+04, 0.88707E+04, 0.89410E+04, 0.90118E+04, 0.90830E+04, 0.91547E+04 &
, 0.92268E+04, 0.92994E+04, 0.93723E+04, 0.94458E+04, 0.95197E+04, 0.95940E+04, 0.96688E+04, 0.97440E+04, 0.98198E+04, 0.98959E+04 &
, 0.99726E+04, 0.10050E+05, 0.10127E+05, 0.10205E+05, 0.10284E+05, 0.10363E+05, 0.10442E+05, 0.10522E+05, 0.10603E+05, 0.10683E+05 &
, 0.10765E+05, 0.10847E+05, 0.10929E+05, 0.11012E+05, 0.11095E+05, 0.11179E+05, 0.11264E+05, 0.11348E+05, 0.11434E+05, 0.11520E+05 &
, 0.11606E+05, 0.11693E+05, 0.11780E+05, 0.11868E+05, 0.11957E+05, 0.12046E+05, 0.12135E+05, 0.12225E+05, 0.12316E+05, 0.12407E+05 &
, 0.12499E+05, 0.12591E+05, 0.12684E+05, 0.12777E+05, 0.12871E+05, 0.12965E+05, 0.13060E+05, 0.13155E+05, 0.13251E+05, 0.13348E+05 &
, 0.13445E+05, 0.13543E+05, 0.13641E+05, 0.13740E+05, 0.13839E+05, 0.13939E+05, 0.14040E+05, 0.14141E+05, 0.14243E+05, 0.14345E+05 &
, 0.14448E+05, 0.14552E+05, 0.14656E+05, 0.14760E+05, 0.14866E+05, 0.14972E+05, 0.15078E+05, 0.15185E+05, 0.15293E+05, 0.15401E+05 &
, 0.15510E+05, 0.15620E+05, 0.15730E+05, 0.15841E+05, 0.15952E+05, 0.16064E+05, 0.16177E+05, 0.16290E+05, 0.16404E+05, 0.16519E+05 &
, 0.16634E+05, 0.16750E+05, 0.16867E+05, 0.16984E+05, 0.17102E+05, 0.17220E+05, 0.17340E+05, 0.17459E+05, 0.17580E+05, 0.17701E+05 &
, 0.17823E+05, 0.17946E+05, 0.18069E+05, 0.18193E+05, 0.18317E+05, 0.18443E+05, 0.18569E+05, 0.18695E+05, 0.18823E+05, 0.18951E+05 &
, 0.19080E+05, 0.19209E+05, 0.19339E+05, 0.19470E+05, 0.19602E+05, 0.19734E+05, 0.19867E+05, 0.20001E+05, 0.20136E+05, 0.20271E+05 &
, 0.20407E+05, 0.20544E+05, 0.20681E+05, 0.20820E+05, 0.20959E+05, 0.21098E+05, 0.21239E+05, 0.21380E+05, 0.21522E+05, 0.21665E+05 &
, 0.21808E+05, 0.21953E+05, 0.22098E+05, 0.22244E+05, 0.22390E+05, 0.22538E+05, 0.22686E+05, 0.22835E+05, 0.22985E+05, 0.23135E+05 &
, 0.23287E+05, 0.23439E+05, 0.23592E+05, 0.23746E+05, 0.23900E+05, 0.24056E+05, 0.24212E+05, 0.24369E+05, 0.24527E+05, 0.24686E+05 &
, 0.24846E+05, 0.25006E+05, 0.25167E+05, 0.25329E+05, 0.25492E+05, 0.25656E+05, 0.25821E+05, 0.25987E+05, 0.26153E+05, 0.26320E+05 &
, 0.26488E+05, 0.26657E+05, 0.26827E+05, 0.26998E+05, 0.27170E+05, 0.27342E+05, 0.27516E+05, 0.27690E+05, 0.27865E+05, 0.28041E+05 &
, 0.28219E+05, 0.28397E+05, 0.28575E+05, 0.28755E+05, 0.28936E+05, 0.29118E+05, 0.29300E+05, 0.29484E+05, 0.29668E+05, 0.29853E+05 &
, 0.30040E+05, 0.30227E+05, 0.30415E+05, 0.30604E+05, 0.30795E+05, 0.30986E+05, 0.31178E+05, 0.31371E+05, 0.31565E+05, 0.31760E+05 &
, 0.31956E+05, 0.32153E+05, 0.32351E+05, 0.32550E+05, 0.32750E+05, 0.32950E+05, 0.33152E+05, 0.33355E+05, 0.33559E+05, 0.33764E+05 &
, 0.33970E+05, 0.34177E+05, 0.34385E+05, 0.34595E+05, 0.34805E+05, 0.35016E+05, 0.35228E+05, 0.35441E+05, 0.35656E+05, 0.35871E+05 &
, 0.36087E+05, 0.36305E+05, 0.36524E+05, 0.36743E+05, 0.36964E+05, 0.37186E+05, 0.37409E+05, 0.37633E+05, 0.37858E+05, 0.38084E+05 &
, 0.38311E+05, 0.38540E+05, 0.38769E+05, 0.39000E+05, 0.39231E+05, 0.39464E+05, 0.39698E+05, 0.39933E+05, 0.40170E+05, 0.40407E+05 &
, 0.40646E+05, 0.40885E+05, 0.41126E+05, 0.41368E+05, 0.41611E+05, 0.41856E+05, 0.42101E+05, 0.42348E+05, 0.42596E+05, 0.42845E+05 &
, 0.43095E+05, 0.43347E+05, 0.43599E+05, 0.43853E+05, 0.44108E+05, 0.44365E+05, 0.44622E+05, 0.44881E+05, 0.45141E+05, 0.45402E+05 &
, 0.45664E+05, 0.45928E+05, 0.46193E+05, 0.46459E+05, 0.46726E+05, 0.46995E+05, 0.47265E+05, 0.47536E+05, 0.47808E+05, 0.48082E+05 &
, 0.48357E+05, 0.48633E+05, 0.48911E+05, 0.49190E+05, 0.49470E+05, 0.49752E+05, 0.50034E+05, 0.50318E+05, 0.50604E+05, 0.50891E+05 &
, 0.51179E+05, 0.51468E+05, 0.51759E+05, 0.52051E+05, 0.52344E+05, 0.52639E+05, 0.52935E+05, 0.53233E+05, 0.53531E+05, 0.53832E+05 &
, 0.54133E+05, 0.54436E+05, 0.54741E+05, 0.55046E+05, 0.55354E+05, 0.55662E+05, 0.55972E+05, 0.56283E+05, 0.56596E+05, 0.56910E+05 &
, 0.57226E+05, 0.57543E+05, 0.57862E+05, 0.58182E+05, 0.58503E+05, 0.58826E+05, 0.59150E+05, 0.59476E+05, 0.59803E+05, 0.60132E+05 /
 !            4           2
 DATA(QofT(          42 ,J),J=1,501)/   201.139999389648       &
, 0.22106E+03, 0.24098E+03, 0.26089E+03, 0.28081E+03, 0.30073E+03, 0.32066E+03, 0.34058E+03, 0.36050E+03, 0.38043E+03, 0.40035E+03 &
, 0.42028E+03, 0.44021E+03, 0.46014E+03, 0.48007E+03, 0.50000E+03, 0.51993E+03, 0.53986E+03, 0.55979E+03, 0.57973E+03, 0.59966E+03 &
, 0.61960E+03, 0.63953E+03, 0.65948E+03, 0.67942E+03, 0.69936E+03, 0.71930E+03, 0.73925E+03, 0.75920E+03, 0.77916E+03, 0.79911E+03 &
, 0.81907E+03, 0.83904E+03, 0.85901E+03, 0.87899E+03, 0.89898E+03, 0.91897E+03, 0.93898E+03, 0.95899E+03, 0.97902E+03, 0.99906E+03 &
, 0.10191E+04, 0.10392E+04, 0.10593E+04, 0.10794E+04, 0.10995E+04, 0.11197E+04, 0.11399E+04, 0.11601E+04, 0.11803E+04, 0.12006E+04 &
, 0.12209E+04, 0.12412E+04, 0.12616E+04, 0.12820E+04, 0.13024E+04, 0.13229E+04, 0.13435E+04, 0.13641E+04, 0.13847E+04, 0.14054E+04 &
, 0.14261E+04, 0.14470E+04, 0.14678E+04, 0.14888E+04, 0.15098E+04, 0.15308E+04, 0.15520E+04, 0.15732E+04, 0.15945E+04, 0.16158E+04 &
, 0.16373E+04, 0.16588E+04, 0.16804E+04, 0.17021E+04, 0.17239E+04, 0.17458E+04, 0.17677E+04, 0.17898E+04, 0.18120E+04, 0.18342E+04 &
, 0.18566E+04, 0.18790E+04, 0.19016E+04, 0.19243E+04, 0.19471E+04, 0.19700E+04, 0.19930E+04, 0.20162E+04, 0.20394E+04, 0.20628E+04 &
, 0.20863E+04, 0.21099E+04, 0.21337E+04, 0.21575E+04, 0.21815E+04, 0.22057E+04, 0.22300E+04, 0.22544E+04, 0.22789E+04, 0.23036E+04 &
, 0.23285E+04, 0.23534E+04, 0.23785E+04, 0.24038E+04, 0.24292E+04, 0.24548E+04, 0.24805E+04, 0.25064E+04, 0.25324E+04, 0.25586E+04 &
, 0.25850E+04, 0.26115E+04, 0.26381E+04, 0.26650E+04, 0.26919E+04, 0.27191E+04, 0.27465E+04, 0.27740E+04, 0.28016E+04, 0.28295E+04 &
, 0.28575E+04, 0.28857E+04, 0.29140E+04, 0.29426E+04, 0.29713E+04, 0.30003E+04, 0.30294E+04, 0.30586E+04, 0.30881E+04, 0.31178E+04 &
, 0.31476E+04, 0.31777E+04, 0.32079E+04, 0.32383E+04, 0.32689E+04, 0.32998E+04, 0.33308E+04, 0.33620E+04, 0.33934E+04, 0.34251E+04 &
, 0.34569E+04, 0.34889E+04, 0.35211E+04, 0.35536E+04, 0.35862E+04, 0.36191E+04, 0.36522E+04, 0.36855E+04, 0.37190E+04, 0.37528E+04 &
, 0.37867E+04, 0.38209E+04, 0.38553E+04, 0.38899E+04, 0.39247E+04, 0.39598E+04, 0.39951E+04, 0.40306E+04, 0.40663E+04, 0.41024E+04 &
, 0.41386E+04, 0.41750E+04, 0.42117E+04, 0.42486E+04, 0.42858E+04, 0.43232E+04, 0.43608E+04, 0.43987E+04, 0.44368E+04, 0.44752E+04 &
, 0.45138E+04, 0.45527E+04, 0.45919E+04, 0.46312E+04, 0.46708E+04, 0.47107E+04, 0.47509E+04, 0.47913E+04, 0.48319E+04, 0.48728E+04 &
, 0.49140E+04, 0.49554E+04, 0.49972E+04, 0.50391E+04, 0.50814E+04, 0.51239E+04, 0.51667E+04, 0.52097E+04, 0.52531E+04, 0.52966E+04 &
, 0.53405E+04, 0.53847E+04, 0.54292E+04, 0.54738E+04, 0.55189E+04, 0.55641E+04, 0.56097E+04, 0.56556E+04, 0.57017E+04, 0.57482E+04 &
, 0.57949E+04, 0.58419E+04, 0.58893E+04, 0.59369E+04, 0.59848E+04, 0.60330E+04, 0.60815E+04, 0.61304E+04, 0.61795E+04, 0.62289E+04 &
, 0.62786E+04, 0.63287E+04, 0.63791E+04, 0.64297E+04, 0.64807E+04, 0.65320E+04, 0.65836E+04, 0.66356E+04, 0.66879E+04, 0.67404E+04 &
, 0.67933E+04, 0.68466E+04, 0.69001E+04, 0.69540E+04, 0.70082E+04, 0.70628E+04, 0.71177E+04, 0.71729E+04, 0.72284E+04, 0.72843E+04 &
, 0.73406E+04, 0.73972E+04, 0.74541E+04, 0.75113E+04, 0.75690E+04, 0.76269E+04, 0.76852E+04, 0.77439E+04, 0.78029E+04, 0.78623E+04 &
, 0.79220E+04, 0.79821E+04, 0.80426E+04, 0.81034E+04, 0.81646E+04, 0.82262E+04, 0.82881E+04, 0.83503E+04, 0.84130E+04, 0.84761E+04 &
, 0.85395E+04, 0.86032E+04, 0.86674E+04, 0.87320E+04, 0.87969E+04, 0.88622E+04, 0.89279E+04, 0.89940E+04, 0.90605E+04, 0.91274E+04 &
, 0.91946E+04, 0.92622E+04, 0.93303E+04, 0.93988E+04, 0.94676E+04, 0.95369E+04, 0.96065E+04, 0.96766E+04, 0.97471E+04, 0.98180E+04 &
, 0.98893E+04, 0.99610E+04, 0.10033E+05, 0.10106E+05, 0.10179E+05, 0.10252E+05, 0.10326E+05, 0.10400E+05, 0.10475E+05, 0.10550E+05 &
, 0.10625E+05, 0.10701E+05, 0.10778E+05, 0.10855E+05, 0.10932E+05, 0.11010E+05, 0.11088E+05, 0.11166E+05, 0.11245E+05, 0.11325E+05 &
, 0.11405E+05, 0.11485E+05, 0.11566E+05, 0.11648E+05, 0.11729E+05, 0.11812E+05, 0.11894E+05, 0.11978E+05, 0.12061E+05, 0.12145E+05 &
, 0.12230E+05, 0.12315E+05, 0.12401E+05, 0.12487E+05, 0.12573E+05, 0.12660E+05, 0.12747E+05, 0.12835E+05, 0.12924E+05, 0.13013E+05 &
, 0.13102E+05, 0.13192E+05, 0.13283E+05, 0.13374E+05, 0.13465E+05, 0.13557E+05, 0.13649E+05, 0.13742E+05, 0.13836E+05, 0.13930E+05 &
, 0.14024E+05, 0.14119E+05, 0.14215E+05, 0.14311E+05, 0.14407E+05, 0.14505E+05, 0.14602E+05, 0.14700E+05, 0.14799E+05, 0.14898E+05 &
, 0.14998E+05, 0.15098E+05, 0.15199E+05, 0.15301E+05, 0.15403E+05, 0.15505E+05, 0.15608E+05, 0.15712E+05, 0.15816E+05, 0.15921E+05 &
, 0.16026E+05, 0.16132E+05, 0.16238E+05, 0.16345E+05, 0.16453E+05, 0.16561E+05, 0.16670E+05, 0.16779E+05, 0.16889E+05, 0.16999E+05 &
, 0.17110E+05, 0.17222E+05, 0.17334E+05, 0.17447E+05, 0.17560E+05, 0.17674E+05, 0.17789E+05, 0.17904E+05, 0.18020E+05, 0.18136E+05 &
, 0.18253E+05, 0.18371E+05, 0.18489E+05, 0.18608E+05, 0.18728E+05, 0.18848E+05, 0.18969E+05, 0.19090E+05, 0.19212E+05, 0.19335E+05 &
, 0.19458E+05, 0.19582E+05, 0.19706E+05, 0.19832E+05, 0.19957E+05, 0.20084E+05, 0.20211E+05, 0.20339E+05, 0.20467E+05, 0.20596E+05 &
, 0.20726E+05, 0.20857E+05, 0.20988E+05, 0.21120E+05, 0.21252E+05, 0.21385E+05, 0.21519E+05, 0.21653E+05, 0.21789E+05, 0.21924E+05 &
, 0.22061E+05, 0.22198E+05, 0.22336E+05, 0.22475E+05, 0.22614E+05, 0.22754E+05, 0.22895E+05, 0.23036E+05, 0.23178E+05, 0.23321E+05 &
, 0.23465E+05, 0.23609E+05, 0.23754E+05, 0.23900E+05, 0.24046E+05, 0.24194E+05, 0.24342E+05, 0.24490E+05, 0.24640E+05, 0.24790E+05 &
, 0.24941E+05, 0.25092E+05, 0.25245E+05, 0.25398E+05, 0.25552E+05, 0.25706E+05, 0.25862E+05, 0.26018E+05, 0.26175E+05, 0.26333E+05 &
, 0.26491E+05, 0.26651E+05, 0.26811E+05, 0.26971E+05, 0.27133E+05, 0.27296E+05, 0.27459E+05, 0.27623E+05, 0.27788E+05, 0.27953E+05 &
, 0.28120E+05, 0.28287E+05, 0.28455E+05, 0.28624E+05, 0.28793E+05, 0.28964E+05, 0.29135E+05, 0.29307E+05, 0.29480E+05, 0.29654E+05 &
, 0.29828E+05, 0.30004E+05, 0.30180E+05, 0.30357E+05, 0.30535E+05, 0.30714E+05, 0.30894E+05, 0.31074E+05, 0.31256E+05, 0.31438E+05 &
, 0.31621E+05, 0.31805E+05, 0.31990E+05, 0.32176E+05, 0.32362E+05, 0.32550E+05, 0.32738E+05, 0.32927E+05, 0.33118E+05, 0.33309E+05 &
, 0.33501E+05, 0.33693E+05, 0.33887E+05, 0.34082E+05, 0.34278E+05, 0.34474E+05, 0.34671E+05, 0.34870E+05, 0.35069E+05, 0.35269E+05 &
, 0.35470E+05, 0.35672E+05, 0.35875E+05, 0.36079E+05, 0.36284E+05, 0.36490E+05, 0.36697E+05, 0.36904E+05, 0.37113E+05, 0.37323E+05 &
, 0.37533E+05, 0.37745E+05, 0.37957E+05, 0.38171E+05, 0.38385E+05, 0.38601E+05, 0.38817E+05, 0.39035E+05, 0.39253E+05, 0.39473E+05 &
, 0.39693E+05, 0.39915E+05, 0.40137E+05, 0.40361E+05, 0.40585E+05, 0.40811E+05, 0.41037E+05, 0.41265E+05, 0.41493E+05, 0.41723E+05 /
 !            4           3
 DATA(QofT(          43 ,J),J=1,501)/   208.089996337891       &
, 0.22870E+03, 0.24931E+03, 0.26992E+03, 0.29054E+03, 0.31116E+03, 0.33177E+03, 0.35239E+03, 0.37301E+03, 0.39363E+03, 0.41425E+03 &
, 0.43487E+03, 0.45549E+03, 0.47612E+03, 0.49674E+03, 0.51737E+03, 0.53800E+03, 0.55862E+03, 0.57925E+03, 0.59988E+03, 0.62051E+03 &
, 0.64114E+03, 0.66177E+03, 0.68241E+03, 0.70305E+03, 0.72368E+03, 0.74433E+03, 0.76497E+03, 0.78561E+03, 0.80626E+03, 0.82691E+03 &
, 0.84757E+03, 0.86823E+03, 0.88889E+03, 0.90956E+03, 0.93024E+03, 0.95093E+03, 0.97162E+03, 0.99233E+03, 0.10130E+04, 0.10338E+04 &
, 0.10545E+04, 0.10753E+04, 0.10960E+04, 0.11168E+04, 0.11377E+04, 0.11585E+04, 0.11793E+04, 0.12002E+04, 0.12211E+04, 0.12421E+04 &
, 0.12631E+04, 0.12841E+04, 0.13051E+04, 0.13262E+04, 0.13473E+04, 0.13685E+04, 0.13897E+04, 0.14109E+04, 0.14322E+04, 0.14536E+04 &
, 0.14750E+04, 0.14965E+04, 0.15180E+04, 0.15396E+04, 0.15612E+04, 0.15830E+04, 0.16047E+04, 0.16266E+04, 0.16485E+04, 0.16705E+04 &
, 0.16926E+04, 0.17148E+04, 0.17370E+04, 0.17593E+04, 0.17818E+04, 0.18043E+04, 0.18269E+04, 0.18495E+04, 0.18723E+04, 0.18952E+04 &
, 0.19182E+04, 0.19413E+04, 0.19645E+04, 0.19878E+04, 0.20112E+04, 0.20347E+04, 0.20584E+04, 0.20821E+04, 0.21060E+04, 0.21300E+04 &
, 0.21541E+04, 0.21783E+04, 0.22027E+04, 0.22272E+04, 0.22518E+04, 0.22766E+04, 0.23014E+04, 0.23265E+04, 0.23516E+04, 0.23769E+04 &
, 0.24024E+04, 0.24279E+04, 0.24537E+04, 0.24795E+04, 0.25056E+04, 0.25317E+04, 0.25581E+04, 0.25845E+04, 0.26112E+04, 0.26380E+04 &
, 0.26649E+04, 0.26920E+04, 0.27193E+04, 0.27468E+04, 0.27744E+04, 0.28021E+04, 0.28301E+04, 0.28582E+04, 0.28865E+04, 0.29149E+04 &
, 0.29435E+04, 0.29723E+04, 0.30013E+04, 0.30305E+04, 0.30598E+04, 0.30894E+04, 0.31191E+04, 0.31490E+04, 0.31791E+04, 0.32094E+04 &
, 0.32398E+04, 0.32705E+04, 0.33013E+04, 0.33324E+04, 0.33636E+04, 0.33951E+04, 0.34267E+04, 0.34586E+04, 0.34906E+04, 0.35229E+04 &
, 0.35553E+04, 0.35880E+04, 0.36209E+04, 0.36540E+04, 0.36873E+04, 0.37208E+04, 0.37545E+04, 0.37885E+04, 0.38226E+04, 0.38570E+04 &
, 0.38916E+04, 0.39264E+04, 0.39615E+04, 0.39968E+04, 0.40322E+04, 0.40680E+04, 0.41039E+04, 0.41401E+04, 0.41765E+04, 0.42132E+04 &
, 0.42501E+04, 0.42872E+04, 0.43246E+04, 0.43622E+04, 0.44000E+04, 0.44381E+04, 0.44764E+04, 0.45150E+04, 0.45538E+04, 0.45929E+04 &
, 0.46322E+04, 0.46718E+04, 0.47116E+04, 0.47517E+04, 0.47920E+04, 0.48326E+04, 0.48734E+04, 0.49146E+04, 0.49559E+04, 0.49976E+04 &
, 0.50394E+04, 0.50816E+04, 0.51241E+04, 0.51668E+04, 0.52098E+04, 0.52530E+04, 0.52965E+04, 0.53403E+04, 0.53844E+04, 0.54287E+04 &
, 0.54733E+04, 0.55183E+04, 0.55635E+04, 0.56089E+04, 0.56547E+04, 0.57008E+04, 0.57471E+04, 0.57937E+04, 0.58406E+04, 0.58879E+04 &
, 0.59354E+04, 0.59832E+04, 0.60313E+04, 0.60797E+04, 0.61284E+04, 0.61774E+04, 0.62267E+04, 0.62764E+04, 0.63263E+04, 0.63766E+04 &
, 0.64271E+04, 0.64780E+04, 0.65291E+04, 0.65807E+04, 0.66324E+04, 0.66846E+04, 0.67370E+04, 0.67898E+04, 0.68429E+04, 0.68963E+04 &
, 0.69501E+04, 0.70042E+04, 0.70586E+04, 0.71133E+04, 0.71684E+04, 0.72238E+04, 0.72796E+04, 0.73357E+04, 0.73921E+04, 0.74489E+04 &
, 0.75060E+04, 0.75635E+04, 0.76213E+04, 0.76794E+04, 0.77379E+04, 0.77968E+04, 0.78560E+04, 0.79156E+04, 0.79755E+04, 0.80358E+04 &
, 0.80965E+04, 0.81575E+04, 0.82189E+04, 0.82807E+04, 0.83428E+04, 0.84052E+04, 0.84681E+04, 0.85313E+04, 0.85950E+04, 0.86590E+04 &
, 0.87233E+04, 0.87881E+04, 0.88532E+04, 0.89187E+04, 0.89846E+04, 0.90509E+04, 0.91176E+04, 0.91847E+04, 0.92521E+04, 0.93200E+04 &
, 0.93883E+04, 0.94570E+04, 0.95260E+04, 0.95955E+04, 0.96653E+04, 0.97356E+04, 0.98063E+04, 0.98774E+04, 0.99489E+04, 0.10021E+05 &
, 0.10093E+05, 0.10166E+05, 0.10239E+05, 0.10313E+05, 0.10387E+05, 0.10461E+05, 0.10536E+05, 0.10611E+05, 0.10687E+05, 0.10763E+05 &
, 0.10840E+05, 0.10917E+05, 0.10994E+05, 0.11072E+05, 0.11151E+05, 0.11230E+05, 0.11309E+05, 0.11389E+05, 0.11469E+05, 0.11549E+05 &
, 0.11631E+05, 0.11712E+05, 0.11794E+05, 0.11877E+05, 0.11960E+05, 0.12043E+05, 0.12127E+05, 0.12211E+05, 0.12296E+05, 0.12381E+05 &
, 0.12467E+05, 0.12553E+05, 0.12640E+05, 0.12727E+05, 0.12815E+05, 0.12903E+05, 0.12992E+05, 0.13081E+05, 0.13171E+05, 0.13261E+05 &
, 0.13351E+05, 0.13443E+05, 0.13534E+05, 0.13626E+05, 0.13719E+05, 0.13812E+05, 0.13906E+05, 0.14000E+05, 0.14095E+05, 0.14190E+05 &
, 0.14286E+05, 0.14382E+05, 0.14479E+05, 0.14576E+05, 0.14674E+05, 0.14772E+05, 0.14871E+05, 0.14971E+05, 0.15071E+05, 0.15171E+05 &
, 0.15273E+05, 0.15374E+05, 0.15476E+05, 0.15579E+05, 0.15682E+05, 0.15786E+05, 0.15891E+05, 0.15996E+05, 0.16101E+05, 0.16207E+05 &
, 0.16314E+05, 0.16421E+05, 0.16529E+05, 0.16637E+05, 0.16746E+05, 0.16856E+05, 0.16966E+05, 0.17077E+05, 0.17188E+05, 0.17300E+05 &
, 0.17412E+05, 0.17525E+05, 0.17639E+05, 0.17753E+05, 0.17868E+05, 0.17983E+05, 0.18099E+05, 0.18216E+05, 0.18333E+05, 0.18451E+05 &
, 0.18570E+05, 0.18689E+05, 0.18809E+05, 0.18929E+05, 0.19050E+05, 0.19172E+05, 0.19294E+05, 0.19417E+05, 0.19540E+05, 0.19664E+05 &
, 0.19789E+05, 0.19915E+05, 0.20041E+05, 0.20167E+05, 0.20295E+05, 0.20423E+05, 0.20552E+05, 0.20681E+05, 0.20811E+05, 0.20942E+05 &
, 0.21073E+05, 0.21205E+05, 0.21338E+05, 0.21471E+05, 0.21605E+05, 0.21740E+05, 0.21875E+05, 0.22011E+05, 0.22148E+05, 0.22286E+05 &
, 0.22424E+05, 0.22563E+05, 0.22702E+05, 0.22842E+05, 0.22984E+05, 0.23125E+05, 0.23268E+05, 0.23411E+05, 0.23554E+05, 0.23699E+05 &
, 0.23844E+05, 0.23990E+05, 0.24137E+05, 0.24285E+05, 0.24433E+05, 0.24582E+05, 0.24731E+05, 0.24882E+05, 0.25033E+05, 0.25185E+05 &
, 0.25337E+05, 0.25491E+05, 0.25645E+05, 0.25800E+05, 0.25956E+05, 0.26112E+05, 0.26269E+05, 0.26427E+05, 0.26586E+05, 0.26745E+05 &
, 0.26906E+05, 0.27067E+05, 0.27229E+05, 0.27391E+05, 0.27555E+05, 0.27719E+05, 0.27884E+05, 0.28050E+05, 0.28217E+05, 0.28384E+05 &
, 0.28552E+05, 0.28721E+05, 0.28891E+05, 0.29062E+05, 0.29234E+05, 0.29406E+05, 0.29579E+05, 0.29753E+05, 0.29928E+05, 0.30104E+05 &
, 0.30280E+05, 0.30458E+05, 0.30636E+05, 0.30815E+05, 0.30995E+05, 0.31176E+05, 0.31357E+05, 0.31540E+05, 0.31723E+05, 0.31907E+05 &
, 0.32093E+05, 0.32279E+05, 0.32465E+05, 0.32653E+05, 0.32842E+05, 0.33031E+05, 0.33222E+05, 0.33413E+05, 0.33605E+05, 0.33798E+05 &
, 0.33992E+05, 0.34187E+05, 0.34383E+05, 0.34580E+05, 0.34777E+05, 0.34976E+05, 0.35176E+05, 0.35376E+05, 0.35577E+05, 0.35780E+05 &
, 0.35983E+05, 0.36187E+05, 0.36392E+05, 0.36598E+05, 0.36805E+05, 0.37013E+05, 0.37222E+05, 0.37432E+05, 0.37643E+05, 0.37855E+05 &
, 0.38068E+05, 0.38281E+05, 0.38496E+05, 0.38712E+05, 0.38928E+05, 0.39146E+05, 0.39365E+05, 0.39585E+05, 0.39805E+05, 0.40027E+05 &
, 0.40250E+05, 0.40473E+05, 0.40698E+05, 0.40924E+05, 0.41151E+05, 0.41378E+05, 0.41607E+05, 0.41837E+05, 0.42068E+05, 0.42300E+05 /
 !            4           4
 DATA(QofT(          44 ,J),J=1,501)/   319.380004882812       &
, 0.35102E+03, 0.38267E+03, 0.41432E+03, 0.44596E+03, 0.47761E+03, 0.50926E+03, 0.54092E+03, 0.57257E+03, 0.60422E+03, 0.63588E+03 &
, 0.66754E+03, 0.69920E+03, 0.73086E+03, 0.76252E+03, 0.79419E+03, 0.82585E+03, 0.85752E+03, 0.88919E+03, 0.92086E+03, 0.95253E+03 &
, 0.98421E+03, 0.10159E+04, 0.10476E+04, 0.10792E+04, 0.11109E+04, 0.11426E+04, 0.11743E+04, 0.12060E+04, 0.12377E+04, 0.12694E+04 &
, 0.13011E+04, 0.13328E+04, 0.13645E+04, 0.13963E+04, 0.14280E+04, 0.14598E+04, 0.14916E+04, 0.15234E+04, 0.15552E+04, 0.15870E+04 &
, 0.16188E+04, 0.16507E+04, 0.16826E+04, 0.17145E+04, 0.17465E+04, 0.17785E+04, 0.18105E+04, 0.18426E+04, 0.18747E+04, 0.19068E+04 &
, 0.19390E+04, 0.19713E+04, 0.20036E+04, 0.20360E+04, 0.20684E+04, 0.21009E+04, 0.21335E+04, 0.21661E+04, 0.21988E+04, 0.22316E+04 &
, 0.22645E+04, 0.22975E+04, 0.23305E+04, 0.23637E+04, 0.23969E+04, 0.24303E+04, 0.24637E+04, 0.24973E+04, 0.25310E+04, 0.25648E+04 &
, 0.25987E+04, 0.26328E+04, 0.26669E+04, 0.27012E+04, 0.27357E+04, 0.27703E+04, 0.28050E+04, 0.28398E+04, 0.28749E+04, 0.29100E+04 &
, 0.29453E+04, 0.29808E+04, 0.30164E+04, 0.30523E+04, 0.30882E+04, 0.31244E+04, 0.31607E+04, 0.31972E+04, 0.32339E+04, 0.32708E+04 &
, 0.33079E+04, 0.33451E+04, 0.33826E+04, 0.34202E+04, 0.34580E+04, 0.34961E+04, 0.35344E+04, 0.35728E+04, 0.36115E+04, 0.36504E+04 &
, 0.36895E+04, 0.37289E+04, 0.37684E+04, 0.38082E+04, 0.38482E+04, 0.38885E+04, 0.39290E+04, 0.39697E+04, 0.40106E+04, 0.40519E+04 &
, 0.40933E+04, 0.41350E+04, 0.41770E+04, 0.42192E+04, 0.42616E+04, 0.43044E+04, 0.43473E+04, 0.43906E+04, 0.44341E+04, 0.44779E+04 &
, 0.45219E+04, 0.45662E+04, 0.46109E+04, 0.46558E+04, 0.47009E+04, 0.47464E+04, 0.47921E+04, 0.48381E+04, 0.48844E+04, 0.49310E+04 &
, 0.49779E+04, 0.50251E+04, 0.50726E+04, 0.51204E+04, 0.51685E+04, 0.52169E+04, 0.52657E+04, 0.53147E+04, 0.53641E+04, 0.54138E+04 &
, 0.54637E+04, 0.55141E+04, 0.55647E+04, 0.56156E+04, 0.56669E+04, 0.57185E+04, 0.57705E+04, 0.58228E+04, 0.58754E+04, 0.59284E+04 &
, 0.59817E+04, 0.60353E+04, 0.60893E+04, 0.61436E+04, 0.61984E+04, 0.62534E+04, 0.63088E+04, 0.63646E+04, 0.64207E+04, 0.64772E+04 &
, 0.65340E+04, 0.65912E+04, 0.66489E+04, 0.67068E+04, 0.67651E+04, 0.68238E+04, 0.68829E+04, 0.69424E+04, 0.70022E+04, 0.70625E+04 &
, 0.71231E+04, 0.71841E+04, 0.72455E+04, 0.73073E+04, 0.73695E+04, 0.74321E+04, 0.74951E+04, 0.75585E+04, 0.76223E+04, 0.76865E+04 &
, 0.77511E+04, 0.78162E+04, 0.78816E+04, 0.79475E+04, 0.80138E+04, 0.80805E+04, 0.81476E+04, 0.82152E+04, 0.82832E+04, 0.83516E+04 &
, 0.84205E+04, 0.84898E+04, 0.85595E+04, 0.86297E+04, 0.87003E+04, 0.87714E+04, 0.88429E+04, 0.89149E+04, 0.89873E+04, 0.90602E+04 &
, 0.91335E+04, 0.92073E+04, 0.92816E+04, 0.93563E+04, 0.94315E+04, 0.95072E+04, 0.95833E+04, 0.96599E+04, 0.97370E+04, 0.98146E+04 &
, 0.98926E+04, 0.99711E+04, 0.10050E+05, 0.10130E+05, 0.10210E+05, 0.10290E+05, 0.10371E+05, 0.10453E+05, 0.10535E+05, 0.10617E+05 &
, 0.10700E+05, 0.10784E+05, 0.10868E+05, 0.10952E+05, 0.11037E+05, 0.11123E+05, 0.11209E+05, 0.11296E+05, 0.11383E+05, 0.11471E+05 &
, 0.11559E+05, 0.11648E+05, 0.11737E+05, 0.11827E+05, 0.11917E+05, 0.12008E+05, 0.12100E+05, 0.12192E+05, 0.12285E+05, 0.12378E+05 &
, 0.12471E+05, 0.12566E+05, 0.12661E+05, 0.12756E+05, 0.12852E+05, 0.12949E+05, 0.13046E+05, 0.13144E+05, 0.13242E+05, 0.13341E+05 &
, 0.13440E+05, 0.13540E+05, 0.13641E+05, 0.13742E+05, 0.13844E+05, 0.13947E+05, 0.14050E+05, 0.14154E+05, 0.14258E+05, 0.14363E+05 &
, 0.14468E+05, 0.14575E+05, 0.14681E+05, 0.14789E+05, 0.14897E+05, 0.15005E+05, 0.15115E+05, 0.15225E+05, 0.15335E+05, 0.15447E+05 &
, 0.15558E+05, 0.15671E+05, 0.15784E+05, 0.15898E+05, 0.16013E+05, 0.16128E+05, 0.16244E+05, 0.16360E+05, 0.16477E+05, 0.16595E+05 &
, 0.16714E+05, 0.16833E+05, 0.16953E+05, 0.17073E+05, 0.17194E+05, 0.17317E+05, 0.17439E+05, 0.17563E+05, 0.17687E+05, 0.17811E+05 &
, 0.17937E+05, 0.18063E+05, 0.18190E+05, 0.18317E+05, 0.18446E+05, 0.18575E+05, 0.18705E+05, 0.18835E+05, 0.18966E+05, 0.19098E+05 &
, 0.19231E+05, 0.19365E+05, 0.19499E+05, 0.19634E+05, 0.19769E+05, 0.19906E+05, 0.20043E+05, 0.20181E+05, 0.20320E+05, 0.20460E+05 &
, 0.20600E+05, 0.20741E+05, 0.20883E+05, 0.21026E+05, 0.21169E+05, 0.21313E+05, 0.21458E+05, 0.21604E+05, 0.21751E+05, 0.21898E+05 &
, 0.22046E+05, 0.22196E+05, 0.22345E+05, 0.22496E+05, 0.22648E+05, 0.22800E+05, 0.22953E+05, 0.23107E+05, 0.23262E+05, 0.23418E+05 &
, 0.23574E+05, 0.23732E+05, 0.23890E+05, 0.24049E+05, 0.24209E+05, 0.24370E+05, 0.24531E+05, 0.24694E+05, 0.24857E+05, 0.25022E+05 &
, 0.25187E+05, 0.25353E+05, 0.25520E+05, 0.25688E+05, 0.25856E+05, 0.26026E+05, 0.26196E+05, 0.26368E+05, 0.26540E+05, 0.26714E+05 &
, 0.26888E+05, 0.27063E+05, 0.27239E+05, 0.27416E+05, 0.27594E+05, 0.27773E+05, 0.27952E+05, 0.28133E+05, 0.28315E+05, 0.28497E+05 &
, 0.28681E+05, 0.28865E+05, 0.29051E+05, 0.29237E+05, 0.29425E+05, 0.29613E+05, 0.29802E+05, 0.29993E+05, 0.30184E+05, 0.30377E+05 &
, 0.30570E+05, 0.30764E+05, 0.30960E+05, 0.31156E+05, 0.31353E+05, 0.31552E+05, 0.31751E+05, 0.31952E+05, 0.32153E+05, 0.32356E+05 &
, 0.32559E+05, 0.32764E+05, 0.32969E+05, 0.33176E+05, 0.33384E+05, 0.33592E+05, 0.33802E+05, 0.34013E+05, 0.34225E+05, 0.34438E+05 &
, 0.34652E+05, 0.34867E+05, 0.35084E+05, 0.35301E+05, 0.35520E+05, 0.35739E+05, 0.35960E+05, 0.36182E+05, 0.36405E+05, 0.36629E+05 &
, 0.36854E+05, 0.37080E+05, 0.37307E+05, 0.37536E+05, 0.37766E+05, 0.37996E+05, 0.38228E+05, 0.38462E+05, 0.38696E+05, 0.38931E+05 &
, 0.39168E+05, 0.39406E+05, 0.39645E+05, 0.39885E+05, 0.40126E+05, 0.40369E+05, 0.40612E+05, 0.40857E+05, 0.41103E+05, 0.41351E+05 &
, 0.41599E+05, 0.41849E+05, 0.42100E+05, 0.42352E+05, 0.42605E+05, 0.42860E+05, 0.43116E+05, 0.43373E+05, 0.43631E+05, 0.43891E+05 &
, 0.44152E+05, 0.44414E+05, 0.44677E+05, 0.44942E+05, 0.45208E+05, 0.45475E+05, 0.45744E+05, 0.46014E+05, 0.46285E+05, 0.46557E+05 &
, 0.46831E+05, 0.47106E+05, 0.47382E+05, 0.47660E+05, 0.47939E+05, 0.48219E+05, 0.48501E+05, 0.48784E+05, 0.49068E+05, 0.49354E+05 &
, 0.49641E+05, 0.49929E+05, 0.50219E+05, 0.50510E+05, 0.50803E+05, 0.51097E+05, 0.51392E+05, 0.51689E+05, 0.51987E+05, 0.52286E+05 &
, 0.52587E+05, 0.52889E+05, 0.53193E+05, 0.53498E+05, 0.53805E+05, 0.54112E+05, 0.54422E+05, 0.54733E+05, 0.55045E+05, 0.55359E+05 &
, 0.55674E+05, 0.55991E+05, 0.56309E+05, 0.56628E+05, 0.56949E+05, 0.57272E+05, 0.57596E+05, 0.57921E+05, 0.58249E+05, 0.58577E+05 &
, 0.58907E+05, 0.59239E+05, 0.59572E+05, 0.59906E+05, 0.60242E+05, 0.60580E+05, 0.60919E+05, 0.61260E+05, 0.61602E+05, 0.61946E+05 &
, 0.62292E+05, 0.62639E+05, 0.62987E+05, 0.63337E+05, 0.63689E+05, 0.64043E+05, 0.64398E+05, 0.64754E+05, 0.65112E+05, 0.65472E+05 /
 !            4           5
 DATA(QofT(          45 ,J),J=1,501)/   1864.50000000000       &
, 0.20492E+04, 0.22339E+04, 0.24186E+04, 0.26033E+04, 0.27880E+04, 0.29727E+04, 0.31574E+04, 0.33422E+04, 0.35269E+04, 0.37117E+04 &
, 0.38965E+04, 0.40812E+04, 0.42660E+04, 0.44508E+04, 0.46356E+04, 0.48204E+04, 0.50053E+04, 0.51901E+04, 0.53749E+04, 0.55598E+04 &
, 0.57446E+04, 0.59295E+04, 0.61144E+04, 0.62993E+04, 0.64842E+04, 0.66692E+04, 0.68541E+04, 0.70390E+04, 0.72240E+04, 0.74091E+04 &
, 0.75942E+04, 0.77793E+04, 0.79644E+04, 0.81496E+04, 0.83349E+04, 0.85202E+04, 0.87056E+04, 0.88912E+04, 0.90768E+04, 0.92624E+04 &
, 0.94483E+04, 0.96343E+04, 0.98204E+04, 0.10007E+05, 0.10193E+05, 0.10380E+05, 0.10567E+05, 0.10754E+05, 0.10941E+05, 0.11129E+05 &
, 0.11317E+05, 0.11505E+05, 0.11693E+05, 0.11882E+05, 0.12071E+05, 0.12261E+05, 0.12451E+05, 0.12641E+05, 0.12832E+05, 0.13024E+05 &
, 0.13215E+05, 0.13408E+05, 0.13600E+05, 0.13794E+05, 0.13988E+05, 0.14182E+05, 0.14377E+05, 0.14573E+05, 0.14769E+05, 0.14966E+05 &
, 0.15164E+05, 0.15363E+05, 0.15562E+05, 0.15762E+05, 0.15962E+05, 0.16164E+05, 0.16366E+05, 0.16570E+05, 0.16774E+05, 0.16979E+05 &
, 0.17184E+05, 0.17391E+05, 0.17599E+05, 0.17807E+05, 0.18017E+05, 0.18228E+05, 0.18439E+05, 0.18652E+05, 0.18866E+05, 0.19081E+05 &
, 0.19296E+05, 0.19513E+05, 0.19731E+05, 0.19951E+05, 0.20171E+05, 0.20393E+05, 0.20616E+05, 0.20839E+05, 0.21065E+05, 0.21291E+05 &
, 0.21519E+05, 0.21748E+05, 0.21978E+05, 0.22210E+05, 0.22443E+05, 0.22677E+05, 0.22912E+05, 0.23150E+05, 0.23388E+05, 0.23628E+05 &
, 0.23869E+05, 0.24112E+05, 0.24356E+05, 0.24601E+05, 0.24848E+05, 0.25097E+05, 0.25347E+05, 0.25599E+05, 0.25852E+05, 0.26106E+05 &
, 0.26362E+05, 0.26620E+05, 0.26880E+05, 0.27141E+05, 0.27403E+05, 0.27668E+05, 0.27934E+05, 0.28201E+05, 0.28470E+05, 0.28741E+05 &
, 0.29014E+05, 0.29288E+05, 0.29565E+05, 0.29842E+05, 0.30122E+05, 0.30404E+05, 0.30687E+05, 0.30972E+05, 0.31258E+05, 0.31547E+05 &
, 0.31838E+05, 0.32130E+05, 0.32424E+05, 0.32720E+05, 0.33018E+05, 0.33318E+05, 0.33620E+05, 0.33924E+05, 0.34229E+05, 0.34537E+05 &
, 0.34847E+05, 0.35158E+05, 0.35472E+05, 0.35787E+05, 0.36105E+05, 0.36425E+05, 0.36746E+05, 0.37070E+05, 0.37396E+05, 0.37724E+05 &
, 0.38054E+05, 0.38386E+05, 0.38721E+05, 0.39057E+05, 0.39396E+05, 0.39736E+05, 0.40079E+05, 0.40425E+05, 0.40772E+05, 0.41121E+05 &
, 0.41473E+05, 0.41827E+05, 0.42184E+05, 0.42542E+05, 0.42903E+05, 0.43266E+05, 0.43632E+05, 0.44000E+05, 0.44370E+05, 0.44743E+05 &
, 0.45117E+05, 0.45495E+05, 0.45874E+05, 0.46256E+05, 0.46641E+05, 0.47028E+05, 0.47417E+05, 0.47809E+05, 0.48204E+05, 0.48600E+05 &
, 0.49000E+05, 0.49401E+05, 0.49806E+05, 0.50213E+05, 0.50622E+05, 0.51034E+05, 0.51449E+05, 0.51866E+05, 0.52286E+05, 0.52708E+05 &
, 0.53133E+05, 0.53561E+05, 0.53992E+05, 0.54425E+05, 0.54861E+05, 0.55299E+05, 0.55740E+05, 0.56184E+05, 0.56631E+05, 0.57081E+05 &
, 0.57533E+05, 0.57988E+05, 0.58446E+05, 0.58907E+05, 0.59370E+05, 0.59837E+05, 0.60306E+05, 0.60778E+05, 0.61253E+05, 0.61731E+05 &
, 0.62212E+05, 0.62696E+05, 0.63183E+05, 0.63672E+05, 0.64165E+05, 0.64661E+05, 0.65159E+05, 0.65661E+05, 0.66166E+05, 0.66674E+05 &
, 0.67185E+05, 0.67699E+05, 0.68216E+05, 0.68737E+05, 0.69260E+05, 0.69787E+05, 0.70317E+05, 0.70850E+05, 0.71386E+05, 0.71925E+05 &
, 0.72468E+05, 0.73014E+05, 0.73563E+05, 0.74115E+05, 0.74671E+05, 0.75230E+05, 0.75792E+05, 0.76358E+05, 0.76927E+05, 0.77499E+05 &
, 0.78075E+05, 0.78655E+05, 0.79237E+05, 0.79823E+05, 0.80413E+05, 0.81006E+05, 0.81602E+05, 0.82202E+05, 0.82806E+05, 0.83413E+05 &
, 0.84023E+05, 0.84638E+05, 0.85255E+05, 0.85877E+05, 0.86502E+05, 0.87130E+05, 0.87763E+05, 0.88398E+05, 0.89038E+05, 0.89681E+05 &
, 0.90329E+05, 0.90979E+05, 0.91634E+05, 0.92293E+05, 0.92955E+05, 0.93621E+05, 0.94290E+05, 0.94964E+05, 0.95642E+05, 0.96323E+05 &
, 0.97008E+05, 0.97698E+05, 0.98391E+05, 0.99088E+05, 0.99789E+05, 0.10049E+06, 0.10120E+06, 0.10192E+06, 0.10263E+06, 0.10335E+06 &
, 0.10408E+06, 0.10481E+06, 0.10554E+06, 0.10628E+06, 0.10702E+06, 0.10777E+06, 0.10852E+06, 0.10927E+06, 0.11003E+06, 0.11079E+06 &
, 0.11156E+06, 0.11233E+06, 0.11311E+06, 0.11389E+06, 0.11467E+06, 0.11546E+06, 0.11625E+06, 0.11705E+06, 0.11785E+06, 0.11866E+06 &
, 0.11947E+06, 0.12028E+06, 0.12110E+06, 0.12193E+06, 0.12276E+06, 0.12359E+06, 0.12443E+06, 0.12527E+06, 0.12612E+06, 0.12697E+06 &
, 0.12783E+06, 0.12869E+06, 0.12955E+06, 0.13042E+06, 0.13130E+06, 0.13218E+06, 0.13306E+06, 0.13395E+06, 0.13485E+06, 0.13574E+06 &
, 0.13665E+06, 0.13756E+06, 0.13847E+06, 0.13939E+06, 0.14031E+06, 0.14124E+06, 0.14218E+06, 0.14311E+06, 0.14406E+06, 0.14501E+06 &
, 0.14596E+06, 0.14692E+06, 0.14788E+06, 0.14885E+06, 0.14983E+06, 0.15080E+06, 0.15179E+06, 0.15278E+06, 0.15377E+06, 0.15477E+06 &
, 0.15578E+06, 0.15679E+06, 0.15781E+06, 0.15883E+06, 0.15985E+06, 0.16089E+06, 0.16192E+06, 0.16297E+06, 0.16401E+06, 0.16507E+06 &
, 0.16613E+06, 0.16719E+06, 0.16826E+06, 0.16934E+06, 0.17042E+06, 0.17151E+06, 0.17260E+06, 0.17370E+06, 0.17480E+06, 0.17591E+06 &
, 0.17703E+06, 0.17815E+06, 0.17928E+06, 0.18041E+06, 0.18155E+06, 0.18269E+06, 0.18384E+06, 0.18500E+06, 0.18616E+06, 0.18733E+06 &
, 0.18851E+06, 0.18969E+06, 0.19087E+06, 0.19206E+06, 0.19326E+06, 0.19447E+06, 0.19568E+06, 0.19689E+06, 0.19812E+06, 0.19935E+06 &
, 0.20058E+06, 0.20182E+06, 0.20307E+06, 0.20432E+06, 0.20558E+06, 0.20685E+06, 0.20812E+06, 0.20940E+06, 0.21069E+06, 0.21198E+06 &
, 0.21328E+06, 0.21458E+06, 0.21589E+06, 0.21721E+06, 0.21854E+06, 0.21987E+06, 0.22120E+06, 0.22255E+06, 0.22390E+06, 0.22526E+06 &
, 0.22662E+06, 0.22799E+06, 0.22937E+06, 0.23075E+06, 0.23215E+06, 0.23354E+06, 0.23495E+06, 0.23636E+06, 0.23778E+06, 0.23920E+06 &
, 0.24064E+06, 0.24208E+06, 0.24352E+06, 0.24498E+06, 0.24644E+06, 0.24791E+06, 0.24938E+06, 0.25086E+06, 0.25235E+06, 0.25385E+06 &
, 0.25535E+06, 0.25686E+06, 0.25838E+06, 0.25991E+06, 0.26144E+06, 0.26298E+06, 0.26453E+06, 0.26608E+06, 0.26764E+06, 0.26922E+06 &
, 0.27079E+06, 0.27238E+06, 0.27397E+06, 0.27557E+06, 0.27718E+06, 0.27879E+06, 0.28042E+06, 0.28205E+06, 0.28369E+06, 0.28533E+06 &
, 0.28699E+06, 0.28865E+06, 0.29032E+06, 0.29199E+06, 0.29368E+06, 0.29537E+06, 0.29707E+06, 0.29878E+06, 0.30050E+06, 0.30223E+06 &
, 0.30396E+06, 0.30570E+06, 0.30745E+06, 0.30921E+06, 0.31097E+06, 0.31275E+06, 0.31453E+06, 0.31632E+06, 0.31812E+06, 0.31993E+06 &
, 0.32174E+06, 0.32357E+06, 0.32540E+06, 0.32724E+06, 0.32909E+06, 0.33095E+06, 0.33281E+06, 0.33469E+06, 0.33657E+06, 0.33847E+06 &
, 0.34037E+06, 0.34228E+06, 0.34419E+06, 0.34612E+06, 0.34806E+06, 0.35000E+06, 0.35196E+06, 0.35392E+06, 0.35589E+06, 0.35787E+06 &
, 0.35986E+06, 0.36186E+06, 0.36386E+06, 0.36588E+06, 0.36791E+06, 0.36994E+06, 0.37199E+06, 0.37404E+06, 0.37610E+06, 0.37817E+06 /
 !            5           1
 DATA(QofT(          46 ,J),J=1,501)/   7.57350015640259       &
, 0.82957E+01, 0.90181E+01, 0.97406E+01, 0.10463E+02, 0.11186E+02, 0.11909E+02, 0.12631E+02, 0.13354E+02, 0.14077E+02, 0.14800E+02 &
, 0.15523E+02, 0.16246E+02, 0.16969E+02, 0.17692E+02, 0.18415E+02, 0.19138E+02, 0.19861E+02, 0.20584E+02, 0.21307E+02, 0.22031E+02 &
, 0.22754E+02, 0.23477E+02, 0.24200E+02, 0.24923E+02, 0.25646E+02, 0.26370E+02, 0.27093E+02, 0.27816E+02, 0.28539E+02, 0.29262E+02 &
, 0.29986E+02, 0.30709E+02, 0.31432E+02, 0.32156E+02, 0.32879E+02, 0.33602E+02, 0.34325E+02, 0.35049E+02, 0.35772E+02, 0.36495E+02 &
, 0.37219E+02, 0.37942E+02, 0.38665E+02, 0.39389E+02, 0.40112E+02, 0.40836E+02, 0.41559E+02, 0.42282E+02, 0.43006E+02, 0.43729E+02 &
, 0.44453E+02, 0.45176E+02, 0.45899E+02, 0.46623E+02, 0.47346E+02, 0.48070E+02, 0.48793E+02, 0.49517E+02, 0.50240E+02, 0.50964E+02 &
, 0.51687E+02, 0.52411E+02, 0.53134E+02, 0.53858E+02, 0.54581E+02, 0.55305E+02, 0.56028E+02, 0.56752E+02, 0.57475E+02, 0.58199E+02 &
, 0.58922E+02, 0.59646E+02, 0.60370E+02, 0.61093E+02, 0.61817E+02, 0.62540E+02, 0.63264E+02, 0.63988E+02, 0.64711E+02, 0.65435E+02 &
, 0.66158E+02, 0.66882E+02, 0.67606E+02, 0.68329E+02, 0.69053E+02, 0.69777E+02, 0.70500E+02, 0.71224E+02, 0.71948E+02, 0.72671E+02 &
, 0.73395E+02, 0.74119E+02, 0.74842E+02, 0.75566E+02, 0.76290E+02, 0.77014E+02, 0.77737E+02, 0.78461E+02, 0.79185E+02, 0.79909E+02 &
, 0.80632E+02, 0.81356E+02, 0.82080E+02, 0.82804E+02, 0.83528E+02, 0.84252E+02, 0.84975E+02, 0.85699E+02, 0.86423E+02, 0.87147E+02 &
, 0.87871E+02, 0.88595E+02, 0.89318E+02, 0.90042E+02, 0.90766E+02, 0.91490E+02, 0.92214E+02, 0.92938E+02, 0.93662E+02, 0.94386E+02 &
, 0.95110E+02, 0.95834E+02, 0.96558E+02, 0.97282E+02, 0.98006E+02, 0.98730E+02, 0.99454E+02, 0.10018E+03, 0.10090E+03, 0.10163E+03 &
, 0.10235E+03, 0.10307E+03, 0.10380E+03, 0.10452E+03, 0.10525E+03, 0.10597E+03, 0.10670E+03, 0.10742E+03, 0.10814E+03, 0.10887E+03 &
, 0.10959E+03, 0.11032E+03, 0.11104E+03, 0.11177E+03, 0.11249E+03, 0.11321E+03, 0.11394E+03, 0.11466E+03, 0.11539E+03, 0.11611E+03 &
, 0.11684E+03, 0.11756E+03, 0.11829E+03, 0.11901E+03, 0.11974E+03, 0.12046E+03, 0.12119E+03, 0.12191E+03, 0.12264E+03, 0.12336E+03 &
, 0.12409E+03, 0.12481E+03, 0.12554E+03, 0.12626E+03, 0.12699E+03, 0.12771E+03, 0.12844E+03, 0.12916E+03, 0.12989E+03, 0.13061E+03 &
, 0.13134E+03, 0.13207E+03, 0.13279E+03, 0.13352E+03, 0.13424E+03, 0.13497E+03, 0.13569E+03, 0.13642E+03, 0.13715E+03, 0.13787E+03 &
, 0.13860E+03, 0.13933E+03, 0.14005E+03, 0.14078E+03, 0.14151E+03, 0.14223E+03, 0.14296E+03, 0.14369E+03, 0.14441E+03, 0.14514E+03 &
, 0.14587E+03, 0.14660E+03, 0.14732E+03, 0.14805E+03, 0.14878E+03, 0.14951E+03, 0.15023E+03, 0.15096E+03, 0.15169E+03, 0.15242E+03 &
, 0.15315E+03, 0.15388E+03, 0.15460E+03, 0.15533E+03, 0.15606E+03, 0.15679E+03, 0.15752E+03, 0.15825E+03, 0.15898E+03, 0.15971E+03 &
, 0.16044E+03, 0.16117E+03, 0.16190E+03, 0.16263E+03, 0.16336E+03, 0.16409E+03, 0.16482E+03, 0.16555E+03, 0.16628E+03, 0.16702E+03 &
, 0.16775E+03, 0.16848E+03, 0.16921E+03, 0.16994E+03, 0.17068E+03, 0.17141E+03, 0.17214E+03, 0.17287E+03, 0.17361E+03, 0.17434E+03 &
, 0.17507E+03, 0.17581E+03, 0.17654E+03, 0.17728E+03, 0.17801E+03, 0.17875E+03, 0.17948E+03, 0.18022E+03, 0.18095E+03, 0.18169E+03 &
, 0.18242E+03, 0.18316E+03, 0.18389E+03, 0.18463E+03, 0.18537E+03, 0.18611E+03, 0.18684E+03, 0.18758E+03, 0.18832E+03, 0.18906E+03 &
, 0.18979E+03, 0.19053E+03, 0.19127E+03, 0.19201E+03, 0.19275E+03, 0.19349E+03, 0.19423E+03, 0.19497E+03, 0.19571E+03, 0.19645E+03 &
, 0.19720E+03, 0.19794E+03, 0.19868E+03, 0.19942E+03, 0.20016E+03, 0.20091E+03, 0.20165E+03, 0.20239E+03, 0.20314E+03, 0.20388E+03 &
, 0.20463E+03, 0.20537E+03, 0.20612E+03, 0.20686E+03, 0.20761E+03, 0.20835E+03, 0.20910E+03, 0.20985E+03, 0.21060E+03, 0.21134E+03 &
, 0.21209E+03, 0.21284E+03, 0.21359E+03, 0.21434E+03, 0.21509E+03, 0.21584E+03, 0.21659E+03, 0.21734E+03, 0.21809E+03, 0.21884E+03 &
, 0.21959E+03, 0.22034E+03, 0.22110E+03, 0.22185E+03, 0.22260E+03, 0.22336E+03, 0.22411E+03, 0.22487E+03, 0.22562E+03, 0.22638E+03 &
, 0.22713E+03, 0.22789E+03, 0.22865E+03, 0.22940E+03, 0.23016E+03, 0.23092E+03, 0.23168E+03, 0.23244E+03, 0.23320E+03, 0.23396E+03 &
, 0.23472E+03, 0.23548E+03, 0.23624E+03, 0.23700E+03, 0.23776E+03, 0.23853E+03, 0.23929E+03, 0.24005E+03, 0.24082E+03, 0.24158E+03 &
, 0.24235E+03, 0.24311E+03, 0.24388E+03, 0.24464E+03, 0.24541E+03, 0.24618E+03, 0.24695E+03, 0.24771E+03, 0.24848E+03, 0.24925E+03 &
, 0.25002E+03, 0.25079E+03, 0.25156E+03, 0.25234E+03, 0.25311E+03, 0.25388E+03, 0.25465E+03, 0.25543E+03, 0.25620E+03, 0.25698E+03 &
, 0.25775E+03, 0.25853E+03, 0.25930E+03, 0.26008E+03, 0.26086E+03, 0.26163E+03, 0.26241E+03, 0.26319E+03, 0.26397E+03, 0.26475E+03 &
, 0.26553E+03, 0.26631E+03, 0.26709E+03, 0.26788E+03, 0.26866E+03, 0.26944E+03, 0.27023E+03, 0.27101E+03, 0.27180E+03, 0.27258E+03 &
, 0.27337E+03, 0.27415E+03, 0.27494E+03, 0.27573E+03, 0.27652E+03, 0.27731E+03, 0.27810E+03, 0.27889E+03, 0.27968E+03, 0.28047E+03 &
, 0.28126E+03, 0.28205E+03, 0.28285E+03, 0.28364E+03, 0.28444E+03, 0.28523E+03, 0.28603E+03, 0.28682E+03, 0.28762E+03, 0.28842E+03 &
, 0.28922E+03, 0.29001E+03, 0.29081E+03, 0.29161E+03, 0.29241E+03, 0.29322E+03, 0.29402E+03, 0.29482E+03, 0.29562E+03, 0.29643E+03 &
, 0.29723E+03, 0.29804E+03, 0.29884E+03, 0.29965E+03, 0.30046E+03, 0.30126E+03, 0.30207E+03, 0.30288E+03, 0.30369E+03, 0.30450E+03 &
, 0.30531E+03, 0.30612E+03, 0.30693E+03, 0.30775E+03, 0.30856E+03, 0.30938E+03, 0.31019E+03, 0.31101E+03, 0.31182E+03, 0.31264E+03 &
, 0.31346E+03, 0.31427E+03, 0.31509E+03, 0.31591E+03, 0.31673E+03, 0.31755E+03, 0.31838E+03, 0.31920E+03, 0.32002E+03, 0.32085E+03 &
, 0.32167E+03, 0.32249E+03, 0.32332E+03, 0.32415E+03, 0.32497E+03, 0.32580E+03, 0.32663E+03, 0.32746E+03, 0.32829E+03, 0.32912E+03 &
, 0.32995E+03, 0.33078E+03, 0.33162E+03, 0.33245E+03, 0.33328E+03, 0.33412E+03, 0.33495E+03, 0.33579E+03, 0.33663E+03, 0.33747E+03 &
, 0.33830E+03, 0.33914E+03, 0.33998E+03, 0.34082E+03, 0.34166E+03, 0.34251E+03, 0.34335E+03, 0.34419E+03, 0.34504E+03, 0.34588E+03 &
, 0.34673E+03, 0.34757E+03, 0.34842E+03, 0.34927E+03, 0.35012E+03, 0.35097E+03, 0.35182E+03, 0.35267E+03, 0.35352E+03, 0.35437E+03 &
, 0.35523E+03, 0.35608E+03, 0.35693E+03, 0.35779E+03, 0.35865E+03, 0.35950E+03, 0.36036E+03, 0.36122E+03, 0.36208E+03, 0.36294E+03 &
, 0.36380E+03, 0.36466E+03, 0.36552E+03, 0.36639E+03, 0.36725E+03, 0.36811E+03, 0.36898E+03, 0.36984E+03, 0.37071E+03, 0.37158E+03 &
, 0.37245E+03, 0.37332E+03, 0.37419E+03, 0.37506E+03, 0.37593E+03, 0.37680E+03, 0.37767E+03, 0.37855E+03, 0.37942E+03, 0.38030E+03 &
, 0.38117E+03, 0.38205E+03, 0.38293E+03, 0.38381E+03, 0.38469E+03, 0.38557E+03, 0.38645E+03, 0.38733E+03, 0.38821E+03, 0.38909E+03 /
 !            5           2
 DATA(QofT(          47 ,J),J=1,501)/   15.8109998703003       &
, 0.17323E+02, 0.18834E+02, 0.20346E+02, 0.21857E+02, 0.23369E+02, 0.24881E+02, 0.26393E+02, 0.27906E+02, 0.29418E+02, 0.30930E+02 &
, 0.32443E+02, 0.33955E+02, 0.35468E+02, 0.36981E+02, 0.38493E+02, 0.40006E+02, 0.41519E+02, 0.43031E+02, 0.44544E+02, 0.46057E+02 &
, 0.47570E+02, 0.49083E+02, 0.50596E+02, 0.52108E+02, 0.53621E+02, 0.55134E+02, 0.56647E+02, 0.58160E+02, 0.59673E+02, 0.61186E+02 &
, 0.62699E+02, 0.64213E+02, 0.65726E+02, 0.67239E+02, 0.68752E+02, 0.70265E+02, 0.71778E+02, 0.73291E+02, 0.74805E+02, 0.76318E+02 &
, 0.77831E+02, 0.79344E+02, 0.80858E+02, 0.82371E+02, 0.83884E+02, 0.85397E+02, 0.86911E+02, 0.88424E+02, 0.89938E+02, 0.91451E+02 &
, 0.92964E+02, 0.94478E+02, 0.95991E+02, 0.97505E+02, 0.99018E+02, 0.10053E+03, 0.10204E+03, 0.10356E+03, 0.10507E+03, 0.10659E+03 &
, 0.10810E+03, 0.10961E+03, 0.11113E+03, 0.11264E+03, 0.11415E+03, 0.11567E+03, 0.11718E+03, 0.11869E+03, 0.12021E+03, 0.12172E+03 &
, 0.12324E+03, 0.12475E+03, 0.12626E+03, 0.12778E+03, 0.12929E+03, 0.13080E+03, 0.13232E+03, 0.13383E+03, 0.13535E+03, 0.13686E+03 &
, 0.13837E+03, 0.13989E+03, 0.14140E+03, 0.14291E+03, 0.14443E+03, 0.14594E+03, 0.14746E+03, 0.14897E+03, 0.15048E+03, 0.15200E+03 &
, 0.15351E+03, 0.15503E+03, 0.15654E+03, 0.15805E+03, 0.15957E+03, 0.16108E+03, 0.16260E+03, 0.16411E+03, 0.16563E+03, 0.16714E+03 &
, 0.16865E+03, 0.17017E+03, 0.17168E+03, 0.17320E+03, 0.17471E+03, 0.17622E+03, 0.17774E+03, 0.17925E+03, 0.18077E+03, 0.18228E+03 &
, 0.18380E+03, 0.18531E+03, 0.18682E+03, 0.18834E+03, 0.18985E+03, 0.19137E+03, 0.19288E+03, 0.19440E+03, 0.19591E+03, 0.19743E+03 &
, 0.19894E+03, 0.20046E+03, 0.20197E+03, 0.20348E+03, 0.20500E+03, 0.20651E+03, 0.20803E+03, 0.20954E+03, 0.21106E+03, 0.21257E+03 &
, 0.21409E+03, 0.21560E+03, 0.21712E+03, 0.21863E+03, 0.22015E+03, 0.22166E+03, 0.22318E+03, 0.22469E+03, 0.22621E+03, 0.22773E+03 &
, 0.22924E+03, 0.23076E+03, 0.23227E+03, 0.23379E+03, 0.23530E+03, 0.23682E+03, 0.23833E+03, 0.23985E+03, 0.24137E+03, 0.24288E+03 &
, 0.24440E+03, 0.24592E+03, 0.24743E+03, 0.24895E+03, 0.25046E+03, 0.25198E+03, 0.25350E+03, 0.25501E+03, 0.25653E+03, 0.25805E+03 &
, 0.25957E+03, 0.26108E+03, 0.26260E+03, 0.26412E+03, 0.26564E+03, 0.26715E+03, 0.26867E+03, 0.27019E+03, 0.27171E+03, 0.27323E+03 &
, 0.27474E+03, 0.27626E+03, 0.27778E+03, 0.27930E+03, 0.28082E+03, 0.28234E+03, 0.28386E+03, 0.28538E+03, 0.28690E+03, 0.28842E+03 &
, 0.28994E+03, 0.29146E+03, 0.29298E+03, 0.29450E+03, 0.29602E+03, 0.29754E+03, 0.29906E+03, 0.30059E+03, 0.30211E+03, 0.30363E+03 &
, 0.30515E+03, 0.30668E+03, 0.30820E+03, 0.30972E+03, 0.31125E+03, 0.31277E+03, 0.31429E+03, 0.31582E+03, 0.31734E+03, 0.31887E+03 &
, 0.32039E+03, 0.32192E+03, 0.32344E+03, 0.32497E+03, 0.32650E+03, 0.32802E+03, 0.32955E+03, 0.33108E+03, 0.33261E+03, 0.33413E+03 &
, 0.33566E+03, 0.33719E+03, 0.33872E+03, 0.34025E+03, 0.34178E+03, 0.34331E+03, 0.34484E+03, 0.34637E+03, 0.34790E+03, 0.34944E+03 &
, 0.35097E+03, 0.35250E+03, 0.35404E+03, 0.35557E+03, 0.35710E+03, 0.35864E+03, 0.36017E+03, 0.36171E+03, 0.36325E+03, 0.36478E+03 &
, 0.36632E+03, 0.36786E+03, 0.36939E+03, 0.37093E+03, 0.37247E+03, 0.37401E+03, 0.37555E+03, 0.37709E+03, 0.37863E+03, 0.38018E+03 &
, 0.38172E+03, 0.38326E+03, 0.38480E+03, 0.38635E+03, 0.38789E+03, 0.38944E+03, 0.39098E+03, 0.39253E+03, 0.39408E+03, 0.39563E+03 &
, 0.39717E+03, 0.39872E+03, 0.40027E+03, 0.40182E+03, 0.40337E+03, 0.40492E+03, 0.40648E+03, 0.40803E+03, 0.40958E+03, 0.41114E+03 &
, 0.41269E+03, 0.41425E+03, 0.41580E+03, 0.41736E+03, 0.41892E+03, 0.42048E+03, 0.42204E+03, 0.42360E+03, 0.42516E+03, 0.42672E+03 &
, 0.42828E+03, 0.42984E+03, 0.43141E+03, 0.43297E+03, 0.43454E+03, 0.43610E+03, 0.43767E+03, 0.43924E+03, 0.44080E+03, 0.44237E+03 &
, 0.44394E+03, 0.44551E+03, 0.44709E+03, 0.44866E+03, 0.45023E+03, 0.45181E+03, 0.45338E+03, 0.45496E+03, 0.45653E+03, 0.45811E+03 &
, 0.45969E+03, 0.46127E+03, 0.46285E+03, 0.46443E+03, 0.46601E+03, 0.46760E+03, 0.46918E+03, 0.47077E+03, 0.47235E+03, 0.47394E+03 &
, 0.47553E+03, 0.47711E+03, 0.47870E+03, 0.48030E+03, 0.48189E+03, 0.48348E+03, 0.48507E+03, 0.48667E+03, 0.48826E+03, 0.48986E+03 &
, 0.49146E+03, 0.49306E+03, 0.49465E+03, 0.49626E+03, 0.49786E+03, 0.49946E+03, 0.50106E+03, 0.50267E+03, 0.50427E+03, 0.50588E+03 &
, 0.50749E+03, 0.50910E+03, 0.51071E+03, 0.51232E+03, 0.51393E+03, 0.51554E+03, 0.51716E+03, 0.51877E+03, 0.52039E+03, 0.52201E+03 &
, 0.52363E+03, 0.52525E+03, 0.52687E+03, 0.52849E+03, 0.53011E+03, 0.53174E+03, 0.53336E+03, 0.53499E+03, 0.53662E+03, 0.53825E+03 &
, 0.53988E+03, 0.54151E+03, 0.54314E+03, 0.54478E+03, 0.54641E+03, 0.54805E+03, 0.54969E+03, 0.55132E+03, 0.55296E+03, 0.55461E+03 &
, 0.55625E+03, 0.55789E+03, 0.55954E+03, 0.56118E+03, 0.56283E+03, 0.56448E+03, 0.56613E+03, 0.56778E+03, 0.56943E+03, 0.57108E+03 &
, 0.57274E+03, 0.57439E+03, 0.57605E+03, 0.57771E+03, 0.57937E+03, 0.58103E+03, 0.58269E+03, 0.58436E+03, 0.58602E+03, 0.58769E+03 &
, 0.58936E+03, 0.59103E+03, 0.59270E+03, 0.59437E+03, 0.59604E+03, 0.59772E+03, 0.59939E+03, 0.60107E+03, 0.60275E+03, 0.60443E+03 &
, 0.60611E+03, 0.60779E+03, 0.60947E+03, 0.61116E+03, 0.61285E+03, 0.61453E+03, 0.61622E+03, 0.61791E+03, 0.61961E+03, 0.62130E+03 &
, 0.62299E+03, 0.62469E+03, 0.62639E+03, 0.62809E+03, 0.62979E+03, 0.63149E+03, 0.63319E+03, 0.63490E+03, 0.63660E+03, 0.63831E+03 &
, 0.64002E+03, 0.64173E+03, 0.64344E+03, 0.64516E+03, 0.64687E+03, 0.64859E+03, 0.65030E+03, 0.65202E+03, 0.65374E+03, 0.65547E+03 &
, 0.65719E+03, 0.65892E+03, 0.66064E+03, 0.66237E+03, 0.66410E+03, 0.66583E+03, 0.66756E+03, 0.66930E+03, 0.67103E+03, 0.67277E+03 &
, 0.67451E+03, 0.67625E+03, 0.67799E+03, 0.67973E+03, 0.68148E+03, 0.68322E+03, 0.68497E+03, 0.68672E+03, 0.68847E+03, 0.69022E+03 &
, 0.69198E+03, 0.69373E+03, 0.69549E+03, 0.69725E+03, 0.69901E+03, 0.70077E+03, 0.70253E+03, 0.70430E+03, 0.70606E+03, 0.70783E+03 &
, 0.70960E+03, 0.71137E+03, 0.71314E+03, 0.71492E+03, 0.71669E+03, 0.71847E+03, 0.72025E+03, 0.72203E+03, 0.72381E+03, 0.72559E+03 &
, 0.72738E+03, 0.72917E+03, 0.73095E+03, 0.73274E+03, 0.73453E+03, 0.73633E+03, 0.73812E+03, 0.73992E+03, 0.74172E+03, 0.74352E+03 &
, 0.74532E+03, 0.74712E+03, 0.74893E+03, 0.75073E+03, 0.75254E+03, 0.75435E+03, 0.75616E+03, 0.75797E+03, 0.75979E+03, 0.76160E+03 &
, 0.76342E+03, 0.76524E+03, 0.76706E+03, 0.76888E+03, 0.77071E+03, 0.77253E+03, 0.77436E+03, 0.77619E+03, 0.77802E+03, 0.77985E+03 &
, 0.78169E+03, 0.78352E+03, 0.78536E+03, 0.78720E+03, 0.78904E+03, 0.79088E+03, 0.79273E+03, 0.79457E+03, 0.79642E+03, 0.79827E+03 &
, 0.80012E+03, 0.80197E+03, 0.80383E+03, 0.80568E+03, 0.80754E+03, 0.80940E+03, 0.81126E+03, 0.81313E+03, 0.81499E+03, 0.81686E+03 /
 !            5           3
 DATA(QofT(          48 ,J),J=1,501)/   7.93459987640381       &
, 0.86930E+01, 0.94516E+01, 0.10210E+02, 0.10969E+02, 0.11728E+02, 0.12487E+02, 0.13246E+02, 0.14005E+02, 0.14764E+02, 0.15523E+02 &
, 0.16282E+02, 0.17041E+02, 0.17800E+02, 0.18560E+02, 0.19319E+02, 0.20078E+02, 0.20837E+02, 0.21597E+02, 0.22356E+02, 0.23115E+02 &
, 0.23874E+02, 0.24634E+02, 0.25393E+02, 0.26152E+02, 0.26912E+02, 0.27671E+02, 0.28431E+02, 0.29190E+02, 0.29949E+02, 0.30709E+02 &
, 0.31468E+02, 0.32228E+02, 0.32987E+02, 0.33746E+02, 0.34506E+02, 0.35265E+02, 0.36025E+02, 0.36784E+02, 0.37544E+02, 0.38303E+02 &
, 0.39063E+02, 0.39822E+02, 0.40582E+02, 0.41341E+02, 0.42101E+02, 0.42860E+02, 0.43620E+02, 0.44380E+02, 0.45139E+02, 0.45899E+02 &
, 0.46658E+02, 0.47418E+02, 0.48178E+02, 0.48937E+02, 0.49697E+02, 0.50456E+02, 0.51216E+02, 0.51976E+02, 0.52735E+02, 0.53495E+02 &
, 0.54255E+02, 0.55014E+02, 0.55774E+02, 0.56534E+02, 0.57293E+02, 0.58053E+02, 0.58813E+02, 0.59573E+02, 0.60332E+02, 0.61092E+02 &
, 0.61852E+02, 0.62611E+02, 0.63371E+02, 0.64131E+02, 0.64891E+02, 0.65651E+02, 0.66410E+02, 0.67170E+02, 0.67930E+02, 0.68690E+02 &
, 0.69450E+02, 0.70209E+02, 0.70969E+02, 0.71729E+02, 0.72489E+02, 0.73249E+02, 0.74009E+02, 0.74768E+02, 0.75528E+02, 0.76288E+02 &
, 0.77048E+02, 0.77808E+02, 0.78568E+02, 0.79328E+02, 0.80088E+02, 0.80848E+02, 0.81608E+02, 0.82367E+02, 0.83127E+02, 0.83887E+02 &
, 0.84647E+02, 0.85407E+02, 0.86167E+02, 0.86927E+02, 0.87687E+02, 0.88447E+02, 0.89207E+02, 0.89967E+02, 0.90727E+02, 0.91488E+02 &
, 0.92248E+02, 0.93008E+02, 0.93768E+02, 0.94528E+02, 0.95288E+02, 0.96048E+02, 0.96808E+02, 0.97568E+02, 0.98329E+02, 0.99089E+02 &
, 0.99849E+02, 0.10061E+03, 0.10137E+03, 0.10213E+03, 0.10289E+03, 0.10365E+03, 0.10441E+03, 0.10517E+03, 0.10593E+03, 0.10669E+03 &
, 0.10745E+03, 0.10821E+03, 0.10897E+03, 0.10973E+03, 0.11049E+03, 0.11125E+03, 0.11201E+03, 0.11277E+03, 0.11354E+03, 0.11430E+03 &
, 0.11506E+03, 0.11582E+03, 0.11658E+03, 0.11734E+03, 0.11810E+03, 0.11886E+03, 0.11962E+03, 0.12038E+03, 0.12114E+03, 0.12190E+03 &
, 0.12266E+03, 0.12343E+03, 0.12419E+03, 0.12495E+03, 0.12571E+03, 0.12647E+03, 0.12723E+03, 0.12799E+03, 0.12875E+03, 0.12952E+03 &
, 0.13028E+03, 0.13104E+03, 0.13180E+03, 0.13256E+03, 0.13332E+03, 0.13409E+03, 0.13485E+03, 0.13561E+03, 0.13637E+03, 0.13713E+03 &
, 0.13790E+03, 0.13866E+03, 0.13942E+03, 0.14018E+03, 0.14095E+03, 0.14171E+03, 0.14247E+03, 0.14323E+03, 0.14400E+03, 0.14476E+03 &
, 0.14552E+03, 0.14629E+03, 0.14705E+03, 0.14781E+03, 0.14858E+03, 0.14934E+03, 0.15010E+03, 0.15087E+03, 0.15163E+03, 0.15239E+03 &
, 0.15316E+03, 0.15392E+03, 0.15469E+03, 0.15545E+03, 0.15622E+03, 0.15698E+03, 0.15775E+03, 0.15851E+03, 0.15928E+03, 0.16004E+03 &
, 0.16081E+03, 0.16157E+03, 0.16234E+03, 0.16311E+03, 0.16387E+03, 0.16464E+03, 0.16541E+03, 0.16617E+03, 0.16694E+03, 0.16771E+03 &
, 0.16847E+03, 0.16924E+03, 0.17001E+03, 0.17078E+03, 0.17154E+03, 0.17231E+03, 0.17308E+03, 0.17385E+03, 0.17462E+03, 0.17539E+03 &
, 0.17616E+03, 0.17693E+03, 0.17770E+03, 0.17847E+03, 0.17924E+03, 0.18001E+03, 0.18078E+03, 0.18155E+03, 0.18232E+03, 0.18309E+03 &
, 0.18386E+03, 0.18463E+03, 0.18541E+03, 0.18618E+03, 0.18695E+03, 0.18772E+03, 0.18850E+03, 0.18927E+03, 0.19004E+03, 0.19082E+03 &
, 0.19159E+03, 0.19237E+03, 0.19314E+03, 0.19392E+03, 0.19469E+03, 0.19547E+03, 0.19624E+03, 0.19702E+03, 0.19780E+03, 0.19857E+03 &
, 0.19935E+03, 0.20013E+03, 0.20091E+03, 0.20168E+03, 0.20246E+03, 0.20324E+03, 0.20402E+03, 0.20480E+03, 0.20558E+03, 0.20636E+03 &
, 0.20714E+03, 0.20792E+03, 0.20870E+03, 0.20948E+03, 0.21027E+03, 0.21105E+03, 0.21183E+03, 0.21262E+03, 0.21340E+03, 0.21418E+03 &
, 0.21497E+03, 0.21575E+03, 0.21654E+03, 0.21732E+03, 0.21811E+03, 0.21889E+03, 0.21968E+03, 0.22047E+03, 0.22125E+03, 0.22204E+03 &
, 0.22283E+03, 0.22362E+03, 0.22441E+03, 0.22520E+03, 0.22599E+03, 0.22678E+03, 0.22757E+03, 0.22836E+03, 0.22915E+03, 0.22994E+03 &
, 0.23074E+03, 0.23153E+03, 0.23232E+03, 0.23312E+03, 0.23391E+03, 0.23471E+03, 0.23550E+03, 0.23630E+03, 0.23709E+03, 0.23789E+03 &
, 0.23869E+03, 0.23948E+03, 0.24028E+03, 0.24108E+03, 0.24188E+03, 0.24268E+03, 0.24348E+03, 0.24428E+03, 0.24508E+03, 0.24588E+03 &
, 0.24669E+03, 0.24749E+03, 0.24829E+03, 0.24910E+03, 0.24990E+03, 0.25070E+03, 0.25151E+03, 0.25232E+03, 0.25312E+03, 0.25393E+03 &
, 0.25474E+03, 0.25554E+03, 0.25635E+03, 0.25716E+03, 0.25797E+03, 0.25878E+03, 0.25959E+03, 0.26040E+03, 0.26121E+03, 0.26203E+03 &
, 0.26284E+03, 0.26365E+03, 0.26447E+03, 0.26528E+03, 0.26610E+03, 0.26691E+03, 0.26773E+03, 0.26855E+03, 0.26936E+03, 0.27018E+03 &
, 0.27100E+03, 0.27182E+03, 0.27264E+03, 0.27346E+03, 0.27428E+03, 0.27510E+03, 0.27592E+03, 0.27675E+03, 0.27757E+03, 0.27839E+03 &
, 0.27922E+03, 0.28004E+03, 0.28087E+03, 0.28170E+03, 0.28252E+03, 0.28335E+03, 0.28418E+03, 0.28501E+03, 0.28584E+03, 0.28667E+03 &
, 0.28750E+03, 0.28833E+03, 0.28916E+03, 0.29000E+03, 0.29083E+03, 0.29166E+03, 0.29250E+03, 0.29334E+03, 0.29417E+03, 0.29501E+03 &
, 0.29585E+03, 0.29668E+03, 0.29752E+03, 0.29836E+03, 0.29920E+03, 0.30004E+03, 0.30088E+03, 0.30173E+03, 0.30257E+03, 0.30341E+03 &
, 0.30426E+03, 0.30510E+03, 0.30595E+03, 0.30679E+03, 0.30764E+03, 0.30849E+03, 0.30934E+03, 0.31019E+03, 0.31104E+03, 0.31189E+03 &
, 0.31274E+03, 0.31359E+03, 0.31444E+03, 0.31530E+03, 0.31615E+03, 0.31700E+03, 0.31786E+03, 0.31872E+03, 0.31957E+03, 0.32043E+03 &
, 0.32129E+03, 0.32215E+03, 0.32301E+03, 0.32387E+03, 0.32473E+03, 0.32559E+03, 0.32645E+03, 0.32732E+03, 0.32818E+03, 0.32905E+03 &
, 0.32991E+03, 0.33078E+03, 0.33165E+03, 0.33251E+03, 0.33338E+03, 0.33425E+03, 0.33512E+03, 0.33599E+03, 0.33687E+03, 0.33774E+03 &
, 0.33861E+03, 0.33948E+03, 0.34036E+03, 0.34123E+03, 0.34211E+03, 0.34299E+03, 0.34387E+03, 0.34474E+03, 0.34562E+03, 0.34650E+03 &
, 0.34738E+03, 0.34827E+03, 0.34915E+03, 0.35003E+03, 0.35092E+03, 0.35180E+03, 0.35269E+03, 0.35357E+03, 0.35446E+03, 0.35535E+03 &
, 0.35624E+03, 0.35713E+03, 0.35802E+03, 0.35891E+03, 0.35980E+03, 0.36069E+03, 0.36158E+03, 0.36248E+03, 0.36337E+03, 0.36427E+03 &
, 0.36517E+03, 0.36606E+03, 0.36696E+03, 0.36786E+03, 0.36876E+03, 0.36966E+03, 0.37056E+03, 0.37147E+03, 0.37237E+03, 0.37327E+03 &
, 0.37418E+03, 0.37508E+03, 0.37599E+03, 0.37690E+03, 0.37780E+03, 0.37871E+03, 0.37962E+03, 0.38053E+03, 0.38144E+03, 0.38236E+03 &
, 0.38327E+03, 0.38418E+03, 0.38510E+03, 0.38601E+03, 0.38693E+03, 0.38785E+03, 0.38877E+03, 0.38968E+03, 0.39060E+03, 0.39152E+03 &
, 0.39245E+03, 0.39337E+03, 0.39429E+03, 0.39522E+03, 0.39614E+03, 0.39707E+03, 0.39799E+03, 0.39892E+03, 0.39985E+03, 0.40078E+03 &
, 0.40171E+03, 0.40264E+03, 0.40357E+03, 0.40450E+03, 0.40543E+03, 0.40637E+03, 0.40730E+03, 0.40824E+03, 0.40918E+03, 0.41011E+03 /
 !            5           4
 DATA(QofT(          49 ,J),J=1,501)/   46.5639991760254       &
, 0.51010E+02, 0.55457E+02, 0.59904E+02, 0.64352E+02, 0.68801E+02, 0.73250E+02, 0.77699E+02, 0.82149E+02, 0.86599E+02, 0.91049E+02 &
, 0.95499E+02, 0.99949E+02, 0.10440E+03, 0.10885E+03, 0.11330E+03, 0.11775E+03, 0.12220E+03, 0.12665E+03, 0.13111E+03, 0.13556E+03 &
, 0.14001E+03, 0.14446E+03, 0.14891E+03, 0.15336E+03, 0.15781E+03, 0.16227E+03, 0.16672E+03, 0.17117E+03, 0.17562E+03, 0.18007E+03 &
, 0.18452E+03, 0.18898E+03, 0.19343E+03, 0.19788E+03, 0.20233E+03, 0.20679E+03, 0.21124E+03, 0.21569E+03, 0.22014E+03, 0.22459E+03 &
, 0.22905E+03, 0.23350E+03, 0.23795E+03, 0.24241E+03, 0.24686E+03, 0.25131E+03, 0.25576E+03, 0.26022E+03, 0.26467E+03, 0.26912E+03 &
, 0.27358E+03, 0.27803E+03, 0.28248E+03, 0.28693E+03, 0.29139E+03, 0.29584E+03, 0.30029E+03, 0.30475E+03, 0.30920E+03, 0.31365E+03 &
, 0.31811E+03, 0.32256E+03, 0.32701E+03, 0.33147E+03, 0.33592E+03, 0.34038E+03, 0.34483E+03, 0.34928E+03, 0.35374E+03, 0.35819E+03 &
, 0.36264E+03, 0.36710E+03, 0.37155E+03, 0.37601E+03, 0.38046E+03, 0.38491E+03, 0.38937E+03, 0.39382E+03, 0.39828E+03, 0.40273E+03 &
, 0.40719E+03, 0.41164E+03, 0.41609E+03, 0.42055E+03, 0.42500E+03, 0.42946E+03, 0.43391E+03, 0.43837E+03, 0.44282E+03, 0.44728E+03 &
, 0.45173E+03, 0.45619E+03, 0.46064E+03, 0.46510E+03, 0.46955E+03, 0.47401E+03, 0.47846E+03, 0.48292E+03, 0.48737E+03, 0.49183E+03 &
, 0.49628E+03, 0.50074E+03, 0.50519E+03, 0.50965E+03, 0.51410E+03, 0.51856E+03, 0.52301E+03, 0.52747E+03, 0.53192E+03, 0.53638E+03 &
, 0.54084E+03, 0.54529E+03, 0.54975E+03, 0.55420E+03, 0.55866E+03, 0.56312E+03, 0.56757E+03, 0.57203E+03, 0.57648E+03, 0.58094E+03 &
, 0.58540E+03, 0.58985E+03, 0.59431E+03, 0.59877E+03, 0.60322E+03, 0.60768E+03, 0.61214E+03, 0.61660E+03, 0.62105E+03, 0.62551E+03 &
, 0.62997E+03, 0.63443E+03, 0.63888E+03, 0.64334E+03, 0.64780E+03, 0.65226E+03, 0.65671E+03, 0.66117E+03, 0.66563E+03, 0.67009E+03 &
, 0.67455E+03, 0.67901E+03, 0.68347E+03, 0.68793E+03, 0.69239E+03, 0.69685E+03, 0.70131E+03, 0.70577E+03, 0.71023E+03, 0.71469E+03 &
, 0.71915E+03, 0.72361E+03, 0.72807E+03, 0.73253E+03, 0.73699E+03, 0.74146E+03, 0.74592E+03, 0.75038E+03, 0.75484E+03, 0.75931E+03 &
, 0.76377E+03, 0.76823E+03, 0.77270E+03, 0.77716E+03, 0.78163E+03, 0.78609E+03, 0.79056E+03, 0.79503E+03, 0.79949E+03, 0.80396E+03 &
, 0.80843E+03, 0.81289E+03, 0.81736E+03, 0.82183E+03, 0.82630E+03, 0.83077E+03, 0.83524E+03, 0.83971E+03, 0.84418E+03, 0.84865E+03 &
, 0.85312E+03, 0.85760E+03, 0.86207E+03, 0.86654E+03, 0.87102E+03, 0.87549E+03, 0.87997E+03, 0.88445E+03, 0.88892E+03, 0.89340E+03 &
, 0.89788E+03, 0.90236E+03, 0.90684E+03, 0.91132E+03, 0.91580E+03, 0.92028E+03, 0.92477E+03, 0.92925E+03, 0.93373E+03, 0.93822E+03 &
, 0.94271E+03, 0.94719E+03, 0.95168E+03, 0.95617E+03, 0.96066E+03, 0.96515E+03, 0.96964E+03, 0.97413E+03, 0.97863E+03, 0.98312E+03 &
, 0.98762E+03, 0.99211E+03, 0.99661E+03, 0.10011E+04, 0.10056E+04, 0.10101E+04, 0.10146E+04, 0.10191E+04, 0.10236E+04, 0.10281E+04 &
, 0.10326E+04, 0.10371E+04, 0.10416E+04, 0.10462E+04, 0.10507E+04, 0.10552E+04, 0.10597E+04, 0.10642E+04, 0.10687E+04, 0.10732E+04 &
, 0.10778E+04, 0.10823E+04, 0.10868E+04, 0.10913E+04, 0.10959E+04, 0.11004E+04, 0.11049E+04, 0.11094E+04, 0.11140E+04, 0.11185E+04 &
, 0.11230E+04, 0.11276E+04, 0.11321E+04, 0.11367E+04, 0.11412E+04, 0.11457E+04, 0.11503E+04, 0.11548E+04, 0.11594E+04, 0.11639E+04 &
, 0.11685E+04, 0.11730E+04, 0.11776E+04, 0.11821E+04, 0.11867E+04, 0.11913E+04, 0.11958E+04, 0.12004E+04, 0.12050E+04, 0.12095E+04 &
, 0.12141E+04, 0.12187E+04, 0.12232E+04, 0.12278E+04, 0.12324E+04, 0.12370E+04, 0.12416E+04, 0.12461E+04, 0.12507E+04, 0.12553E+04 &
, 0.12599E+04, 0.12645E+04, 0.12691E+04, 0.12737E+04, 0.12783E+04, 0.12829E+04, 0.12875E+04, 0.12921E+04, 0.12967E+04, 0.13013E+04 &
, 0.13059E+04, 0.13105E+04, 0.13152E+04, 0.13198E+04, 0.13244E+04, 0.13290E+04, 0.13337E+04, 0.13383E+04, 0.13429E+04, 0.13476E+04 &
, 0.13522E+04, 0.13568E+04, 0.13615E+04, 0.13661E+04, 0.13708E+04, 0.13754E+04, 0.13801E+04, 0.13847E+04, 0.13894E+04, 0.13940E+04 &
, 0.13987E+04, 0.14034E+04, 0.14080E+04, 0.14127E+04, 0.14174E+04, 0.14221E+04, 0.14267E+04, 0.14314E+04, 0.14361E+04, 0.14408E+04 &
, 0.14455E+04, 0.14502E+04, 0.14549E+04, 0.14596E+04, 0.14643E+04, 0.14690E+04, 0.14737E+04, 0.14784E+04, 0.14831E+04, 0.14879E+04 &
, 0.14926E+04, 0.14973E+04, 0.15020E+04, 0.15068E+04, 0.15115E+04, 0.15162E+04, 0.15210E+04, 0.15257E+04, 0.15305E+04, 0.15352E+04 &
, 0.15400E+04, 0.15447E+04, 0.15495E+04, 0.15542E+04, 0.15590E+04, 0.15638E+04, 0.15685E+04, 0.15733E+04, 0.15781E+04, 0.15829E+04 &
, 0.15877E+04, 0.15925E+04, 0.15973E+04, 0.16020E+04, 0.16068E+04, 0.16117E+04, 0.16165E+04, 0.16213E+04, 0.16261E+04, 0.16309E+04 &
, 0.16357E+04, 0.16405E+04, 0.16454E+04, 0.16502E+04, 0.16550E+04, 0.16599E+04, 0.16647E+04, 0.16696E+04, 0.16744E+04, 0.16793E+04 &
, 0.16841E+04, 0.16890E+04, 0.16938E+04, 0.16987E+04, 0.17036E+04, 0.17084E+04, 0.17133E+04, 0.17182E+04, 0.17231E+04, 0.17280E+04 &
, 0.17329E+04, 0.17378E+04, 0.17427E+04, 0.17476E+04, 0.17525E+04, 0.17574E+04, 0.17623E+04, 0.17672E+04, 0.17722E+04, 0.17771E+04 &
, 0.17820E+04, 0.17870E+04, 0.17919E+04, 0.17968E+04, 0.18018E+04, 0.18067E+04, 0.18117E+04, 0.18167E+04, 0.18216E+04, 0.18266E+04 &
, 0.18316E+04, 0.18365E+04, 0.18415E+04, 0.18465E+04, 0.18515E+04, 0.18565E+04, 0.18615E+04, 0.18665E+04, 0.18715E+04, 0.18765E+04 &
, 0.18815E+04, 0.18865E+04, 0.18915E+04, 0.18966E+04, 0.19016E+04, 0.19066E+04, 0.19117E+04, 0.19167E+04, 0.19217E+04, 0.19268E+04 &
, 0.19318E+04, 0.19369E+04, 0.19420E+04, 0.19470E+04, 0.19521E+04, 0.19572E+04, 0.19623E+04, 0.19674E+04, 0.19724E+04, 0.19775E+04 &
, 0.19826E+04, 0.19877E+04, 0.19928E+04, 0.19979E+04, 0.20031E+04, 0.20082E+04, 0.20133E+04, 0.20184E+04, 0.20236E+04, 0.20287E+04 &
, 0.20338E+04, 0.20390E+04, 0.20441E+04, 0.20493E+04, 0.20545E+04, 0.20596E+04, 0.20648E+04, 0.20700E+04, 0.20751E+04, 0.20803E+04 &
, 0.20855E+04, 0.20907E+04, 0.20959E+04, 0.21011E+04, 0.21063E+04, 0.21115E+04, 0.21167E+04, 0.21219E+04, 0.21272E+04, 0.21324E+04 &
, 0.21376E+04, 0.21429E+04, 0.21481E+04, 0.21533E+04, 0.21586E+04, 0.21639E+04, 0.21691E+04, 0.21744E+04, 0.21796E+04, 0.21849E+04 &
, 0.21902E+04, 0.21955E+04, 0.22008E+04, 0.22061E+04, 0.22114E+04, 0.22167E+04, 0.22220E+04, 0.22273E+04, 0.22326E+04, 0.22379E+04 &
, 0.22432E+04, 0.22486E+04, 0.22539E+04, 0.22593E+04, 0.22646E+04, 0.22699E+04, 0.22753E+04, 0.22807E+04, 0.22860E+04, 0.22914E+04 &
, 0.22968E+04, 0.23021E+04, 0.23075E+04, 0.23129E+04, 0.23183E+04, 0.23237E+04, 0.23291E+04, 0.23345E+04, 0.23399E+04, 0.23454E+04 &
, 0.23508E+04, 0.23562E+04, 0.23616E+04, 0.23671E+04, 0.23725E+04, 0.23780E+04, 0.23834E+04, 0.23889E+04, 0.23943E+04, 0.23998E+04 /
 !            5           5
 DATA(QofT(          50 ,J),J=1,501)/   16.6040000915527       &
, 0.18194E+02, 0.19785E+02, 0.21376E+02, 0.22967E+02, 0.24558E+02, 0.26150E+02, 0.27741E+02, 0.29333E+02, 0.30924E+02, 0.32516E+02 &
, 0.34108E+02, 0.35700E+02, 0.37292E+02, 0.38884E+02, 0.40476E+02, 0.42068E+02, 0.43660E+02, 0.45252E+02, 0.46844E+02, 0.48436E+02 &
, 0.50028E+02, 0.51621E+02, 0.53213E+02, 0.54805E+02, 0.56397E+02, 0.57990E+02, 0.59582E+02, 0.61174E+02, 0.62767E+02, 0.64359E+02 &
, 0.65951E+02, 0.67544E+02, 0.69136E+02, 0.70729E+02, 0.72321E+02, 0.73914E+02, 0.75506E+02, 0.77099E+02, 0.78691E+02, 0.80284E+02 &
, 0.81876E+02, 0.83469E+02, 0.85062E+02, 0.86654E+02, 0.88247E+02, 0.89840E+02, 0.91432E+02, 0.93025E+02, 0.94618E+02, 0.96210E+02 &
, 0.97803E+02, 0.99396E+02, 0.10099E+03, 0.10258E+03, 0.10417E+03, 0.10577E+03, 0.10736E+03, 0.10895E+03, 0.11055E+03, 0.11214E+03 &
, 0.11373E+03, 0.11532E+03, 0.11692E+03, 0.11851E+03, 0.12010E+03, 0.12170E+03, 0.12329E+03, 0.12488E+03, 0.12648E+03, 0.12807E+03 &
, 0.12966E+03, 0.13125E+03, 0.13285E+03, 0.13444E+03, 0.13603E+03, 0.13763E+03, 0.13922E+03, 0.14081E+03, 0.14241E+03, 0.14400E+03 &
, 0.14559E+03, 0.14719E+03, 0.14878E+03, 0.15037E+03, 0.15197E+03, 0.15356E+03, 0.15515E+03, 0.15675E+03, 0.15834E+03, 0.15993E+03 &
, 0.16153E+03, 0.16312E+03, 0.16471E+03, 0.16631E+03, 0.16790E+03, 0.16949E+03, 0.17109E+03, 0.17268E+03, 0.17427E+03, 0.17587E+03 &
, 0.17746E+03, 0.17905E+03, 0.18065E+03, 0.18224E+03, 0.18384E+03, 0.18543E+03, 0.18702E+03, 0.18862E+03, 0.19021E+03, 0.19180E+03 &
, 0.19340E+03, 0.19499E+03, 0.19659E+03, 0.19818E+03, 0.19977E+03, 0.20137E+03, 0.20296E+03, 0.20455E+03, 0.20615E+03, 0.20774E+03 &
, 0.20934E+03, 0.21093E+03, 0.21252E+03, 0.21412E+03, 0.21571E+03, 0.21731E+03, 0.21890E+03, 0.22050E+03, 0.22209E+03, 0.22369E+03 &
, 0.22528E+03, 0.22687E+03, 0.22847E+03, 0.23006E+03, 0.23166E+03, 0.23325E+03, 0.23485E+03, 0.23644E+03, 0.23804E+03, 0.23963E+03 &
, 0.24123E+03, 0.24282E+03, 0.24442E+03, 0.24601E+03, 0.24761E+03, 0.24920E+03, 0.25080E+03, 0.25240E+03, 0.25399E+03, 0.25559E+03 &
, 0.25718E+03, 0.25878E+03, 0.26038E+03, 0.26197E+03, 0.26357E+03, 0.26516E+03, 0.26676E+03, 0.26836E+03, 0.26996E+03, 0.27155E+03 &
, 0.27315E+03, 0.27475E+03, 0.27634E+03, 0.27794E+03, 0.27954E+03, 0.28114E+03, 0.28274E+03, 0.28433E+03, 0.28593E+03, 0.28753E+03 &
, 0.28913E+03, 0.29073E+03, 0.29233E+03, 0.29393E+03, 0.29553E+03, 0.29713E+03, 0.29873E+03, 0.30033E+03, 0.30193E+03, 0.30353E+03 &
, 0.30513E+03, 0.30673E+03, 0.30833E+03, 0.30993E+03, 0.31154E+03, 0.31314E+03, 0.31474E+03, 0.31634E+03, 0.31795E+03, 0.31955E+03 &
, 0.32115E+03, 0.32276E+03, 0.32436E+03, 0.32597E+03, 0.32757E+03, 0.32918E+03, 0.33078E+03, 0.33239E+03, 0.33399E+03, 0.33560E+03 &
, 0.33721E+03, 0.33881E+03, 0.34042E+03, 0.34203E+03, 0.34364E+03, 0.34525E+03, 0.34686E+03, 0.34847E+03, 0.35008E+03, 0.35169E+03 &
, 0.35330E+03, 0.35491E+03, 0.35652E+03, 0.35813E+03, 0.35975E+03, 0.36136E+03, 0.36297E+03, 0.36459E+03, 0.36620E+03, 0.36782E+03 &
, 0.36943E+03, 0.37105E+03, 0.37266E+03, 0.37428E+03, 0.37590E+03, 0.37752E+03, 0.37914E+03, 0.38076E+03, 0.38237E+03, 0.38400E+03 &
, 0.38562E+03, 0.38724E+03, 0.38886E+03, 0.39048E+03, 0.39211E+03, 0.39373E+03, 0.39535E+03, 0.39698E+03, 0.39860E+03, 0.40023E+03 &
, 0.40186E+03, 0.40349E+03, 0.40511E+03, 0.40674E+03, 0.40837E+03, 0.41000E+03, 0.41163E+03, 0.41327E+03, 0.41490E+03, 0.41653E+03 &
, 0.41817E+03, 0.41980E+03, 0.42144E+03, 0.42307E+03, 0.42471E+03, 0.42635E+03, 0.42798E+03, 0.42962E+03, 0.43126E+03, 0.43290E+03 &
, 0.43455E+03, 0.43619E+03, 0.43783E+03, 0.43947E+03, 0.44112E+03, 0.44276E+03, 0.44441E+03, 0.44606E+03, 0.44771E+03, 0.44936E+03 &
, 0.45101E+03, 0.45266E+03, 0.45431E+03, 0.45596E+03, 0.45761E+03, 0.45927E+03, 0.46092E+03, 0.46258E+03, 0.46424E+03, 0.46589E+03 &
, 0.46755E+03, 0.46921E+03, 0.47087E+03, 0.47253E+03, 0.47420E+03, 0.47586E+03, 0.47753E+03, 0.47919E+03, 0.48086E+03, 0.48253E+03 &
, 0.48419E+03, 0.48586E+03, 0.48753E+03, 0.48921E+03, 0.49088E+03, 0.49255E+03, 0.49423E+03, 0.49590E+03, 0.49758E+03, 0.49926E+03 &
, 0.50094E+03, 0.50262E+03, 0.50430E+03, 0.50598E+03, 0.50766E+03, 0.50935E+03, 0.51103E+03, 0.51272E+03, 0.51441E+03, 0.51610E+03 &
, 0.51779E+03, 0.51948E+03, 0.52117E+03, 0.52286E+03, 0.52456E+03, 0.52625E+03, 0.52795E+03, 0.52965E+03, 0.53135E+03, 0.53305E+03 &
, 0.53475E+03, 0.53645E+03, 0.53816E+03, 0.53986E+03, 0.54157E+03, 0.54328E+03, 0.54499E+03, 0.54670E+03, 0.54841E+03, 0.55012E+03 &
, 0.55184E+03, 0.55355E+03, 0.55527E+03, 0.55699E+03, 0.55870E+03, 0.56042E+03, 0.56215E+03, 0.56387E+03, 0.56559E+03, 0.56732E+03 &
, 0.56905E+03, 0.57077E+03, 0.57250E+03, 0.57424E+03, 0.57597E+03, 0.57770E+03, 0.57944E+03, 0.58117E+03, 0.58291E+03, 0.58465E+03 &
, 0.58639E+03, 0.58813E+03, 0.58987E+03, 0.59162E+03, 0.59336E+03, 0.59511E+03, 0.59686E+03, 0.59861E+03, 0.60036E+03, 0.60211E+03 &
, 0.60387E+03, 0.60562E+03, 0.60738E+03, 0.60914E+03, 0.61090E+03, 0.61266E+03, 0.61442E+03, 0.61619E+03, 0.61795E+03, 0.61972E+03 &
, 0.62149E+03, 0.62326E+03, 0.62503E+03, 0.62680E+03, 0.62858E+03, 0.63036E+03, 0.63213E+03, 0.63391E+03, 0.63569E+03, 0.63747E+03 &
, 0.63926E+03, 0.64104E+03, 0.64283E+03, 0.64462E+03, 0.64641E+03, 0.64820E+03, 0.64999E+03, 0.65178E+03, 0.65358E+03, 0.65538E+03 &
, 0.65718E+03, 0.65898E+03, 0.66078E+03, 0.66258E+03, 0.66439E+03, 0.66619E+03, 0.66800E+03, 0.66981E+03, 0.67162E+03, 0.67344E+03 &
, 0.67525E+03, 0.67707E+03, 0.67888E+03, 0.68070E+03, 0.68253E+03, 0.68435E+03, 0.68617E+03, 0.68800E+03, 0.68982E+03, 0.69165E+03 &
, 0.69348E+03, 0.69532E+03, 0.69715E+03, 0.69899E+03, 0.70082E+03, 0.70266E+03, 0.70450E+03, 0.70634E+03, 0.70819E+03, 0.71003E+03 &
, 0.71188E+03, 0.71373E+03, 0.71558E+03, 0.71743E+03, 0.71929E+03, 0.72114E+03, 0.72300E+03, 0.72486E+03, 0.72672E+03, 0.72858E+03 &
, 0.73044E+03, 0.73231E+03, 0.73418E+03, 0.73605E+03, 0.73792E+03, 0.73979E+03, 0.74166E+03, 0.74354E+03, 0.74542E+03, 0.74730E+03 &
, 0.74918E+03, 0.75106E+03, 0.75294E+03, 0.75483E+03, 0.75672E+03, 0.75861E+03, 0.76050E+03, 0.76239E+03, 0.76429E+03, 0.76618E+03 &
, 0.76808E+03, 0.76998E+03, 0.77188E+03, 0.77379E+03, 0.77569E+03, 0.77760E+03, 0.77951E+03, 0.78142E+03, 0.78333E+03, 0.78525E+03 &
, 0.78716E+03, 0.78908E+03, 0.79100E+03, 0.79292E+03, 0.79485E+03, 0.79677E+03, 0.79870E+03, 0.80063E+03, 0.80256E+03, 0.80449E+03 &
, 0.80642E+03, 0.80836E+03, 0.81030E+03, 0.81224E+03, 0.81418E+03, 0.81612E+03, 0.81807E+03, 0.82001E+03, 0.82196E+03, 0.82391E+03 &
, 0.82587E+03, 0.82782E+03, 0.82978E+03, 0.83173E+03, 0.83369E+03, 0.83565E+03, 0.83762E+03, 0.83958E+03, 0.84155E+03, 0.84352E+03 &
, 0.84549E+03, 0.84746E+03, 0.84944E+03, 0.85141E+03, 0.85339E+03, 0.85537E+03, 0.85735E+03, 0.85934E+03, 0.86132E+03, 0.86331E+03 /
 !            5           6
 DATA(QofT(          51 ,J),J=1,501)/   97.3300018310547       &
, 0.10664E+03, 0.11596E+03, 0.12527E+03, 0.13459E+03, 0.14391E+03, 0.15323E+03, 0.16255E+03, 0.17187E+03, 0.18119E+03, 0.19051E+03 &
, 0.19983E+03, 0.20915E+03, 0.21847E+03, 0.22780E+03, 0.23712E+03, 0.24644E+03, 0.25576E+03, 0.26509E+03, 0.27441E+03, 0.28373E+03 &
, 0.29306E+03, 0.30238E+03, 0.31170E+03, 0.32103E+03, 0.33035E+03, 0.33968E+03, 0.34900E+03, 0.35832E+03, 0.36765E+03, 0.37697E+03 &
, 0.38630E+03, 0.39562E+03, 0.40495E+03, 0.41427E+03, 0.42360E+03, 0.43292E+03, 0.44225E+03, 0.45158E+03, 0.46090E+03, 0.47023E+03 &
, 0.47955E+03, 0.48888E+03, 0.49821E+03, 0.50753E+03, 0.51686E+03, 0.52618E+03, 0.53551E+03, 0.54484E+03, 0.55416E+03, 0.56349E+03 &
, 0.57282E+03, 0.58214E+03, 0.59147E+03, 0.60080E+03, 0.61013E+03, 0.61945E+03, 0.62878E+03, 0.63811E+03, 0.64744E+03, 0.65676E+03 &
, 0.66609E+03, 0.67542E+03, 0.68475E+03, 0.69408E+03, 0.70340E+03, 0.71273E+03, 0.72206E+03, 0.73139E+03, 0.74072E+03, 0.75005E+03 &
, 0.75938E+03, 0.76870E+03, 0.77803E+03, 0.78736E+03, 0.79669E+03, 0.80602E+03, 0.81535E+03, 0.82468E+03, 0.83401E+03, 0.84334E+03 &
, 0.85267E+03, 0.86200E+03, 0.87133E+03, 0.88066E+03, 0.88999E+03, 0.89932E+03, 0.90865E+03, 0.91798E+03, 0.92731E+03, 0.93664E+03 &
, 0.94597E+03, 0.95530E+03, 0.96463E+03, 0.97396E+03, 0.98329E+03, 0.99262E+03, 0.10020E+04, 0.10113E+04, 0.10206E+04, 0.10299E+04 &
, 0.10393E+04, 0.10486E+04, 0.10579E+04, 0.10673E+04, 0.10766E+04, 0.10859E+04, 0.10953E+04, 0.11046E+04, 0.11139E+04, 0.11233E+04 &
, 0.11326E+04, 0.11419E+04, 0.11513E+04, 0.11606E+04, 0.11699E+04, 0.11793E+04, 0.11886E+04, 0.11979E+04, 0.12073E+04, 0.12166E+04 &
, 0.12259E+04, 0.12353E+04, 0.12446E+04, 0.12539E+04, 0.12633E+04, 0.12726E+04, 0.12820E+04, 0.12913E+04, 0.13006E+04, 0.13100E+04 &
, 0.13193E+04, 0.13286E+04, 0.13380E+04, 0.13473E+04, 0.13566E+04, 0.13660E+04, 0.13753E+04, 0.13847E+04, 0.13940E+04, 0.14033E+04 &
, 0.14127E+04, 0.14220E+04, 0.14314E+04, 0.14407E+04, 0.14500E+04, 0.14594E+04, 0.14687E+04, 0.14781E+04, 0.14874E+04, 0.14968E+04 &
, 0.15061E+04, 0.15155E+04, 0.15248E+04, 0.15341E+04, 0.15435E+04, 0.15528E+04, 0.15622E+04, 0.15715E+04, 0.15809E+04, 0.15902E+04 &
, 0.15996E+04, 0.16089E+04, 0.16183E+04, 0.16276E+04, 0.16370E+04, 0.16464E+04, 0.16557E+04, 0.16651E+04, 0.16744E+04, 0.16838E+04 &
, 0.16932E+04, 0.17025E+04, 0.17119E+04, 0.17212E+04, 0.17306E+04, 0.17400E+04, 0.17493E+04, 0.17587E+04, 0.17681E+04, 0.17774E+04 &
, 0.17868E+04, 0.17962E+04, 0.18056E+04, 0.18149E+04, 0.18243E+04, 0.18337E+04, 0.18431E+04, 0.18525E+04, 0.18618E+04, 0.18712E+04 &
, 0.18806E+04, 0.18900E+04, 0.18994E+04, 0.19088E+04, 0.19182E+04, 0.19276E+04, 0.19370E+04, 0.19464E+04, 0.19558E+04, 0.19652E+04 &
, 0.19746E+04, 0.19840E+04, 0.19934E+04, 0.20028E+04, 0.20122E+04, 0.20216E+04, 0.20311E+04, 0.20405E+04, 0.20499E+04, 0.20593E+04 &
, 0.20687E+04, 0.20782E+04, 0.20876E+04, 0.20970E+04, 0.21065E+04, 0.21159E+04, 0.21254E+04, 0.21348E+04, 0.21443E+04, 0.21537E+04 &
, 0.21632E+04, 0.21726E+04, 0.21821E+04, 0.21915E+04, 0.22010E+04, 0.22105E+04, 0.22199E+04, 0.22294E+04, 0.22389E+04, 0.22484E+04 &
, 0.22578E+04, 0.22673E+04, 0.22768E+04, 0.22863E+04, 0.22958E+04, 0.23053E+04, 0.23148E+04, 0.23243E+04, 0.23338E+04, 0.23433E+04 &
, 0.23529E+04, 0.23624E+04, 0.23719E+04, 0.23814E+04, 0.23910E+04, 0.24005E+04, 0.24100E+04, 0.24196E+04, 0.24291E+04, 0.24387E+04 &
, 0.24482E+04, 0.24578E+04, 0.24674E+04, 0.24769E+04, 0.24865E+04, 0.24961E+04, 0.25056E+04, 0.25152E+04, 0.25248E+04, 0.25344E+04 &
, 0.25440E+04, 0.25536E+04, 0.25632E+04, 0.25728E+04, 0.25824E+04, 0.25921E+04, 0.26017E+04, 0.26113E+04, 0.26210E+04, 0.26306E+04 &
, 0.26402E+04, 0.26499E+04, 0.26595E+04, 0.26692E+04, 0.26789E+04, 0.26885E+04, 0.26982E+04, 0.27079E+04, 0.27176E+04, 0.27273E+04 &
, 0.27370E+04, 0.27467E+04, 0.27564E+04, 0.27661E+04, 0.27758E+04, 0.27855E+04, 0.27952E+04, 0.28050E+04, 0.28147E+04, 0.28245E+04 &
, 0.28342E+04, 0.28440E+04, 0.28537E+04, 0.28635E+04, 0.28733E+04, 0.28831E+04, 0.28928E+04, 0.29026E+04, 0.29124E+04, 0.29222E+04 &
, 0.29320E+04, 0.29419E+04, 0.29517E+04, 0.29615E+04, 0.29713E+04, 0.29812E+04, 0.29910E+04, 0.30009E+04, 0.30107E+04, 0.30206E+04 &
, 0.30305E+04, 0.30403E+04, 0.30502E+04, 0.30601E+04, 0.30700E+04, 0.30799E+04, 0.30898E+04, 0.30998E+04, 0.31097E+04, 0.31196E+04 &
, 0.31296E+04, 0.31395E+04, 0.31494E+04, 0.31594E+04, 0.31694E+04, 0.31793E+04, 0.31893E+04, 0.31993E+04, 0.32093E+04, 0.32193E+04 &
, 0.32293E+04, 0.32393E+04, 0.32494E+04, 0.32594E+04, 0.32694E+04, 0.32795E+04, 0.32895E+04, 0.32996E+04, 0.33096E+04, 0.33197E+04 &
, 0.33298E+04, 0.33399E+04, 0.33500E+04, 0.33601E+04, 0.33702E+04, 0.33803E+04, 0.33904E+04, 0.34006E+04, 0.34107E+04, 0.34209E+04 &
, 0.34310E+04, 0.34412E+04, 0.34514E+04, 0.34615E+04, 0.34717E+04, 0.34819E+04, 0.34921E+04, 0.35023E+04, 0.35126E+04, 0.35228E+04 &
, 0.35330E+04, 0.35433E+04, 0.35535E+04, 0.35638E+04, 0.35741E+04, 0.35843E+04, 0.35946E+04, 0.36049E+04, 0.36152E+04, 0.36255E+04 &
, 0.36358E+04, 0.36462E+04, 0.36565E+04, 0.36668E+04, 0.36772E+04, 0.36876E+04, 0.36979E+04, 0.37083E+04, 0.37187E+04, 0.37291E+04 &
, 0.37395E+04, 0.37499E+04, 0.37603E+04, 0.37708E+04, 0.37812E+04, 0.37916E+04, 0.38021E+04, 0.38126E+04, 0.38230E+04, 0.38335E+04 &
, 0.38440E+04, 0.38545E+04, 0.38650E+04, 0.38755E+04, 0.38861E+04, 0.38966E+04, 0.39071E+04, 0.39177E+04, 0.39283E+04, 0.39388E+04 &
, 0.39494E+04, 0.39600E+04, 0.39706E+04, 0.39812E+04, 0.39918E+04, 0.40025E+04, 0.40131E+04, 0.40237E+04, 0.40344E+04, 0.40451E+04 &
, 0.40557E+04, 0.40664E+04, 0.40771E+04, 0.40878E+04, 0.40985E+04, 0.41092E+04, 0.41200E+04, 0.41307E+04, 0.41414E+04, 0.41522E+04 &
, 0.41630E+04, 0.41737E+04, 0.41845E+04, 0.41953E+04, 0.42061E+04, 0.42169E+04, 0.42278E+04, 0.42386E+04, 0.42495E+04, 0.42603E+04 &
, 0.42712E+04, 0.42820E+04, 0.42929E+04, 0.43038E+04, 0.43147E+04, 0.43256E+04, 0.43366E+04, 0.43475E+04, 0.43584E+04, 0.43694E+04 &
, 0.43803E+04, 0.43913E+04, 0.44023E+04, 0.44133E+04, 0.44243E+04, 0.44353E+04, 0.44463E+04, 0.44573E+04, 0.44684E+04, 0.44794E+04 &
, 0.44905E+04, 0.45016E+04, 0.45127E+04, 0.45237E+04, 0.45348E+04, 0.45460E+04, 0.45571E+04, 0.45682E+04, 0.45794E+04, 0.45905E+04 &
, 0.46017E+04, 0.46128E+04, 0.46240E+04, 0.46352E+04, 0.46464E+04, 0.46576E+04, 0.46689E+04, 0.46801E+04, 0.46913E+04, 0.47026E+04 &
, 0.47139E+04, 0.47251E+04, 0.47364E+04, 0.47477E+04, 0.47590E+04, 0.47704E+04, 0.47817E+04, 0.47930E+04, 0.48044E+04, 0.48157E+04 &
, 0.48271E+04, 0.48385E+04, 0.48499E+04, 0.48613E+04, 0.48727E+04, 0.48841E+04, 0.48955E+04, 0.49070E+04, 0.49184E+04, 0.49299E+04 &
, 0.49414E+04, 0.49529E+04, 0.49644E+04, 0.49759E+04, 0.49874E+04, 0.49989E+04, 0.50105E+04, 0.50220E+04, 0.50336E+04, 0.50452E+04 /
 !            5           7
 DATA(QofT(          52 ,J),J=1,501)/   8.21469974517822       &
, 0.90012E+01, 0.97879E+01, 0.10575E+02, 0.11361E+02, 0.12148E+02, 0.12935E+02, 0.13722E+02, 0.14510E+02, 0.15297E+02, 0.16084E+02 &
, 0.16871E+02, 0.17658E+02, 0.18445E+02, 0.19233E+02, 0.20020E+02, 0.20807E+02, 0.21594E+02, 0.22382E+02, 0.23169E+02, 0.23956E+02 &
, 0.24744E+02, 0.25531E+02, 0.26319E+02, 0.27106E+02, 0.27893E+02, 0.28681E+02, 0.29468E+02, 0.30256E+02, 0.31043E+02, 0.31831E+02 &
, 0.32618E+02, 0.33406E+02, 0.34193E+02, 0.34981E+02, 0.35768E+02, 0.36556E+02, 0.37343E+02, 0.38131E+02, 0.38918E+02, 0.39706E+02 &
, 0.40493E+02, 0.41281E+02, 0.42069E+02, 0.42856E+02, 0.43644E+02, 0.44431E+02, 0.45219E+02, 0.46007E+02, 0.46794E+02, 0.47582E+02 &
, 0.48370E+02, 0.49157E+02, 0.49945E+02, 0.50733E+02, 0.51520E+02, 0.52308E+02, 0.53096E+02, 0.53883E+02, 0.54671E+02, 0.55459E+02 &
, 0.56247E+02, 0.57034E+02, 0.57822E+02, 0.58610E+02, 0.59398E+02, 0.60185E+02, 0.60973E+02, 0.61761E+02, 0.62549E+02, 0.63336E+02 &
, 0.64124E+02, 0.64912E+02, 0.65700E+02, 0.66488E+02, 0.67276E+02, 0.68063E+02, 0.68851E+02, 0.69639E+02, 0.70427E+02, 0.71215E+02 &
, 0.72003E+02, 0.72791E+02, 0.73578E+02, 0.74366E+02, 0.75154E+02, 0.75942E+02, 0.76730E+02, 0.77518E+02, 0.78306E+02, 0.79094E+02 &
, 0.79882E+02, 0.80670E+02, 0.81458E+02, 0.82246E+02, 0.83034E+02, 0.83822E+02, 0.84610E+02, 0.85398E+02, 0.86186E+02, 0.86974E+02 &
, 0.87762E+02, 0.88550E+02, 0.89338E+02, 0.90126E+02, 0.90914E+02, 0.91702E+02, 0.92491E+02, 0.93279E+02, 0.94067E+02, 0.94855E+02 &
, 0.95643E+02, 0.96431E+02, 0.97219E+02, 0.98008E+02, 0.98796E+02, 0.99584E+02, 0.10037E+03, 0.10116E+03, 0.10195E+03, 0.10274E+03 &
, 0.10353E+03, 0.10431E+03, 0.10510E+03, 0.10589E+03, 0.10668E+03, 0.10747E+03, 0.10826E+03, 0.10904E+03, 0.10983E+03, 0.11062E+03 &
, 0.11141E+03, 0.11220E+03, 0.11299E+03, 0.11378E+03, 0.11456E+03, 0.11535E+03, 0.11614E+03, 0.11693E+03, 0.11772E+03, 0.11851E+03 &
, 0.11930E+03, 0.12008E+03, 0.12087E+03, 0.12166E+03, 0.12245E+03, 0.12324E+03, 0.12403E+03, 0.12482E+03, 0.12561E+03, 0.12640E+03 &
, 0.12719E+03, 0.12798E+03, 0.12876E+03, 0.12955E+03, 0.13034E+03, 0.13113E+03, 0.13192E+03, 0.13271E+03, 0.13350E+03, 0.13429E+03 &
, 0.13508E+03, 0.13587E+03, 0.13666E+03, 0.13745E+03, 0.13824E+03, 0.13903E+03, 0.13982E+03, 0.14061E+03, 0.14140E+03, 0.14219E+03 &
, 0.14298E+03, 0.14377E+03, 0.14456E+03, 0.14536E+03, 0.14615E+03, 0.14694E+03, 0.14773E+03, 0.14852E+03, 0.14931E+03, 0.15010E+03 &
, 0.15089E+03, 0.15169E+03, 0.15248E+03, 0.15327E+03, 0.15406E+03, 0.15485E+03, 0.15565E+03, 0.15644E+03, 0.15723E+03, 0.15802E+03 &
, 0.15882E+03, 0.15961E+03, 0.16040E+03, 0.16120E+03, 0.16199E+03, 0.16278E+03, 0.16358E+03, 0.16437E+03, 0.16517E+03, 0.16596E+03 &
, 0.16675E+03, 0.16755E+03, 0.16834E+03, 0.16914E+03, 0.16993E+03, 0.17073E+03, 0.17153E+03, 0.17232E+03, 0.17312E+03, 0.17391E+03 &
, 0.17471E+03, 0.17551E+03, 0.17630E+03, 0.17710E+03, 0.17790E+03, 0.17869E+03, 0.17949E+03, 0.18029E+03, 0.18109E+03, 0.18189E+03 &
, 0.18269E+03, 0.18348E+03, 0.18428E+03, 0.18508E+03, 0.18588E+03, 0.18668E+03, 0.18748E+03, 0.18828E+03, 0.18908E+03, 0.18988E+03 &
, 0.19069E+03, 0.19149E+03, 0.19229E+03, 0.19309E+03, 0.19389E+03, 0.19470E+03, 0.19550E+03, 0.19630E+03, 0.19711E+03, 0.19791E+03 &
, 0.19871E+03, 0.19952E+03, 0.20032E+03, 0.20113E+03, 0.20193E+03, 0.20274E+03, 0.20355E+03, 0.20435E+03, 0.20516E+03, 0.20597E+03 &
, 0.20677E+03, 0.20758E+03, 0.20839E+03, 0.20920E+03, 0.21001E+03, 0.21082E+03, 0.21163E+03, 0.21244E+03, 0.21325E+03, 0.21406E+03 &
, 0.21487E+03, 0.21568E+03, 0.21649E+03, 0.21730E+03, 0.21812E+03, 0.21893E+03, 0.21974E+03, 0.22056E+03, 0.22137E+03, 0.22219E+03 &
, 0.22300E+03, 0.22382E+03, 0.22463E+03, 0.22545E+03, 0.22627E+03, 0.22708E+03, 0.22790E+03, 0.22872E+03, 0.22954E+03, 0.23036E+03 &
, 0.23118E+03, 0.23200E+03, 0.23282E+03, 0.23364E+03, 0.23446E+03, 0.23528E+03, 0.23611E+03, 0.23693E+03, 0.23775E+03, 0.23858E+03 &
, 0.23940E+03, 0.24022E+03, 0.24105E+03, 0.24188E+03, 0.24270E+03, 0.24353E+03, 0.24436E+03, 0.24518E+03, 0.24601E+03, 0.24684E+03 &
, 0.24767E+03, 0.24850E+03, 0.24933E+03, 0.25016E+03, 0.25099E+03, 0.25183E+03, 0.25266E+03, 0.25349E+03, 0.25433E+03, 0.25516E+03 &
, 0.25599E+03, 0.25683E+03, 0.25767E+03, 0.25850E+03, 0.25934E+03, 0.26018E+03, 0.26102E+03, 0.26185E+03, 0.26269E+03, 0.26353E+03 &
, 0.26437E+03, 0.26522E+03, 0.26606E+03, 0.26690E+03, 0.26774E+03, 0.26859E+03, 0.26943E+03, 0.27027E+03, 0.27112E+03, 0.27197E+03 &
, 0.27281E+03, 0.27366E+03, 0.27451E+03, 0.27536E+03, 0.27620E+03, 0.27705E+03, 0.27790E+03, 0.27875E+03, 0.27961E+03, 0.28046E+03 &
, 0.28131E+03, 0.28216E+03, 0.28302E+03, 0.28387E+03, 0.28473E+03, 0.28559E+03, 0.28644E+03, 0.28730E+03, 0.28816E+03, 0.28902E+03 &
, 0.28988E+03, 0.29074E+03, 0.29160E+03, 0.29246E+03, 0.29332E+03, 0.29418E+03, 0.29505E+03, 0.29591E+03, 0.29677E+03, 0.29764E+03 &
, 0.29851E+03, 0.29937E+03, 0.30024E+03, 0.30111E+03, 0.30198E+03, 0.30285E+03, 0.30372E+03, 0.30459E+03, 0.30546E+03, 0.30633E+03 &
, 0.30721E+03, 0.30808E+03, 0.30895E+03, 0.30983E+03, 0.31071E+03, 0.31158E+03, 0.31246E+03, 0.31334E+03, 0.31422E+03, 0.31510E+03 &
, 0.31598E+03, 0.31686E+03, 0.31774E+03, 0.31862E+03, 0.31951E+03, 0.32039E+03, 0.32128E+03, 0.32216E+03, 0.32305E+03, 0.32393E+03 &
, 0.32482E+03, 0.32571E+03, 0.32660E+03, 0.32749E+03, 0.32838E+03, 0.32927E+03, 0.33017E+03, 0.33106E+03, 0.33195E+03, 0.33285E+03 &
, 0.33374E+03, 0.33464E+03, 0.33554E+03, 0.33643E+03, 0.33733E+03, 0.33823E+03, 0.33913E+03, 0.34003E+03, 0.34094E+03, 0.34184E+03 &
, 0.34274E+03, 0.34365E+03, 0.34455E+03, 0.34546E+03, 0.34636E+03, 0.34727E+03, 0.34818E+03, 0.34909E+03, 0.35000E+03, 0.35091E+03 &
, 0.35182E+03, 0.35273E+03, 0.35365E+03, 0.35456E+03, 0.35548E+03, 0.35639E+03, 0.35731E+03, 0.35823E+03, 0.35914E+03, 0.36006E+03 &
, 0.36098E+03, 0.36190E+03, 0.36282E+03, 0.36375E+03, 0.36467E+03, 0.36559E+03, 0.36652E+03, 0.36744E+03, 0.36837E+03, 0.36930E+03 &
, 0.37023E+03, 0.37115E+03, 0.37208E+03, 0.37302E+03, 0.37395E+03, 0.37488E+03, 0.37581E+03, 0.37675E+03, 0.37768E+03, 0.37862E+03 &
, 0.37955E+03, 0.38049E+03, 0.38143E+03, 0.38237E+03, 0.38331E+03, 0.38425E+03, 0.38519E+03, 0.38613E+03, 0.38708E+03, 0.38802E+03 &
, 0.38897E+03, 0.38991E+03, 0.39086E+03, 0.39181E+03, 0.39276E+03, 0.39371E+03, 0.39466E+03, 0.39561E+03, 0.39656E+03, 0.39752E+03 &
, 0.39847E+03, 0.39943E+03, 0.40038E+03, 0.40134E+03, 0.40230E+03, 0.40325E+03, 0.40421E+03, 0.40517E+03, 0.40614E+03, 0.40710E+03 &
, 0.40806E+03, 0.40902E+03, 0.40999E+03, 0.41096E+03, 0.41192E+03, 0.41289E+03, 0.41386E+03, 0.41483E+03, 0.41580E+03, 0.41677E+03 &
, 0.41774E+03, 0.41871E+03, 0.41969E+03, 0.42066E+03, 0.42164E+03, 0.42262E+03, 0.42359E+03, 0.42457E+03, 0.42555E+03, 0.42653E+03 /
 !            5           8
 DATA(QofT(          53 ,J),J=1,501)/   8.64490032196045       &
, 0.94745E+01, 0.10304E+02, 0.11134E+02, 0.11964E+02, 0.12794E+02, 0.13624E+02, 0.14454E+02, 0.15284E+02, 0.16115E+02, 0.16945E+02 &
, 0.17775E+02, 0.18605E+02, 0.19436E+02, 0.20266E+02, 0.21096E+02, 0.21927E+02, 0.22757E+02, 0.23588E+02, 0.24418E+02, 0.25248E+02 &
, 0.26079E+02, 0.26909E+02, 0.27740E+02, 0.28570E+02, 0.29401E+02, 0.30231E+02, 0.31062E+02, 0.31892E+02, 0.32723E+02, 0.33553E+02 &
, 0.34384E+02, 0.35214E+02, 0.36045E+02, 0.36876E+02, 0.37706E+02, 0.38537E+02, 0.39367E+02, 0.40198E+02, 0.41029E+02, 0.41859E+02 &
, 0.42690E+02, 0.43521E+02, 0.44351E+02, 0.45182E+02, 0.46013E+02, 0.46843E+02, 0.47674E+02, 0.48505E+02, 0.49336E+02, 0.50166E+02 &
, 0.50997E+02, 0.51828E+02, 0.52659E+02, 0.53489E+02, 0.54320E+02, 0.55151E+02, 0.55982E+02, 0.56812E+02, 0.57643E+02, 0.58474E+02 &
, 0.59305E+02, 0.60136E+02, 0.60966E+02, 0.61797E+02, 0.62628E+02, 0.63459E+02, 0.64290E+02, 0.65121E+02, 0.65952E+02, 0.66783E+02 &
, 0.67613E+02, 0.68444E+02, 0.69275E+02, 0.70106E+02, 0.70937E+02, 0.71768E+02, 0.72599E+02, 0.73430E+02, 0.74261E+02, 0.75092E+02 &
, 0.75923E+02, 0.76754E+02, 0.77585E+02, 0.78416E+02, 0.79247E+02, 0.80078E+02, 0.80909E+02, 0.81740E+02, 0.82571E+02, 0.83402E+02 &
, 0.84233E+02, 0.85064E+02, 0.85895E+02, 0.86726E+02, 0.87558E+02, 0.88389E+02, 0.89220E+02, 0.90051E+02, 0.90882E+02, 0.91713E+02 &
, 0.92544E+02, 0.93376E+02, 0.94207E+02, 0.95038E+02, 0.95869E+02, 0.96700E+02, 0.97532E+02, 0.98363E+02, 0.99194E+02, 0.10003E+03 &
, 0.10086E+03, 0.10169E+03, 0.10252E+03, 0.10335E+03, 0.10418E+03, 0.10501E+03, 0.10584E+03, 0.10668E+03, 0.10751E+03, 0.10834E+03 &
, 0.10917E+03, 0.11000E+03, 0.11083E+03, 0.11166E+03, 0.11250E+03, 0.11333E+03, 0.11416E+03, 0.11499E+03, 0.11582E+03, 0.11665E+03 &
, 0.11749E+03, 0.11832E+03, 0.11915E+03, 0.11998E+03, 0.12081E+03, 0.12165E+03, 0.12248E+03, 0.12331E+03, 0.12414E+03, 0.12497E+03 &
, 0.12581E+03, 0.12664E+03, 0.12747E+03, 0.12830E+03, 0.12913E+03, 0.12997E+03, 0.13080E+03, 0.13163E+03, 0.13246E+03, 0.13330E+03 &
, 0.13413E+03, 0.13496E+03, 0.13579E+03, 0.13663E+03, 0.13746E+03, 0.13829E+03, 0.13913E+03, 0.13996E+03, 0.14079E+03, 0.14162E+03 &
, 0.14246E+03, 0.14329E+03, 0.14413E+03, 0.14496E+03, 0.14579E+03, 0.14663E+03, 0.14746E+03, 0.14829E+03, 0.14913E+03, 0.14996E+03 &
, 0.15080E+03, 0.15163E+03, 0.15247E+03, 0.15330E+03, 0.15413E+03, 0.15497E+03, 0.15580E+03, 0.15664E+03, 0.15747E+03, 0.15831E+03 &
, 0.15915E+03, 0.15998E+03, 0.16082E+03, 0.16165E+03, 0.16249E+03, 0.16333E+03, 0.16416E+03, 0.16500E+03, 0.16584E+03, 0.16667E+03 &
, 0.16751E+03, 0.16835E+03, 0.16918E+03, 0.17002E+03, 0.17086E+03, 0.17170E+03, 0.17254E+03, 0.17337E+03, 0.17421E+03, 0.17505E+03 &
, 0.17589E+03, 0.17673E+03, 0.17757E+03, 0.17841E+03, 0.17925E+03, 0.18009E+03, 0.18093E+03, 0.18177E+03, 0.18261E+03, 0.18345E+03 &
, 0.18429E+03, 0.18513E+03, 0.18598E+03, 0.18682E+03, 0.18766E+03, 0.18850E+03, 0.18934E+03, 0.19019E+03, 0.19103E+03, 0.19187E+03 &
, 0.19272E+03, 0.19356E+03, 0.19441E+03, 0.19525E+03, 0.19610E+03, 0.19694E+03, 0.19779E+03, 0.19863E+03, 0.19948E+03, 0.20033E+03 &
, 0.20117E+03, 0.20202E+03, 0.20287E+03, 0.20372E+03, 0.20457E+03, 0.20541E+03, 0.20626E+03, 0.20711E+03, 0.20796E+03, 0.20881E+03 &
, 0.20966E+03, 0.21051E+03, 0.21136E+03, 0.21222E+03, 0.21307E+03, 0.21392E+03, 0.21477E+03, 0.21563E+03, 0.21648E+03, 0.21733E+03 &
, 0.21819E+03, 0.21904E+03, 0.21990E+03, 0.22075E+03, 0.22161E+03, 0.22246E+03, 0.22332E+03, 0.22418E+03, 0.22504E+03, 0.22589E+03 &
, 0.22675E+03, 0.22761E+03, 0.22847E+03, 0.22933E+03, 0.23019E+03, 0.23105E+03, 0.23191E+03, 0.23278E+03, 0.23364E+03, 0.23450E+03 &
, 0.23536E+03, 0.23623E+03, 0.23709E+03, 0.23796E+03, 0.23882E+03, 0.23969E+03, 0.24055E+03, 0.24142E+03, 0.24229E+03, 0.24315E+03 &
, 0.24402E+03, 0.24489E+03, 0.24576E+03, 0.24663E+03, 0.24750E+03, 0.24837E+03, 0.24924E+03, 0.25012E+03, 0.25099E+03, 0.25186E+03 &
, 0.25273E+03, 0.25361E+03, 0.25448E+03, 0.25536E+03, 0.25623E+03, 0.25711E+03, 0.25799E+03, 0.25887E+03, 0.25974E+03, 0.26062E+03 &
, 0.26150E+03, 0.26238E+03, 0.26326E+03, 0.26414E+03, 0.26503E+03, 0.26591E+03, 0.26679E+03, 0.26768E+03, 0.26856E+03, 0.26944E+03 &
, 0.27033E+03, 0.27122E+03, 0.27210E+03, 0.27299E+03, 0.27388E+03, 0.27477E+03, 0.27566E+03, 0.27655E+03, 0.27744E+03, 0.27833E+03 &
, 0.27922E+03, 0.28011E+03, 0.28101E+03, 0.28190E+03, 0.28279E+03, 0.28369E+03, 0.28459E+03, 0.28548E+03, 0.28638E+03, 0.28728E+03 &
, 0.28818E+03, 0.28908E+03, 0.28998E+03, 0.29088E+03, 0.29178E+03, 0.29268E+03, 0.29358E+03, 0.29449E+03, 0.29539E+03, 0.29630E+03 &
, 0.29720E+03, 0.29811E+03, 0.29902E+03, 0.29992E+03, 0.30083E+03, 0.30174E+03, 0.30265E+03, 0.30356E+03, 0.30448E+03, 0.30539E+03 &
, 0.30630E+03, 0.30721E+03, 0.30813E+03, 0.30904E+03, 0.30996E+03, 0.31088E+03, 0.31180E+03, 0.31271E+03, 0.31363E+03, 0.31455E+03 &
, 0.31547E+03, 0.31640E+03, 0.31732E+03, 0.31824E+03, 0.31916E+03, 0.32009E+03, 0.32101E+03, 0.32194E+03, 0.32287E+03, 0.32380E+03 &
, 0.32472E+03, 0.32565E+03, 0.32658E+03, 0.32752E+03, 0.32845E+03, 0.32938E+03, 0.33031E+03, 0.33125E+03, 0.33218E+03, 0.33312E+03 &
, 0.33405E+03, 0.33499E+03, 0.33593E+03, 0.33687E+03, 0.33781E+03, 0.33875E+03, 0.33969E+03, 0.34064E+03, 0.34158E+03, 0.34252E+03 &
, 0.34347E+03, 0.34441E+03, 0.34536E+03, 0.34631E+03, 0.34726E+03, 0.34821E+03, 0.34916E+03, 0.35011E+03, 0.35106E+03, 0.35201E+03 &
, 0.35296E+03, 0.35392E+03, 0.35487E+03, 0.35583E+03, 0.35679E+03, 0.35775E+03, 0.35870E+03, 0.35966E+03, 0.36062E+03, 0.36159E+03 &
, 0.36255E+03, 0.36351E+03, 0.36448E+03, 0.36544E+03, 0.36641E+03, 0.36737E+03, 0.36834E+03, 0.36931E+03, 0.37028E+03, 0.37125E+03 &
, 0.37222E+03, 0.37319E+03, 0.37417E+03, 0.37514E+03, 0.37611E+03, 0.37709E+03, 0.37807E+03, 0.37904E+03, 0.38002E+03, 0.38100E+03 &
, 0.38198E+03, 0.38296E+03, 0.38395E+03, 0.38493E+03, 0.38591E+03, 0.38690E+03, 0.38788E+03, 0.38887E+03, 0.38986E+03, 0.39085E+03 &
, 0.39184E+03, 0.39283E+03, 0.39382E+03, 0.39481E+03, 0.39581E+03, 0.39680E+03, 0.39779E+03, 0.39879E+03, 0.39979E+03, 0.40079E+03 &
, 0.40179E+03, 0.40279E+03, 0.40379E+03, 0.40479E+03, 0.40579E+03, 0.40679E+03, 0.40780E+03, 0.40881E+03, 0.40981E+03, 0.41082E+03 &
, 0.41183E+03, 0.41284E+03, 0.41385E+03, 0.41486E+03, 0.41587E+03, 0.41689E+03, 0.41790E+03, 0.41892E+03, 0.41993E+03, 0.42095E+03 &
, 0.42197E+03, 0.42299E+03, 0.42401E+03, 0.42503E+03, 0.42605E+03, 0.42707E+03, 0.42810E+03, 0.42912E+03, 0.43015E+03, 0.43118E+03 &
, 0.43220E+03, 0.43323E+03, 0.43426E+03, 0.43529E+03, 0.43633E+03, 0.43736E+03, 0.43839E+03, 0.43943E+03, 0.44046E+03, 0.44150E+03 &
, 0.44254E+03, 0.44358E+03, 0.44462E+03, 0.44566E+03, 0.44670E+03, 0.44775E+03, 0.44879E+03, 0.44983E+03, 0.45088E+03, 0.45193E+03 /
 !            5           9
 DATA(QofT(          54 ,J),J=1,501)/   50.6230010986328       &
, 0.55476E+02, 0.60330E+02, 0.65184E+02, 0.70039E+02, 0.74894E+02, 0.79750E+02, 0.84606E+02, 0.89462E+02, 0.94318E+02, 0.99175E+02 &
, 0.10403E+03, 0.10889E+03, 0.11375E+03, 0.11860E+03, 0.12346E+03, 0.12832E+03, 0.13317E+03, 0.13803E+03, 0.14289E+03, 0.14775E+03 &
, 0.15261E+03, 0.15746E+03, 0.16232E+03, 0.16718E+03, 0.17204E+03, 0.17690E+03, 0.18175E+03, 0.18661E+03, 0.19147E+03, 0.19633E+03 &
, 0.20119E+03, 0.20605E+03, 0.21091E+03, 0.21576E+03, 0.22062E+03, 0.22548E+03, 0.23034E+03, 0.23520E+03, 0.24006E+03, 0.24492E+03 &
, 0.24978E+03, 0.25464E+03, 0.25950E+03, 0.26435E+03, 0.26921E+03, 0.27407E+03, 0.27893E+03, 0.28379E+03, 0.28865E+03, 0.29351E+03 &
, 0.29837E+03, 0.30323E+03, 0.30809E+03, 0.31295E+03, 0.31781E+03, 0.32267E+03, 0.32753E+03, 0.33239E+03, 0.33725E+03, 0.34211E+03 &
, 0.34697E+03, 0.35183E+03, 0.35669E+03, 0.36155E+03, 0.36641E+03, 0.37127E+03, 0.37613E+03, 0.38099E+03, 0.38585E+03, 0.39071E+03 &
, 0.39557E+03, 0.40043E+03, 0.40529E+03, 0.41016E+03, 0.41502E+03, 0.41988E+03, 0.42474E+03, 0.42960E+03, 0.43446E+03, 0.43932E+03 &
, 0.44418E+03, 0.44904E+03, 0.45390E+03, 0.45877E+03, 0.46363E+03, 0.46849E+03, 0.47335E+03, 0.47821E+03, 0.48307E+03, 0.48793E+03 &
, 0.49279E+03, 0.49766E+03, 0.50252E+03, 0.50738E+03, 0.51224E+03, 0.51710E+03, 0.52197E+03, 0.52683E+03, 0.53169E+03, 0.53655E+03 &
, 0.54141E+03, 0.54628E+03, 0.55114E+03, 0.55600E+03, 0.56086E+03, 0.56572E+03, 0.57059E+03, 0.57545E+03, 0.58031E+03, 0.58517E+03 &
, 0.59004E+03, 0.59490E+03, 0.59976E+03, 0.60463E+03, 0.60949E+03, 0.61435E+03, 0.61922E+03, 0.62408E+03, 0.62894E+03, 0.63381E+03 &
, 0.63867E+03, 0.64353E+03, 0.64840E+03, 0.65326E+03, 0.65813E+03, 0.66299E+03, 0.66785E+03, 0.67272E+03, 0.67758E+03, 0.68245E+03 &
, 0.68731E+03, 0.69218E+03, 0.69704E+03, 0.70191E+03, 0.70678E+03, 0.71164E+03, 0.71651E+03, 0.72137E+03, 0.72624E+03, 0.73111E+03 &
, 0.73598E+03, 0.74084E+03, 0.74571E+03, 0.75058E+03, 0.75545E+03, 0.76031E+03, 0.76518E+03, 0.77005E+03, 0.77492E+03, 0.77979E+03 &
, 0.78466E+03, 0.78953E+03, 0.79440E+03, 0.79927E+03, 0.80414E+03, 0.80901E+03, 0.81389E+03, 0.81876E+03, 0.82363E+03, 0.82851E+03 &
, 0.83338E+03, 0.83825E+03, 0.84313E+03, 0.84800E+03, 0.85288E+03, 0.85775E+03, 0.86263E+03, 0.86751E+03, 0.87239E+03, 0.87727E+03 &
, 0.88214E+03, 0.88702E+03, 0.89190E+03, 0.89678E+03, 0.90167E+03, 0.90655E+03, 0.91143E+03, 0.91631E+03, 0.92120E+03, 0.92608E+03 &
, 0.93097E+03, 0.93586E+03, 0.94074E+03, 0.94563E+03, 0.95052E+03, 0.95541E+03, 0.96030E+03, 0.96519E+03, 0.97009E+03, 0.97498E+03 &
, 0.97987E+03, 0.98477E+03, 0.98966E+03, 0.99456E+03, 0.99946E+03, 0.10044E+04, 0.10093E+04, 0.10142E+04, 0.10191E+04, 0.10240E+04 &
, 0.10289E+04, 0.10338E+04, 0.10387E+04, 0.10436E+04, 0.10485E+04, 0.10534E+04, 0.10583E+04, 0.10632E+04, 0.10682E+04, 0.10731E+04 &
, 0.10780E+04, 0.10829E+04, 0.10878E+04, 0.10927E+04, 0.10977E+04, 0.11026E+04, 0.11075E+04, 0.11124E+04, 0.11174E+04, 0.11223E+04 &
, 0.11272E+04, 0.11322E+04, 0.11371E+04, 0.11420E+04, 0.11470E+04, 0.11519E+04, 0.11569E+04, 0.11618E+04, 0.11668E+04, 0.11717E+04 &
, 0.11767E+04, 0.11816E+04, 0.11866E+04, 0.11915E+04, 0.11965E+04, 0.12014E+04, 0.12064E+04, 0.12113E+04, 0.12163E+04, 0.12213E+04 &
, 0.12262E+04, 0.12312E+04, 0.12362E+04, 0.12412E+04, 0.12461E+04, 0.12511E+04, 0.12561E+04, 0.12611E+04, 0.12661E+04, 0.12711E+04 &
, 0.12760E+04, 0.12810E+04, 0.12860E+04, 0.12910E+04, 0.12960E+04, 0.13010E+04, 0.13060E+04, 0.13110E+04, 0.13160E+04, 0.13211E+04 &
, 0.13261E+04, 0.13311E+04, 0.13361E+04, 0.13411E+04, 0.13462E+04, 0.13512E+04, 0.13562E+04, 0.13612E+04, 0.13663E+04, 0.13713E+04 &
, 0.13764E+04, 0.13814E+04, 0.13864E+04, 0.13915E+04, 0.13965E+04, 0.14016E+04, 0.14066E+04, 0.14117E+04, 0.14168E+04, 0.14218E+04 &
, 0.14269E+04, 0.14320E+04, 0.14370E+04, 0.14421E+04, 0.14472E+04, 0.14523E+04, 0.14574E+04, 0.14625E+04, 0.14676E+04, 0.14727E+04 &
, 0.14778E+04, 0.14829E+04, 0.14880E+04, 0.14931E+04, 0.14982E+04, 0.15033E+04, 0.15084E+04, 0.15135E+04, 0.15187E+04, 0.15238E+04 &
, 0.15289E+04, 0.15341E+04, 0.15392E+04, 0.15443E+04, 0.15495E+04, 0.15546E+04, 0.15598E+04, 0.15649E+04, 0.15701E+04, 0.15753E+04 &
, 0.15804E+04, 0.15856E+04, 0.15908E+04, 0.15959E+04, 0.16011E+04, 0.16063E+04, 0.16115E+04, 0.16167E+04, 0.16219E+04, 0.16271E+04 &
, 0.16323E+04, 0.16375E+04, 0.16427E+04, 0.16479E+04, 0.16531E+04, 0.16583E+04, 0.16636E+04, 0.16688E+04, 0.16740E+04, 0.16793E+04 &
, 0.16845E+04, 0.16897E+04, 0.16950E+04, 0.17003E+04, 0.17055E+04, 0.17108E+04, 0.17160E+04, 0.17213E+04, 0.17266E+04, 0.17318E+04 &
, 0.17371E+04, 0.17424E+04, 0.17477E+04, 0.17530E+04, 0.17583E+04, 0.17636E+04, 0.17689E+04, 0.17742E+04, 0.17795E+04, 0.17848E+04 &
, 0.17902E+04, 0.17955E+04, 0.18008E+04, 0.18062E+04, 0.18115E+04, 0.18168E+04, 0.18222E+04, 0.18275E+04, 0.18329E+04, 0.18383E+04 &
, 0.18436E+04, 0.18490E+04, 0.18544E+04, 0.18597E+04, 0.18651E+04, 0.18705E+04, 0.18759E+04, 0.18813E+04, 0.18867E+04, 0.18921E+04 &
, 0.18975E+04, 0.19029E+04, 0.19084E+04, 0.19138E+04, 0.19192E+04, 0.19246E+04, 0.19301E+04, 0.19355E+04, 0.19410E+04, 0.19464E+04 &
, 0.19519E+04, 0.19573E+04, 0.19628E+04, 0.19683E+04, 0.19738E+04, 0.19792E+04, 0.19847E+04, 0.19902E+04, 0.19957E+04, 0.20012E+04 &
, 0.20067E+04, 0.20122E+04, 0.20177E+04, 0.20233E+04, 0.20288E+04, 0.20343E+04, 0.20398E+04, 0.20454E+04, 0.20509E+04, 0.20565E+04 &
, 0.20620E+04, 0.20676E+04, 0.20731E+04, 0.20787E+04, 0.20843E+04, 0.20899E+04, 0.20954E+04, 0.21010E+04, 0.21066E+04, 0.21122E+04 &
, 0.21178E+04, 0.21234E+04, 0.21290E+04, 0.21347E+04, 0.21403E+04, 0.21459E+04, 0.21515E+04, 0.21572E+04, 0.21628E+04, 0.21685E+04 &
, 0.21741E+04, 0.21798E+04, 0.21855E+04, 0.21911E+04, 0.21968E+04, 0.22025E+04, 0.22082E+04, 0.22139E+04, 0.22195E+04, 0.22252E+04 &
, 0.22310E+04, 0.22367E+04, 0.22424E+04, 0.22481E+04, 0.22538E+04, 0.22596E+04, 0.22653E+04, 0.22710E+04, 0.22768E+04, 0.22825E+04 &
, 0.22883E+04, 0.22941E+04, 0.22998E+04, 0.23056E+04, 0.23114E+04, 0.23172E+04, 0.23230E+04, 0.23288E+04, 0.23346E+04, 0.23404E+04 &
, 0.23462E+04, 0.23520E+04, 0.23578E+04, 0.23637E+04, 0.23695E+04, 0.23753E+04, 0.23812E+04, 0.23870E+04, 0.23929E+04, 0.23987E+04 &
, 0.24046E+04, 0.24105E+04, 0.24164E+04, 0.24222E+04, 0.24281E+04, 0.24340E+04, 0.24399E+04, 0.24458E+04, 0.24518E+04, 0.24577E+04 &
, 0.24636E+04, 0.24695E+04, 0.24755E+04, 0.24814E+04, 0.24873E+04, 0.24933E+04, 0.24993E+04, 0.25052E+04, 0.25112E+04, 0.25172E+04 &
, 0.25231E+04, 0.25291E+04, 0.25351E+04, 0.25411E+04, 0.25471E+04, 0.25531E+04, 0.25591E+04, 0.25652E+04, 0.25712E+04, 0.25772E+04 &
, 0.25832E+04, 0.25893E+04, 0.25953E+04, 0.26014E+04, 0.26074E+04, 0.26135E+04, 0.26196E+04, 0.26257E+04, 0.26317E+04, 0.26378E+04 /
 !            6           1
 DATA(QofT(          55 ,J),J=1,501)/   12.7370004653931       &
, 0.14121E+02, 0.15602E+02, 0.17176E+02, 0.18841E+02, 0.20593E+02, 0.22425E+02, 0.24336E+02, 0.26320E+02, 0.28374E+02, 0.30495E+02 &
, 0.32678E+02, 0.34922E+02, 0.37223E+02, 0.39580E+02, 0.41991E+02, 0.44453E+02, 0.46965E+02, 0.49526E+02, 0.52133E+02, 0.54788E+02 &
, 0.57487E+02, 0.60230E+02, 0.63017E+02, 0.65847E+02, 0.68718E+02, 0.71630E+02, 0.74582E+02, 0.77575E+02, 0.80607E+02, 0.83678E+02 &
, 0.86786E+02, 0.89933E+02, 0.93117E+02, 0.96338E+02, 0.99595E+02, 0.10289E+03, 0.10622E+03, 0.10958E+03, 0.11298E+03, 0.11642E+03 &
, 0.11988E+03, 0.12339E+03, 0.12692E+03, 0.13049E+03, 0.13409E+03, 0.13773E+03, 0.14139E+03, 0.14509E+03, 0.14882E+03, 0.15258E+03 &
, 0.15638E+03, 0.16020E+03, 0.16406E+03, 0.16794E+03, 0.17186E+03, 0.17581E+03, 0.17978E+03, 0.18379E+03, 0.18782E+03, 0.19189E+03 &
, 0.19598E+03, 0.20010E+03, 0.20425E+03, 0.20843E+03, 0.21264E+03, 0.21688E+03, 0.22114E+03, 0.22543E+03, 0.22975E+03, 0.23410E+03 &
, 0.23848E+03, 0.24288E+03, 0.24731E+03, 0.25176E+03, 0.25624E+03, 0.26075E+03, 0.26529E+03, 0.26985E+03, 0.27444E+03, 0.27905E+03 &
, 0.28369E+03, 0.28836E+03, 0.29305E+03, 0.29777E+03, 0.30251E+03, 0.30728E+03, 0.31208E+03, 0.31690E+03, 0.32175E+03, 0.32662E+03 &
, 0.33151E+03, 0.33644E+03, 0.34138E+03, 0.34636E+03, 0.35135E+03, 0.35638E+03, 0.36142E+03, 0.36650E+03, 0.37159E+03, 0.37672E+03 &
, 0.38187E+03, 0.38704E+03, 0.39224E+03, 0.39746E+03, 0.40271E+03, 0.40798E+03, 0.41328E+03, 0.41861E+03, 0.42396E+03, 0.42933E+03 &
, 0.43473E+03, 0.44016E+03, 0.44561E+03, 0.45109E+03, 0.45659E+03, 0.46212E+03, 0.46767E+03, 0.47325E+03, 0.47886E+03, 0.48449E+03 &
, 0.49015E+03, 0.49583E+03, 0.50154E+03, 0.50728E+03, 0.51305E+03, 0.51884E+03, 0.52466E+03, 0.53050E+03, 0.53637E+03, 0.54227E+03 &
, 0.54820E+03, 0.55416E+03, 0.56014E+03, 0.56615E+03, 0.57219E+03, 0.57826E+03, 0.58435E+03, 0.59048E+03, 0.59663E+03, 0.60282E+03 &
, 0.60903E+03, 0.61527E+03, 0.62154E+03, 0.62784E+03, 0.63418E+03, 0.64054E+03, 0.64693E+03, 0.65335E+03, 0.65981E+03, 0.66630E+03 &
, 0.67281E+03, 0.67936E+03, 0.68595E+03, 0.69256E+03, 0.69921E+03, 0.70589E+03, 0.71260E+03, 0.71935E+03, 0.72612E+03, 0.73294E+03 &
, 0.73979E+03, 0.74667E+03, 0.75359E+03, 0.76054E+03, 0.76752E+03, 0.77455E+03, 0.78161E+03, 0.78870E+03, 0.79583E+03, 0.80300E+03 &
, 0.81021E+03, 0.81745E+03, 0.82473E+03, 0.83205E+03, 0.83941E+03, 0.84680E+03, 0.85424E+03, 0.86171E+03, 0.86923E+03, 0.87678E+03 &
, 0.88437E+03, 0.89201E+03, 0.89969E+03, 0.90740E+03, 0.91516E+03, 0.92296E+03, 0.93081E+03, 0.93869E+03, 0.94662E+03, 0.95460E+03 &
, 0.96261E+03, 0.97067E+03, 0.97878E+03, 0.98693E+03, 0.99513E+03, 0.10034E+04, 0.10117E+04, 0.10200E+04, 0.10284E+04, 0.10368E+04 &
, 0.10453E+04, 0.10538E+04, 0.10624E+04, 0.10710E+04, 0.10797E+04, 0.10884E+04, 0.10972E+04, 0.11060E+04, 0.11149E+04, 0.11238E+04 &
, 0.11328E+04, 0.11418E+04, 0.11509E+04, 0.11600E+04, 0.11692E+04, 0.11784E+04, 0.11877E+04, 0.11971E+04, 0.12065E+04, 0.12159E+04 &
, 0.12254E+04, 0.12350E+04, 0.12446E+04, 0.12543E+04, 0.12641E+04, 0.12739E+04, 0.12837E+04, 0.12937E+04, 0.13036E+04, 0.13137E+04 &
, 0.13238E+04, 0.13339E+04, 0.13442E+04, 0.13544E+04, 0.13648E+04, 0.13752E+04, 0.13857E+04, 0.13962E+04, 0.14068E+04, 0.14175E+04 &
, 0.14282E+04, 0.14390E+04, 0.14499E+04, 0.14608E+04, 0.14718E+04, 0.14829E+04, 0.14940E+04, 0.15052E+04, 0.15165E+04, 0.15278E+04 &
, 0.15392E+04, 0.15507E+04, 0.15623E+04, 0.15739E+04, 0.15856E+04, 0.15974E+04, 0.16092E+04, 0.16212E+04, 0.16331E+04, 0.16452E+04 &
, 0.16574E+04, 0.16696E+04, 0.16819E+04, 0.16943E+04, 0.17067E+04, 0.17193E+04, 0.17319E+04, 0.17446E+04, 0.17573E+04, 0.17702E+04 &
, 0.17831E+04, 0.17961E+04, 0.18092E+04, 0.18224E+04, 0.18357E+04, 0.18491E+04, 0.18625E+04, 0.18760E+04, 0.18896E+04, 0.19033E+04 &
, 0.19171E+04, 0.19310E+04, 0.19449E+04, 0.19590E+04, 0.19731E+04, 0.19874E+04, 0.20017E+04, 0.20161E+04, 0.20306E+04, 0.20452E+04 &
, 0.20599E+04, 0.20747E+04, 0.20896E+04, 0.21046E+04, 0.21196E+04, 0.21348E+04, 0.21501E+04, 0.21655E+04, 0.21809E+04, 0.21965E+04 &
, 0.22122E+04, 0.22279E+04, 0.22438E+04, 0.22598E+04, 0.22759E+04, 0.22921E+04, 0.23084E+04, 0.23248E+04, 0.23413E+04, 0.23579E+04 &
, 0.23746E+04, 0.23914E+04, 0.24084E+04, 0.24254E+04, 0.24426E+04, 0.24598E+04, 0.24772E+04, 0.24947E+04, 0.25123E+04, 0.25300E+04 &
, 0.25479E+04, 0.25658E+04, 0.25839E+04, 0.26021E+04, 0.26204E+04, 0.26388E+04, 0.26574E+04, 0.26760E+04, 0.26948E+04, 0.27138E+04 &
, 0.27328E+04, 0.27520E+04, 0.27712E+04, 0.27907E+04, 0.28102E+04, 0.28299E+04, 0.28497E+04, 0.28696E+04, 0.28896E+04, 0.29098E+04 &
, 0.29302E+04, 0.29506E+04, 0.29712E+04, 0.29919E+04, 0.30128E+04, 0.30338E+04, 0.30549E+04, 0.30762E+04, 0.30976E+04, 0.31191E+04 &
, 0.31408E+04, 0.31627E+04, 0.31846E+04, 0.32068E+04, 0.32290E+04, 0.32515E+04, 0.32740E+04, 0.32967E+04, 0.33196E+04, 0.33426E+04 &
, 0.33657E+04, 0.33891E+04, 0.34125E+04, 0.34361E+04, 0.34599E+04, 0.34838E+04, 0.35079E+04, 0.35322E+04, 0.35566E+04, 0.35811E+04 &
, 0.36058E+04, 0.36307E+04, 0.36558E+04, 0.36810E+04, 0.37064E+04, 0.37319E+04, 0.37576E+04, 0.37835E+04, 0.38096E+04, 0.38358E+04 &
, 0.38622E+04, 0.38887E+04, 0.39155E+04, 0.39424E+04, 0.39695E+04, 0.39968E+04, 0.40242E+04, 0.40519E+04, 0.40797E+04, 0.41077E+04 &
, 0.41358E+04, 0.41642E+04, 0.41928E+04, 0.42215E+04, 0.42504E+04, 0.42795E+04, 0.43088E+04, 0.43383E+04, 0.43680E+04, 0.43979E+04 &
, 0.44280E+04, 0.44583E+04, 0.44887E+04, 0.45194E+04, 0.45503E+04, 0.45814E+04, 0.46126E+04, 0.46441E+04, 0.46758E+04, 0.47077E+04 &
, 0.47398E+04, 0.47722E+04, 0.48047E+04, 0.48374E+04, 0.48704E+04, 0.49036E+04, 0.49370E+04, 0.49706E+04, 0.50044E+04, 0.50384E+04 &
, 0.50727E+04, 0.51072E+04, 0.51419E+04, 0.51769E+04, 0.52121E+04, 0.52475E+04, 0.52831E+04, 0.53190E+04, 0.53551E+04, 0.53914E+04 &
, 0.54280E+04, 0.54648E+04, 0.55019E+04, 0.55392E+04, 0.55767E+04, 0.56145E+04, 0.56525E+04, 0.56908E+04, 0.57293E+04, 0.57681E+04 &
, 0.58072E+04, 0.58464E+04, 0.58860E+04, 0.59258E+04, 0.59658E+04, 0.60062E+04, 0.60468E+04, 0.60876E+04, 0.61287E+04, 0.61701E+04 &
, 0.62117E+04, 0.62537E+04, 0.62959E+04, 0.63383E+04, 0.63811E+04, 0.64241E+04, 0.64674E+04, 0.65110E+04, 0.65548E+04, 0.65990E+04 &
, 0.66434E+04, 0.66881E+04, 0.67331E+04, 0.67785E+04, 0.68241E+04, 0.68699E+04, 0.69161E+04, 0.69626E+04, 0.70094E+04, 0.70565E+04 &
, 0.71039E+04, 0.71516E+04, 0.71996E+04, 0.72480E+04, 0.72966E+04, 0.73456E+04, 0.73948E+04, 0.74444E+04, 0.74943E+04, 0.75445E+04 &
, 0.75951E+04, 0.76460E+04, 0.76972E+04, 0.77487E+04, 0.78006E+04, 0.78528E+04, 0.79054E+04, 0.79582E+04, 0.80115E+04, 0.80650E+04 &
, 0.81189E+04, 0.81732E+04, 0.82278E+04, 0.82828E+04, 0.83381E+04, 0.83937E+04, 0.84498E+04, 0.85062E+04, 0.85629E+04, 0.86200E+04 /
 !            6           2
 DATA(QofT(          56 ,J),J=1,501)/   25.4740009307861       &
, 0.28242E+02, 0.31204E+02, 0.34354E+02, 0.37684E+02, 0.41186E+02, 0.44851E+02, 0.48673E+02, 0.52641E+02, 0.56749E+02, 0.60990E+02 &
, 0.65357E+02, 0.69844E+02, 0.74447E+02, 0.79161E+02, 0.83982E+02, 0.88906E+02, 0.93930E+02, 0.99051E+02, 0.10427E+03, 0.10958E+03 &
, 0.11497E+03, 0.12046E+03, 0.12603E+03, 0.13169E+03, 0.13744E+03, 0.14326E+03, 0.14917E+03, 0.15515E+03, 0.16121E+03, 0.16736E+03 &
, 0.17357E+03, 0.17987E+03, 0.18623E+03, 0.19268E+03, 0.19919E+03, 0.20578E+03, 0.21244E+03, 0.21916E+03, 0.22596E+03, 0.23283E+03 &
, 0.23977E+03, 0.24677E+03, 0.25384E+03, 0.26098E+03, 0.26818E+03, 0.27545E+03, 0.28278E+03, 0.29018E+03, 0.29764E+03, 0.30517E+03 &
, 0.31275E+03, 0.32040E+03, 0.32811E+03, 0.33589E+03, 0.34372E+03, 0.35161E+03, 0.35956E+03, 0.36757E+03, 0.37564E+03, 0.38377E+03 &
, 0.39196E+03, 0.40021E+03, 0.40851E+03, 0.41687E+03, 0.42528E+03, 0.43376E+03, 0.44229E+03, 0.45087E+03, 0.45951E+03, 0.46820E+03 &
, 0.47695E+03, 0.48575E+03, 0.49461E+03, 0.50352E+03, 0.51249E+03, 0.52151E+03, 0.53058E+03, 0.53970E+03, 0.54887E+03, 0.55810E+03 &
, 0.56739E+03, 0.57672E+03, 0.58610E+03, 0.59554E+03, 0.60502E+03, 0.61456E+03, 0.62415E+03, 0.63379E+03, 0.64349E+03, 0.65323E+03 &
, 0.66302E+03, 0.67287E+03, 0.68276E+03, 0.69271E+03, 0.70270E+03, 0.71274E+03, 0.72284E+03, 0.73299E+03, 0.74318E+03, 0.75343E+03 &
, 0.76372E+03, 0.77407E+03, 0.78446E+03, 0.79491E+03, 0.80540E+03, 0.81595E+03, 0.82655E+03, 0.83719E+03, 0.84790E+03, 0.85864E+03 &
, 0.86944E+03, 0.88029E+03, 0.89119E+03, 0.90214E+03, 0.91315E+03, 0.92420E+03, 0.93530E+03, 0.94646E+03, 0.95767E+03, 0.96893E+03 &
, 0.98024E+03, 0.99161E+03, 0.10030E+04, 0.10145E+04, 0.10260E+04, 0.10376E+04, 0.10492E+04, 0.10609E+04, 0.10727E+04, 0.10845E+04 &
, 0.10963E+04, 0.11082E+04, 0.11202E+04, 0.11322E+04, 0.11443E+04, 0.11564E+04, 0.11686E+04, 0.11808E+04, 0.11931E+04, 0.12055E+04 &
, 0.12179E+04, 0.12304E+04, 0.12429E+04, 0.12555E+04, 0.12682E+04, 0.12809E+04, 0.12937E+04, 0.13065E+04, 0.13194E+04, 0.13324E+04 &
, 0.13454E+04, 0.13585E+04, 0.13716E+04, 0.13848E+04, 0.13981E+04, 0.14115E+04, 0.14249E+04, 0.14384E+04, 0.14519E+04, 0.14655E+04 &
, 0.14792E+04, 0.14930E+04, 0.15068E+04, 0.15207E+04, 0.15346E+04, 0.15487E+04, 0.15628E+04, 0.15769E+04, 0.15912E+04, 0.16055E+04 &
, 0.16199E+04, 0.16344E+04, 0.16489E+04, 0.16635E+04, 0.16782E+04, 0.16930E+04, 0.17078E+04, 0.17228E+04, 0.17378E+04, 0.17529E+04 &
, 0.17680E+04, 0.17833E+04, 0.17986E+04, 0.18140E+04, 0.18295E+04, 0.18451E+04, 0.18607E+04, 0.18765E+04, 0.18923E+04, 0.19082E+04 &
, 0.19242E+04, 0.19403E+04, 0.19565E+04, 0.19728E+04, 0.19891E+04, 0.20056E+04, 0.20221E+04, 0.20388E+04, 0.20555E+04, 0.20723E+04 &
, 0.20893E+04, 0.21063E+04, 0.21234E+04, 0.21406E+04, 0.21579E+04, 0.21753E+04, 0.21928E+04, 0.22104E+04, 0.22281E+04, 0.22459E+04 &
, 0.22638E+04, 0.22818E+04, 0.22999E+04, 0.23181E+04, 0.23365E+04, 0.23549E+04, 0.23734E+04, 0.23921E+04, 0.24109E+04, 0.24297E+04 &
, 0.24487E+04, 0.24678E+04, 0.24870E+04, 0.25063E+04, 0.25257E+04, 0.25453E+04, 0.25650E+04, 0.25848E+04, 0.26046E+04, 0.26247E+04 &
, 0.26448E+04, 0.26651E+04, 0.26854E+04, 0.27059E+04, 0.27266E+04, 0.27473E+04, 0.27682E+04, 0.27892E+04, 0.28103E+04, 0.28316E+04 &
, 0.28530E+04, 0.28745E+04, 0.28962E+04, 0.29180E+04, 0.29399E+04, 0.29619E+04, 0.29841E+04, 0.30064E+04, 0.30289E+04, 0.30515E+04 &
, 0.30742E+04, 0.30971E+04, 0.31201E+04, 0.31433E+04, 0.31666E+04, 0.31901E+04, 0.32136E+04, 0.32374E+04, 0.32613E+04, 0.32853E+04 &
, 0.33095E+04, 0.33339E+04, 0.33583E+04, 0.33830E+04, 0.34078E+04, 0.34328E+04, 0.34579E+04, 0.34831E+04, 0.35086E+04, 0.35342E+04 &
, 0.35599E+04, 0.35858E+04, 0.36119E+04, 0.36382E+04, 0.36646E+04, 0.36911E+04, 0.37179E+04, 0.37448E+04, 0.37719E+04, 0.37991E+04 &
, 0.38266E+04, 0.38541E+04, 0.38819E+04, 0.39099E+04, 0.39380E+04, 0.39663E+04, 0.39948E+04, 0.40235E+04, 0.40524E+04, 0.40814E+04 &
, 0.41106E+04, 0.41400E+04, 0.41696E+04, 0.41994E+04, 0.42294E+04, 0.42596E+04, 0.42900E+04, 0.43205E+04, 0.43513E+04, 0.43822E+04 &
, 0.44134E+04, 0.44447E+04, 0.44763E+04, 0.45080E+04, 0.45400E+04, 0.45722E+04, 0.46046E+04, 0.46371E+04, 0.46700E+04, 0.47030E+04 &
, 0.47362E+04, 0.47696E+04, 0.48032E+04, 0.48371E+04, 0.48712E+04, 0.49055E+04, 0.49400E+04, 0.49748E+04, 0.50098E+04, 0.50449E+04 &
, 0.50804E+04, 0.51160E+04, 0.51519E+04, 0.51880E+04, 0.52244E+04, 0.52610E+04, 0.52978E+04, 0.53349E+04, 0.53721E+04, 0.54097E+04 &
, 0.54475E+04, 0.54855E+04, 0.55238E+04, 0.55623E+04, 0.56011E+04, 0.56401E+04, 0.56794E+04, 0.57189E+04, 0.57587E+04, 0.57988E+04 &
, 0.58391E+04, 0.58797E+04, 0.59205E+04, 0.59616E+04, 0.60030E+04, 0.60446E+04, 0.60866E+04, 0.61287E+04, 0.61712E+04, 0.62139E+04 &
, 0.62569E+04, 0.63002E+04, 0.63438E+04, 0.63876E+04, 0.64318E+04, 0.64762E+04, 0.65209E+04, 0.65659E+04, 0.66112E+04, 0.66568E+04 &
, 0.67027E+04, 0.67489E+04, 0.67954E+04, 0.68421E+04, 0.68892E+04, 0.69366E+04, 0.69843E+04, 0.70323E+04, 0.70807E+04, 0.71293E+04 &
, 0.71783E+04, 0.72275E+04, 0.72771E+04, 0.73270E+04, 0.73773E+04, 0.74279E+04, 0.74788E+04, 0.75300E+04, 0.75815E+04, 0.76334E+04 &
, 0.76857E+04, 0.77382E+04, 0.77911E+04, 0.78444E+04, 0.78980E+04, 0.79520E+04, 0.80062E+04, 0.80609E+04, 0.81159E+04, 0.81713E+04 &
, 0.82270E+04, 0.82831E+04, 0.83396E+04, 0.83963E+04, 0.84535E+04, 0.85111E+04, 0.85690E+04, 0.86273E+04, 0.86860E+04, 0.87451E+04 &
, 0.88045E+04, 0.88644E+04, 0.89246E+04, 0.89852E+04, 0.90462E+04, 0.91076E+04, 0.91694E+04, 0.92316E+04, 0.92941E+04, 0.93572E+04 &
, 0.94206E+04, 0.94844E+04, 0.95486E+04, 0.96133E+04, 0.96783E+04, 0.97438E+04, 0.98098E+04, 0.98761E+04, 0.99428E+04, 0.10010E+05 &
, 0.10078E+05, 0.10146E+05, 0.10214E+05, 0.10283E+05, 0.10353E+05, 0.10422E+05, 0.10493E+05, 0.10564E+05, 0.10635E+05, 0.10706E+05 &
, 0.10778E+05, 0.10851E+05, 0.10924E+05, 0.10998E+05, 0.11072E+05, 0.11146E+05, 0.11221E+05, 0.11297E+05, 0.11373E+05, 0.11449E+05 &
, 0.11526E+05, 0.11603E+05, 0.11681E+05, 0.11760E+05, 0.11838E+05, 0.11918E+05, 0.11998E+05, 0.12078E+05, 0.12159E+05, 0.12241E+05 &
, 0.12323E+05, 0.12405E+05, 0.12488E+05, 0.12572E+05, 0.12656E+05, 0.12741E+05, 0.12826E+05, 0.12912E+05, 0.12998E+05, 0.13085E+05 &
, 0.13172E+05, 0.13260E+05, 0.13349E+05, 0.13438E+05, 0.13527E+05, 0.13618E+05, 0.13708E+05, 0.13800E+05, 0.13892E+05, 0.13984E+05 &
, 0.14078E+05, 0.14171E+05, 0.14266E+05, 0.14361E+05, 0.14456E+05, 0.14552E+05, 0.14649E+05, 0.14747E+05, 0.14845E+05, 0.14943E+05 &
, 0.15043E+05, 0.15142E+05, 0.15243E+05, 0.15344E+05, 0.15446E+05, 0.15549E+05, 0.15652E+05, 0.15755E+05, 0.15860E+05, 0.15965E+05 &
, 0.16071E+05, 0.16177E+05, 0.16284E+05, 0.16392E+05, 0.16501E+05, 0.16610E+05, 0.16720E+05, 0.16830E+05, 0.16942E+05, 0.17054E+05 /
 !            6           3
 DATA(QofT(          57 ,J),J=1,501)/   91.9010009765625       &
, 0.10451E+03, 0.11774E+03, 0.13156E+03, 0.14596E+03, 0.16091E+03, 0.17639E+03, 0.19239E+03, 0.20889E+03, 0.22587E+03, 0.24332E+03 &
, 0.26122E+03, 0.27956E+03, 0.29832E+03, 0.31750E+03, 0.33709E+03, 0.35707E+03, 0.37744E+03, 0.39818E+03, 0.41930E+03, 0.44079E+03 &
, 0.46263E+03, 0.48483E+03, 0.50737E+03, 0.53026E+03, 0.55348E+03, 0.57703E+03, 0.60090E+03, 0.62510E+03, 0.64962E+03, 0.67445E+03 &
, 0.69959E+03, 0.72503E+03, 0.75078E+03, 0.77683E+03, 0.80317E+03, 0.82980E+03, 0.85672E+03, 0.88393E+03, 0.91142E+03, 0.93920E+03 &
, 0.96725E+03, 0.99557E+03, 0.10242E+04, 0.10530E+04, 0.10822E+04, 0.11116E+04, 0.11412E+04, 0.11711E+04, 0.12013E+04, 0.12318E+04 &
, 0.12624E+04, 0.12934E+04, 0.13246E+04, 0.13560E+04, 0.13877E+04, 0.14196E+04, 0.14518E+04, 0.14842E+04, 0.15168E+04, 0.15497E+04 &
, 0.15829E+04, 0.16162E+04, 0.16498E+04, 0.16836E+04, 0.17177E+04, 0.17520E+04, 0.17865E+04, 0.18212E+04, 0.18562E+04, 0.18914E+04 &
, 0.19268E+04, 0.19624E+04, 0.19982E+04, 0.20343E+04, 0.20706E+04, 0.21071E+04, 0.21438E+04, 0.21808E+04, 0.22179E+04, 0.22553E+04 &
, 0.22929E+04, 0.23307E+04, 0.23687E+04, 0.24070E+04, 0.24454E+04, 0.24840E+04, 0.25229E+04, 0.25620E+04, 0.26013E+04, 0.26408E+04 &
, 0.26805E+04, 0.27204E+04, 0.27605E+04, 0.28009E+04, 0.28415E+04, 0.28822E+04, 0.29232E+04, 0.29644E+04, 0.30058E+04, 0.30474E+04 &
, 0.30892E+04, 0.31313E+04, 0.31735E+04, 0.32160E+04, 0.32587E+04, 0.33016E+04, 0.33447E+04, 0.33881E+04, 0.34316E+04, 0.34754E+04 &
, 0.35194E+04, 0.35636E+04, 0.36080E+04, 0.36527E+04, 0.36976E+04, 0.37427E+04, 0.37880E+04, 0.38336E+04, 0.38793E+04, 0.39254E+04 &
, 0.39716E+04, 0.40181E+04, 0.40648E+04, 0.41117E+04, 0.41589E+04, 0.42064E+04, 0.42540E+04, 0.43019E+04, 0.43501E+04, 0.43985E+04 &
, 0.44471E+04, 0.44960E+04, 0.45452E+04, 0.45946E+04, 0.46442E+04, 0.46941E+04, 0.47443E+04, 0.47947E+04, 0.48454E+04, 0.48964E+04 &
, 0.49476E+04, 0.49991E+04, 0.50509E+04, 0.51029E+04, 0.51552E+04, 0.52078E+04, 0.52607E+04, 0.53138E+04, 0.53673E+04, 0.54210E+04 &
, 0.54750E+04, 0.55294E+04, 0.55840E+04, 0.56389E+04, 0.56941E+04, 0.57496E+04, 0.58054E+04, 0.58616E+04, 0.59180E+04, 0.59748E+04 &
, 0.60319E+04, 0.60893E+04, 0.61470E+04, 0.62050E+04, 0.62634E+04, 0.63221E+04, 0.63811E+04, 0.64405E+04, 0.65002E+04, 0.65603E+04 &
, 0.66207E+04, 0.66815E+04, 0.67426E+04, 0.68041E+04, 0.68659E+04, 0.69281E+04, 0.69907E+04, 0.70536E+04, 0.71169E+04, 0.71806E+04 &
, 0.72447E+04, 0.73091E+04, 0.73739E+04, 0.74391E+04, 0.75048E+04, 0.75708E+04, 0.76372E+04, 0.77040E+04, 0.77712E+04, 0.78388E+04 &
, 0.79069E+04, 0.79753E+04, 0.80442E+04, 0.81135E+04, 0.81833E+04, 0.82534E+04, 0.83240E+04, 0.83951E+04, 0.84666E+04, 0.85385E+04 &
, 0.86109E+04, 0.86838E+04, 0.87571E+04, 0.88308E+04, 0.89051E+04, 0.89798E+04, 0.90549E+04, 0.91306E+04, 0.92068E+04, 0.92834E+04 &
, 0.93605E+04, 0.94381E+04, 0.95162E+04, 0.95948E+04, 0.96740E+04, 0.97536E+04, 0.98338E+04, 0.99144E+04, 0.99956E+04, 0.10077E+05 &
, 0.10160E+05, 0.10242E+05, 0.10326E+05, 0.10410E+05, 0.10494E+05, 0.10579E+05, 0.10665E+05, 0.10751E+05, 0.10837E+05, 0.10925E+05 &
, 0.11013E+05, 0.11101E+05, 0.11190E+05, 0.11280E+05, 0.11370E+05, 0.11461E+05, 0.11552E+05, 0.11644E+05, 0.11737E+05, 0.11830E+05 &
, 0.11924E+05, 0.12018E+05, 0.12113E+05, 0.12209E+05, 0.12306E+05, 0.12403E+05, 0.12500E+05, 0.12599E+05, 0.12698E+05, 0.12798E+05 &
, 0.12898E+05, 0.12999E+05, 0.13101E+05, 0.13203E+05, 0.13307E+05, 0.13411E+05, 0.13515E+05, 0.13620E+05, 0.13726E+05, 0.13833E+05 &
, 0.13941E+05, 0.14049E+05, 0.14158E+05, 0.14267E+05, 0.14378E+05, 0.14489E+05, 0.14601E+05, 0.14714E+05, 0.14827E+05, 0.14942E+05 &
, 0.15057E+05, 0.15173E+05, 0.15289E+05, 0.15407E+05, 0.15525E+05, 0.15644E+05, 0.15764E+05, 0.15885E+05, 0.16006E+05, 0.16129E+05 &
, 0.16252E+05, 0.16376E+05, 0.16501E+05, 0.16627E+05, 0.16754E+05, 0.16881E+05, 0.17010E+05, 0.17139E+05, 0.17270E+05, 0.17401E+05 &
, 0.17533E+05, 0.17666E+05, 0.17800E+05, 0.17935E+05, 0.18071E+05, 0.18207E+05, 0.18345E+05, 0.18484E+05, 0.18623E+05, 0.18764E+05 &
, 0.18906E+05, 0.19048E+05, 0.19192E+05, 0.19336E+05, 0.19482E+05, 0.19628E+05, 0.19776E+05, 0.19925E+05, 0.20074E+05, 0.20225E+05 &
, 0.20377E+05, 0.20530E+05, 0.20684E+05, 0.20839E+05, 0.20995E+05, 0.21152E+05, 0.21310E+05, 0.21470E+05, 0.21630E+05, 0.21792E+05 &
, 0.21954E+05, 0.22118E+05, 0.22283E+05, 0.22450E+05, 0.22617E+05, 0.22785E+05, 0.22955E+05, 0.23126E+05, 0.23298E+05, 0.23471E+05 &
, 0.23646E+05, 0.23822E+05, 0.23999E+05, 0.24177E+05, 0.24356E+05, 0.24537E+05, 0.24719E+05, 0.24902E+05, 0.25087E+05, 0.25272E+05 &
, 0.25459E+05, 0.25648E+05, 0.25838E+05, 0.26029E+05, 0.26221E+05, 0.26415E+05, 0.26610E+05, 0.26806E+05, 0.27004E+05, 0.27203E+05 &
, 0.27404E+05, 0.27606E+05, 0.27809E+05, 0.28014E+05, 0.28221E+05, 0.28428E+05, 0.28637E+05, 0.28848E+05, 0.29060E+05, 0.29274E+05 &
, 0.29489E+05, 0.29706E+05, 0.29924E+05, 0.30143E+05, 0.30364E+05, 0.30587E+05, 0.30811E+05, 0.31037E+05, 0.31265E+05, 0.31494E+05 &
, 0.31724E+05, 0.31956E+05, 0.32190E+05, 0.32426E+05, 0.32663E+05, 0.32902E+05, 0.33142E+05, 0.33384E+05, 0.33628E+05, 0.33873E+05 &
, 0.34120E+05, 0.34369E+05, 0.34620E+05, 0.34872E+05, 0.35127E+05, 0.35382E+05, 0.35640E+05, 0.35900E+05, 0.36161E+05, 0.36424E+05 &
, 0.36689E+05, 0.36956E+05, 0.37224E+05, 0.37495E+05, 0.37767E+05, 0.38042E+05, 0.38318E+05, 0.38596E+05, 0.38876E+05, 0.39158E+05 &
, 0.39442E+05, 0.39728E+05, 0.40016E+05, 0.40305E+05, 0.40597E+05, 0.40891E+05, 0.41187E+05, 0.41485E+05, 0.41785E+05, 0.42087E+05 &
, 0.42392E+05, 0.42698E+05, 0.43006E+05, 0.43317E+05, 0.43630E+05, 0.43945E+05, 0.44262E+05, 0.44581E+05, 0.44902E+05, 0.45226E+05 &
, 0.45552E+05, 0.45880E+05, 0.46211E+05, 0.46543E+05, 0.46878E+05, 0.47216E+05, 0.47555E+05, 0.47897E+05, 0.48241E+05, 0.48588E+05 &
, 0.48937E+05, 0.49289E+05, 0.49643E+05, 0.49999E+05, 0.50358E+05, 0.50719E+05, 0.51083E+05, 0.51449E+05, 0.51818E+05, 0.52189E+05 &
, 0.52563E+05, 0.52939E+05, 0.53318E+05, 0.53700E+05, 0.54084E+05, 0.54471E+05, 0.54860E+05, 0.55252E+05, 0.55647E+05, 0.56045E+05 &
, 0.56445E+05, 0.56848E+05, 0.57254E+05, 0.57662E+05, 0.58073E+05, 0.58487E+05, 0.58904E+05, 0.59324E+05, 0.59747E+05, 0.60172E+05 &
, 0.60601E+05, 0.61032E+05, 0.61466E+05, 0.61903E+05, 0.62344E+05, 0.62787E+05, 0.63233E+05, 0.63682E+05, 0.64135E+05, 0.64590E+05 &
, 0.65048E+05, 0.65510E+05, 0.65975E+05, 0.66443E+05, 0.66914E+05, 0.67388E+05, 0.67865E+05, 0.68346E+05, 0.68830E+05, 0.69317E+05 &
, 0.69808E+05, 0.70302E+05, 0.70799E+05, 0.71299E+05, 0.71803E+05, 0.72310E+05, 0.72821E+05, 0.73335E+05, 0.73853E+05, 0.74374E+05 &
, 0.74899E+05, 0.75427E+05, 0.75959E+05, 0.76494E+05, 0.77033E+05, 0.77576E+05, 0.78122E+05, 0.78672E+05, 0.79225E+05, 0.79783E+05 /
 !            6           4
 DATA(QofT(          58 ,J),J=1,501)/   183.949996948242       &
, 0.20920E+03, 0.23567E+03, 0.26334E+03, 0.29215E+03, 0.32208E+03, 0.35307E+03, 0.38510E+03, 0.41813E+03, 0.45212E+03, 0.48704E+03 &
, 0.52287E+03, 0.55957E+03, 0.59713E+03, 0.63553E+03, 0.67473E+03, 0.71473E+03, 0.75550E+03, 0.79703E+03, 0.83930E+03, 0.88231E+03 &
, 0.92604E+03, 0.97047E+03, 0.10156E+04, 0.10614E+04, 0.11079E+04, 0.11550E+04, 0.12028E+04, 0.12512E+04, 0.13003E+04, 0.13500E+04 &
, 0.14003E+04, 0.14513E+04, 0.15028E+04, 0.15550E+04, 0.16077E+04, 0.16610E+04, 0.17149E+04, 0.17693E+04, 0.18244E+04, 0.18800E+04 &
, 0.19361E+04, 0.19928E+04, 0.20500E+04, 0.21078E+04, 0.21661E+04, 0.22250E+04, 0.22844E+04, 0.23442E+04, 0.24047E+04, 0.24656E+04 &
, 0.25270E+04, 0.25889E+04, 0.26514E+04, 0.27143E+04, 0.27777E+04, 0.28416E+04, 0.29060E+04, 0.29709E+04, 0.30362E+04, 0.31021E+04 &
, 0.31684E+04, 0.32352E+04, 0.33024E+04, 0.33701E+04, 0.34383E+04, 0.35069E+04, 0.35760E+04, 0.36455E+04, 0.37155E+04, 0.37859E+04 &
, 0.38568E+04, 0.39281E+04, 0.39999E+04, 0.40721E+04, 0.41447E+04, 0.42178E+04, 0.42913E+04, 0.43653E+04, 0.44396E+04, 0.45145E+04 &
, 0.45897E+04, 0.46654E+04, 0.47415E+04, 0.48180E+04, 0.48950E+04, 0.49723E+04, 0.50502E+04, 0.51284E+04, 0.52070E+04, 0.52861E+04 &
, 0.53656E+04, 0.54456E+04, 0.55259E+04, 0.56067E+04, 0.56879E+04, 0.57695E+04, 0.58515E+04, 0.59340E+04, 0.60169E+04, 0.61002E+04 &
, 0.61839E+04, 0.62681E+04, 0.63527E+04, 0.64378E+04, 0.65232E+04, 0.66091E+04, 0.66954E+04, 0.67822E+04, 0.68694E+04, 0.69570E+04 &
, 0.70451E+04, 0.71337E+04, 0.72227E+04, 0.73120E+04, 0.74019E+04, 0.74922E+04, 0.75830E+04, 0.76742E+04, 0.77659E+04, 0.78581E+04 &
, 0.79507E+04, 0.80437E+04, 0.81373E+04, 0.82313E+04, 0.83258E+04, 0.84208E+04, 0.85163E+04, 0.86122E+04, 0.87087E+04, 0.88055E+04 &
, 0.89030E+04, 0.90009E+04, 0.90993E+04, 0.91982E+04, 0.92977E+04, 0.93977E+04, 0.94981E+04, 0.95992E+04, 0.97007E+04, 0.98028E+04 &
, 0.99053E+04, 0.10009E+05, 0.10112E+05, 0.10216E+05, 0.10321E+05, 0.10427E+05, 0.10533E+05, 0.10639E+05, 0.10746E+05, 0.10854E+05 &
, 0.10962E+05, 0.11071E+05, 0.11180E+05, 0.11290E+05, 0.11401E+05, 0.11512E+05, 0.11624E+05, 0.11736E+05, 0.11849E+05, 0.11963E+05 &
, 0.12077E+05, 0.12192E+05, 0.12308E+05, 0.12424E+05, 0.12541E+05, 0.12659E+05, 0.12777E+05, 0.12896E+05, 0.13016E+05, 0.13136E+05 &
, 0.13257E+05, 0.13379E+05, 0.13502E+05, 0.13625E+05, 0.13749E+05, 0.13873E+05, 0.13999E+05, 0.14125E+05, 0.14252E+05, 0.14379E+05 &
, 0.14508E+05, 0.14637E+05, 0.14767E+05, 0.14898E+05, 0.15029E+05, 0.15161E+05, 0.15294E+05, 0.15428E+05, 0.15563E+05, 0.15699E+05 &
, 0.15835E+05, 0.15972E+05, 0.16110E+05, 0.16249E+05, 0.16389E+05, 0.16530E+05, 0.16671E+05, 0.16814E+05, 0.16957E+05, 0.17101E+05 &
, 0.17246E+05, 0.17392E+05, 0.17539E+05, 0.17687E+05, 0.17836E+05, 0.17986E+05, 0.18136E+05, 0.18288E+05, 0.18441E+05, 0.18594E+05 &
, 0.18749E+05, 0.18905E+05, 0.19061E+05, 0.19219E+05, 0.19377E+05, 0.19537E+05, 0.19698E+05, 0.19859E+05, 0.20022E+05, 0.20186E+05 &
, 0.20351E+05, 0.20517E+05, 0.20684E+05, 0.20852E+05, 0.21022E+05, 0.21192E+05, 0.21364E+05, 0.21536E+05, 0.21710E+05, 0.21885E+05 &
, 0.22061E+05, 0.22238E+05, 0.22417E+05, 0.22597E+05, 0.22777E+05, 0.22960E+05, 0.23143E+05, 0.23327E+05, 0.23513E+05, 0.23700E+05 &
, 0.23888E+05, 0.24078E+05, 0.24269E+05, 0.24461E+05, 0.24654E+05, 0.24849E+05, 0.25045E+05, 0.25242E+05, 0.25441E+05, 0.25641E+05 &
, 0.25842E+05, 0.26045E+05, 0.26249E+05, 0.26455E+05, 0.26661E+05, 0.26870E+05, 0.27080E+05, 0.27291E+05, 0.27503E+05, 0.27717E+05 &
, 0.27933E+05, 0.28150E+05, 0.28368E+05, 0.28588E+05, 0.28810E+05, 0.29033E+05, 0.29258E+05, 0.29484E+05, 0.29711E+05, 0.29941E+05 &
, 0.30171E+05, 0.30404E+05, 0.30638E+05, 0.30873E+05, 0.31111E+05, 0.31350E+05, 0.31590E+05, 0.31832E+05, 0.32076E+05, 0.32322E+05 &
, 0.32569E+05, 0.32818E+05, 0.33069E+05, 0.33321E+05, 0.33575E+05, 0.33831E+05, 0.34089E+05, 0.34349E+05, 0.34610E+05, 0.34873E+05 &
, 0.35138E+05, 0.35405E+05, 0.35674E+05, 0.35944E+05, 0.36216E+05, 0.36491E+05, 0.36767E+05, 0.37045E+05, 0.37325E+05, 0.37607E+05 &
, 0.37891E+05, 0.38177E+05, 0.38465E+05, 0.38755E+05, 0.39048E+05, 0.39342E+05, 0.39638E+05, 0.39936E+05, 0.40236E+05, 0.40539E+05 &
, 0.40843E+05, 0.41149E+05, 0.41458E+05, 0.41769E+05, 0.42083E+05, 0.42398E+05, 0.42715E+05, 0.43035E+05, 0.43357E+05, 0.43681E+05 &
, 0.44007E+05, 0.44336E+05, 0.44667E+05, 0.45001E+05, 0.45336E+05, 0.45674E+05, 0.46015E+05, 0.46358E+05, 0.46703E+05, 0.47050E+05 &
, 0.47401E+05, 0.47753E+05, 0.48108E+05, 0.48465E+05, 0.48825E+05, 0.49188E+05, 0.49553E+05, 0.49921E+05, 0.50291E+05, 0.50663E+05 &
, 0.51039E+05, 0.51417E+05, 0.51797E+05, 0.52181E+05, 0.52567E+05, 0.52955E+05, 0.53347E+05, 0.53741E+05, 0.54137E+05, 0.54537E+05 &
, 0.54940E+05, 0.55345E+05, 0.55753E+05, 0.56164E+05, 0.56578E+05, 0.56995E+05, 0.57414E+05, 0.57837E+05, 0.58262E+05, 0.58691E+05 &
, 0.59123E+05, 0.59557E+05, 0.59995E+05, 0.60435E+05, 0.60879E+05, 0.61326E+05, 0.61776E+05, 0.62228E+05, 0.62685E+05, 0.63144E+05 &
, 0.63607E+05, 0.64073E+05, 0.64542E+05, 0.65014E+05, 0.65490E+05, 0.65968E+05, 0.66451E+05, 0.66937E+05, 0.67426E+05, 0.67918E+05 &
, 0.68414E+05, 0.68914E+05, 0.69416E+05, 0.69923E+05, 0.70432E+05, 0.70946E+05, 0.71463E+05, 0.71984E+05, 0.72508E+05, 0.73035E+05 &
, 0.73567E+05, 0.74102E+05, 0.74641E+05, 0.75184E+05, 0.75730E+05, 0.76280E+05, 0.76835E+05, 0.77392E+05, 0.77954E+05, 0.78520E+05 &
, 0.79089E+05, 0.79663E+05, 0.80240E+05, 0.80822E+05, 0.81408E+05, 0.81997E+05, 0.82591E+05, 0.83189E+05, 0.83791E+05, 0.84397E+05 &
, 0.85007E+05, 0.85622E+05, 0.86241E+05, 0.86864E+05, 0.87491E+05, 0.88123E+05, 0.88759E+05, 0.89399E+05, 0.90044E+05, 0.90693E+05 &
, 0.91347E+05, 0.92006E+05, 0.92669E+05, 0.93336E+05, 0.94008E+05, 0.94685E+05, 0.95366E+05, 0.96052E+05, 0.96743E+05, 0.97438E+05 &
, 0.98139E+05, 0.98844E+05, 0.99554E+05, 0.10027E+06, 0.10099E+06, 0.10171E+06, 0.10244E+06, 0.10318E+06, 0.10392E+06, 0.10466E+06 &
, 0.10541E+06, 0.10617E+06, 0.10693E+06, 0.10769E+06, 0.10846E+06, 0.10924E+06, 0.11002E+06, 0.11081E+06, 0.11160E+06, 0.11240E+06 &
, 0.11320E+06, 0.11401E+06, 0.11482E+06, 0.11564E+06, 0.11647E+06, 0.11730E+06, 0.11813E+06, 0.11898E+06, 0.11982E+06, 0.12068E+06 &
, 0.12154E+06, 0.12240E+06, 0.12327E+06, 0.12415E+06, 0.12503E+06, 0.12592E+06, 0.12682E+06, 0.12772E+06, 0.12863E+06, 0.12954E+06 &
, 0.13046E+06, 0.13139E+06, 0.13232E+06, 0.13326E+06, 0.13420E+06, 0.13515E+06, 0.13611E+06, 0.13707E+06, 0.13804E+06, 0.13902E+06 &
, 0.14001E+06, 0.14100E+06, 0.14199E+06, 0.14300E+06, 0.14401E+06, 0.14503E+06, 0.14605E+06, 0.14708E+06, 0.14812E+06, 0.14917E+06 &
, 0.15022E+06, 0.15128E+06, 0.15235E+06, 0.15342E+06, 0.15450E+06, 0.15559E+06, 0.15668E+06, 0.15779E+06, 0.15890E+06, 0.16002E+06 /
 !            7           1
 DATA(QofT(          59 ,J),J=1,501)/   15.4119997024536       &
, 0.16849E+02, 0.18288E+02, 0.19729E+02, 0.21171E+02, 0.22615E+02, 0.24059E+02, 0.25504E+02, 0.26950E+02, 0.28396E+02, 0.29843E+02 &
, 0.31290E+02, 0.32737E+02, 0.34185E+02, 0.35633E+02, 0.37082E+02, 0.38530E+02, 0.39979E+02, 0.41428E+02, 0.42877E+02, 0.44326E+02 &
, 0.45775E+02, 0.47224E+02, 0.48674E+02, 0.50123E+02, 0.51573E+02, 0.53023E+02, 0.54473E+02, 0.55923E+02, 0.57373E+02, 0.58823E+02 &
, 0.60273E+02, 0.61723E+02, 0.63174E+02, 0.64624E+02, 0.66074E+02, 0.67525E+02, 0.68975E+02, 0.70426E+02, 0.71877E+02, 0.73327E+02 &
, 0.74778E+02, 0.76229E+02, 0.77679E+02, 0.79130E+02, 0.80581E+02, 0.82032E+02, 0.83483E+02, 0.84934E+02, 0.86385E+02, 0.87836E+02 &
, 0.89287E+02, 0.90738E+02, 0.92189E+02, 0.93640E+02, 0.95091E+02, 0.96543E+02, 0.97994E+02, 0.99445E+02, 0.10090E+03, 0.10235E+03 &
, 0.10380E+03, 0.10525E+03, 0.10670E+03, 0.10815E+03, 0.10960E+03, 0.11106E+03, 0.11251E+03, 0.11396E+03, 0.11541E+03, 0.11686E+03 &
, 0.11831E+03, 0.11977E+03, 0.12122E+03, 0.12267E+03, 0.12412E+03, 0.12557E+03, 0.12702E+03, 0.12848E+03, 0.12993E+03, 0.13138E+03 &
, 0.13283E+03, 0.13428E+03, 0.13574E+03, 0.13719E+03, 0.13864E+03, 0.14009E+03, 0.14154E+03, 0.14300E+03, 0.14445E+03, 0.14590E+03 &
, 0.14735E+03, 0.14881E+03, 0.15026E+03, 0.15171E+03, 0.15316E+03, 0.15462E+03, 0.15607E+03, 0.15752E+03, 0.15897E+03, 0.16043E+03 &
, 0.16188E+03, 0.16333E+03, 0.16479E+03, 0.16624E+03, 0.16769E+03, 0.16915E+03, 0.17060E+03, 0.17205E+03, 0.17351E+03, 0.17496E+03 &
, 0.17641E+03, 0.17787E+03, 0.17932E+03, 0.18078E+03, 0.18223E+03, 0.18369E+03, 0.18514E+03, 0.18660E+03, 0.18805E+03, 0.18951E+03 &
, 0.19096E+03, 0.19242E+03, 0.19387E+03, 0.19533E+03, 0.19678E+03, 0.19824E+03, 0.19970E+03, 0.20115E+03, 0.20261E+03, 0.20407E+03 &
, 0.20553E+03, 0.20698E+03, 0.20844E+03, 0.20990E+03, 0.21136E+03, 0.21282E+03, 0.21428E+03, 0.21573E+03, 0.21719E+03, 0.21865E+03 &
, 0.22011E+03, 0.22157E+03, 0.22304E+03, 0.22450E+03, 0.22596E+03, 0.22742E+03, 0.22888E+03, 0.23035E+03, 0.23181E+03, 0.23327E+03 &
, 0.23474E+03, 0.23620E+03, 0.23767E+03, 0.23913E+03, 0.24060E+03, 0.24206E+03, 0.24353E+03, 0.24500E+03, 0.24646E+03, 0.24793E+03 &
, 0.24940E+03, 0.25087E+03, 0.25234E+03, 0.25381E+03, 0.25528E+03, 0.25675E+03, 0.25822E+03, 0.25970E+03, 0.26117E+03, 0.26264E+03 &
, 0.26412E+03, 0.26559E+03, 0.26707E+03, 0.26854E+03, 0.27002E+03, 0.27150E+03, 0.27298E+03, 0.27445E+03, 0.27593E+03, 0.27741E+03 &
, 0.27889E+03, 0.28038E+03, 0.28186E+03, 0.28334E+03, 0.28482E+03, 0.28631E+03, 0.28779E+03, 0.28928E+03, 0.29077E+03, 0.29225E+03 &
, 0.29374E+03, 0.29523E+03, 0.29672E+03, 0.29821E+03, 0.29970E+03, 0.30120E+03, 0.30269E+03, 0.30418E+03, 0.30568E+03, 0.30717E+03 &
, 0.30867E+03, 0.31017E+03, 0.31167E+03, 0.31317E+03, 0.31467E+03, 0.31617E+03, 0.31767E+03, 0.31917E+03, 0.32068E+03, 0.32218E+03 &
, 0.32369E+03, 0.32519E+03, 0.32670E+03, 0.32821E+03, 0.32972E+03, 0.33123E+03, 0.33274E+03, 0.33426E+03, 0.33577E+03, 0.33729E+03 &
, 0.33880E+03, 0.34032E+03, 0.34184E+03, 0.34335E+03, 0.34487E+03, 0.34640E+03, 0.34792E+03, 0.34944E+03, 0.35097E+03, 0.35249E+03 &
, 0.35402E+03, 0.35555E+03, 0.35707E+03, 0.35860E+03, 0.36013E+03, 0.36167E+03, 0.36320E+03, 0.36473E+03, 0.36627E+03, 0.36781E+03 &
, 0.36934E+03, 0.37088E+03, 0.37242E+03, 0.37396E+03, 0.37551E+03, 0.37705E+03, 0.37859E+03, 0.38014E+03, 0.38169E+03, 0.38324E+03 &
, 0.38479E+03, 0.38634E+03, 0.38789E+03, 0.38944E+03, 0.39099E+03, 0.39255E+03, 0.39411E+03, 0.39566E+03, 0.39722E+03, 0.39878E+03 &
, 0.40035E+03, 0.40191E+03, 0.40347E+03, 0.40504E+03, 0.40660E+03, 0.40817E+03, 0.40974E+03, 0.41131E+03, 0.41288E+03, 0.41446E+03 &
, 0.41603E+03, 0.41760E+03, 0.41918E+03, 0.42076E+03, 0.42234E+03, 0.42392E+03, 0.42550E+03, 0.42708E+03, 0.42867E+03, 0.43025E+03 &
, 0.43184E+03, 0.43343E+03, 0.43502E+03, 0.43661E+03, 0.43820E+03, 0.43979E+03, 0.44138E+03, 0.44298E+03, 0.44458E+03, 0.44618E+03 &
, 0.44778E+03, 0.44938E+03, 0.45098E+03, 0.45258E+03, 0.45419E+03, 0.45579E+03, 0.45740E+03, 0.45901E+03, 0.46062E+03, 0.46223E+03 &
, 0.46384E+03, 0.46546E+03, 0.46707E+03, 0.46869E+03, 0.47030E+03, 0.47192E+03, 0.47354E+03, 0.47517E+03, 0.47679E+03, 0.47841E+03 &
, 0.48004E+03, 0.48167E+03, 0.48329E+03, 0.48492E+03, 0.48655E+03, 0.48819E+03, 0.48982E+03, 0.49145E+03, 0.49309E+03, 0.49473E+03 &
, 0.49637E+03, 0.49801E+03, 0.49965E+03, 0.50129E+03, 0.50293E+03, 0.50458E+03, 0.50623E+03, 0.50787E+03, 0.50952E+03, 0.51117E+03 &
, 0.51282E+03, 0.51448E+03, 0.51613E+03, 0.51779E+03, 0.51944E+03, 0.52110E+03, 0.52276E+03, 0.52442E+03, 0.52609E+03, 0.52775E+03 &
, 0.52941E+03, 0.53108E+03, 0.53275E+03, 0.53442E+03, 0.53609E+03, 0.53776E+03, 0.53943E+03, 0.54110E+03, 0.54278E+03, 0.54445E+03 &
, 0.54613E+03, 0.54781E+03, 0.54949E+03, 0.55117E+03, 0.55285E+03, 0.55454E+03, 0.55622E+03, 0.55791E+03, 0.55960E+03, 0.56129E+03 &
, 0.56298E+03, 0.56467E+03, 0.56636E+03, 0.56806E+03, 0.56975E+03, 0.57145E+03, 0.57315E+03, 0.57485E+03, 0.57655E+03, 0.57825E+03 &
, 0.57995E+03, 0.58166E+03, 0.58336E+03, 0.58507E+03, 0.58678E+03, 0.58849E+03, 0.59020E+03, 0.59191E+03, 0.59362E+03, 0.59533E+03 &
, 0.59705E+03, 0.59877E+03, 0.60049E+03, 0.60220E+03, 0.60392E+03, 0.60565E+03, 0.60737E+03, 0.60909E+03, 0.61082E+03, 0.61255E+03 &
, 0.61427E+03, 0.61600E+03, 0.61773E+03, 0.61946E+03, 0.62120E+03, 0.62293E+03, 0.62467E+03, 0.62640E+03, 0.62814E+03, 0.62988E+03 &
, 0.63162E+03, 0.63336E+03, 0.63510E+03, 0.63684E+03, 0.63859E+03, 0.64033E+03, 0.64208E+03, 0.64383E+03, 0.64558E+03, 0.64733E+03 &
, 0.64908E+03, 0.65083E+03, 0.65259E+03, 0.65434E+03, 0.65610E+03, 0.65786E+03, 0.65962E+03, 0.66138E+03, 0.66314E+03, 0.66490E+03 &
, 0.66666E+03, 0.66843E+03, 0.67019E+03, 0.67196E+03, 0.67373E+03, 0.67550E+03, 0.67727E+03, 0.67904E+03, 0.68081E+03, 0.68259E+03 &
, 0.68436E+03, 0.68614E+03, 0.68791E+03, 0.68969E+03, 0.69147E+03, 0.69325E+03, 0.69503E+03, 0.69682E+03, 0.69860E+03, 0.70039E+03 &
, 0.70217E+03, 0.70396E+03, 0.70575E+03, 0.70754E+03, 0.70933E+03, 0.71112E+03, 0.71291E+03, 0.71470E+03, 0.71650E+03, 0.71830E+03 &
, 0.72009E+03, 0.72189E+03, 0.72369E+03, 0.72549E+03, 0.72729E+03, 0.72909E+03, 0.73090E+03, 0.73270E+03, 0.73451E+03, 0.73631E+03 &
, 0.73812E+03, 0.73993E+03, 0.74174E+03, 0.74355E+03, 0.74536E+03, 0.74718E+03, 0.74899E+03, 0.75081E+03, 0.75262E+03, 0.75444E+03 &
, 0.75626E+03, 0.75808E+03, 0.75990E+03, 0.76172E+03, 0.76354E+03, 0.76537E+03, 0.76719E+03, 0.76902E+03, 0.77084E+03, 0.77267E+03 &
, 0.77450E+03, 0.77633E+03, 0.77816E+03, 0.77999E+03, 0.78182E+03, 0.78366E+03, 0.78549E+03, 0.78733E+03, 0.78917E+03, 0.79100E+03 &
, 0.79284E+03, 0.79468E+03, 0.79652E+03, 0.79836E+03, 0.80021E+03, 0.80205E+03, 0.80390E+03, 0.80574E+03, 0.80759E+03, 0.80944E+03 /
 !            7           2
 DATA(QofT(          60 ,J),J=1,501)/   30.9209995269775       &
, 0.33970E+02, 0.37023E+02, 0.40078E+02, 0.43136E+02, 0.46196E+02, 0.49257E+02, 0.52319E+02, 0.55382E+02, 0.58447E+02, 0.61512E+02 &
, 0.64577E+02, 0.67644E+02, 0.70711E+02, 0.73778E+02, 0.76846E+02, 0.79914E+02, 0.82982E+02, 0.86051E+02, 0.89120E+02, 0.92189E+02 &
, 0.95259E+02, 0.98328E+02, 0.10140E+03, 0.10447E+03, 0.10754E+03, 0.11061E+03, 0.11368E+03, 0.11675E+03, 0.11982E+03, 0.12289E+03 &
, 0.12596E+03, 0.12903E+03, 0.13211E+03, 0.13518E+03, 0.13825E+03, 0.14132E+03, 0.14439E+03, 0.14746E+03, 0.15054E+03, 0.15361E+03 &
, 0.15668E+03, 0.15975E+03, 0.16282E+03, 0.16590E+03, 0.16897E+03, 0.17204E+03, 0.17511E+03, 0.17819E+03, 0.18126E+03, 0.18433E+03 &
, 0.18741E+03, 0.19048E+03, 0.19355E+03, 0.19662E+03, 0.19970E+03, 0.20277E+03, 0.20584E+03, 0.20892E+03, 0.21199E+03, 0.21506E+03 &
, 0.21814E+03, 0.22121E+03, 0.22428E+03, 0.22736E+03, 0.23043E+03, 0.23351E+03, 0.23658E+03, 0.23965E+03, 0.24273E+03, 0.24580E+03 &
, 0.24888E+03, 0.25195E+03, 0.25502E+03, 0.25810E+03, 0.26117E+03, 0.26425E+03, 0.26732E+03, 0.27040E+03, 0.27347E+03, 0.27654E+03 &
, 0.27962E+03, 0.28269E+03, 0.28577E+03, 0.28884E+03, 0.29192E+03, 0.29499E+03, 0.29807E+03, 0.30114E+03, 0.30422E+03, 0.30730E+03 &
, 0.31037E+03, 0.31345E+03, 0.31652E+03, 0.31960E+03, 0.32267E+03, 0.32575E+03, 0.32883E+03, 0.33190E+03, 0.33498E+03, 0.33806E+03 &
, 0.34114E+03, 0.34421E+03, 0.34729E+03, 0.35037E+03, 0.35345E+03, 0.35652E+03, 0.35960E+03, 0.36268E+03, 0.36576E+03, 0.36884E+03 &
, 0.37192E+03, 0.37500E+03, 0.37808E+03, 0.38116E+03, 0.38424E+03, 0.38732E+03, 0.39040E+03, 0.39348E+03, 0.39657E+03, 0.39965E+03 &
, 0.40273E+03, 0.40582E+03, 0.40890E+03, 0.41199E+03, 0.41507E+03, 0.41816E+03, 0.42124E+03, 0.42433E+03, 0.42742E+03, 0.43050E+03 &
, 0.43359E+03, 0.43668E+03, 0.43977E+03, 0.44286E+03, 0.44595E+03, 0.44904E+03, 0.45214E+03, 0.45523E+03, 0.45832E+03, 0.46142E+03 &
, 0.46451E+03, 0.46761E+03, 0.47071E+03, 0.47381E+03, 0.47691E+03, 0.48001E+03, 0.48311E+03, 0.48621E+03, 0.48931E+03, 0.49242E+03 &
, 0.49552E+03, 0.49863E+03, 0.50173E+03, 0.50484E+03, 0.50795E+03, 0.51106E+03, 0.51417E+03, 0.51729E+03, 0.52040E+03, 0.52352E+03 &
, 0.52663E+03, 0.52975E+03, 0.53287E+03, 0.53599E+03, 0.53911E+03, 0.54224E+03, 0.54536E+03, 0.54849E+03, 0.55162E+03, 0.55475E+03 &
, 0.55788E+03, 0.56101E+03, 0.56414E+03, 0.56728E+03, 0.57042E+03, 0.57356E+03, 0.57670E+03, 0.57984E+03, 0.58299E+03, 0.58613E+03 &
, 0.58928E+03, 0.59243E+03, 0.59558E+03, 0.59874E+03, 0.60189E+03, 0.60505E+03, 0.60821E+03, 0.61137E+03, 0.61453E+03, 0.61770E+03 &
, 0.62087E+03, 0.62404E+03, 0.62721E+03, 0.63038E+03, 0.63356E+03, 0.63674E+03, 0.63992E+03, 0.64310E+03, 0.64629E+03, 0.64947E+03 &
, 0.65266E+03, 0.65585E+03, 0.65905E+03, 0.66225E+03, 0.66545E+03, 0.66865E+03, 0.67185E+03, 0.67506E+03, 0.67827E+03, 0.68148E+03 &
, 0.68469E+03, 0.68791E+03, 0.69113E+03, 0.69435E+03, 0.69758E+03, 0.70081E+03, 0.70404E+03, 0.70727E+03, 0.71051E+03, 0.71374E+03 &
, 0.71699E+03, 0.72023E+03, 0.72348E+03, 0.72673E+03, 0.72998E+03, 0.73324E+03, 0.73649E+03, 0.73976E+03, 0.74302E+03, 0.74629E+03 &
, 0.74956E+03, 0.75283E+03, 0.75611E+03, 0.75939E+03, 0.76267E+03, 0.76596E+03, 0.76925E+03, 0.77254E+03, 0.77584E+03, 0.77914E+03 &
, 0.78244E+03, 0.78575E+03, 0.78905E+03, 0.79237E+03, 0.79568E+03, 0.79900E+03, 0.80232E+03, 0.80565E+03, 0.80898E+03, 0.81231E+03 &
, 0.81565E+03, 0.81899E+03, 0.82233E+03, 0.82567E+03, 0.82902E+03, 0.83238E+03, 0.83573E+03, 0.83910E+03, 0.84246E+03, 0.84583E+03 &
, 0.84920E+03, 0.85257E+03, 0.85595E+03, 0.85933E+03, 0.86272E+03, 0.86611E+03, 0.86950E+03, 0.87290E+03, 0.87630E+03, 0.87971E+03 &
, 0.88312E+03, 0.88653E+03, 0.88994E+03, 0.89336E+03, 0.89679E+03, 0.90022E+03, 0.90365E+03, 0.90708E+03, 0.91052E+03, 0.91397E+03 &
, 0.91741E+03, 0.92087E+03, 0.92432E+03, 0.92778E+03, 0.93124E+03, 0.93471E+03, 0.93818E+03, 0.94166E+03, 0.94514E+03, 0.94862E+03 &
, 0.95211E+03, 0.95560E+03, 0.95910E+03, 0.96260E+03, 0.96610E+03, 0.96961E+03, 0.97312E+03, 0.97664E+03, 0.98016E+03, 0.98369E+03 &
, 0.98722E+03, 0.99075E+03, 0.99429E+03, 0.99783E+03, 0.10014E+04, 0.10049E+04, 0.10085E+04, 0.10120E+04, 0.10156E+04, 0.10192E+04 &
, 0.10228E+04, 0.10263E+04, 0.10299E+04, 0.10335E+04, 0.10371E+04, 0.10407E+04, 0.10443E+04, 0.10479E+04, 0.10515E+04, 0.10551E+04 &
, 0.10587E+04, 0.10623E+04, 0.10660E+04, 0.10696E+04, 0.10732E+04, 0.10769E+04, 0.10805E+04, 0.10842E+04, 0.10878E+04, 0.10915E+04 &
, 0.10951E+04, 0.10988E+04, 0.11025E+04, 0.11062E+04, 0.11098E+04, 0.11135E+04, 0.11172E+04, 0.11209E+04, 0.11246E+04, 0.11283E+04 &
, 0.11320E+04, 0.11357E+04, 0.11395E+04, 0.11432E+04, 0.11469E+04, 0.11506E+04, 0.11544E+04, 0.11581E+04, 0.11619E+04, 0.11656E+04 &
, 0.11694E+04, 0.11731E+04, 0.11769E+04, 0.11807E+04, 0.11844E+04, 0.11882E+04, 0.11920E+04, 0.11958E+04, 0.11996E+04, 0.12034E+04 &
, 0.12072E+04, 0.12110E+04, 0.12148E+04, 0.12186E+04, 0.12225E+04, 0.12263E+04, 0.12301E+04, 0.12340E+04, 0.12378E+04, 0.12417E+04 &
, 0.12455E+04, 0.12494E+04, 0.12532E+04, 0.12571E+04, 0.12610E+04, 0.12648E+04, 0.12687E+04, 0.12726E+04, 0.12765E+04, 0.12804E+04 &
, 0.12843E+04, 0.12882E+04, 0.12921E+04, 0.12961E+04, 0.13000E+04, 0.13039E+04, 0.13078E+04, 0.13118E+04, 0.13157E+04, 0.13197E+04 &
, 0.13236E+04, 0.13276E+04, 0.13315E+04, 0.13355E+04, 0.13395E+04, 0.13435E+04, 0.13475E+04, 0.13514E+04, 0.13554E+04, 0.13594E+04 &
, 0.13634E+04, 0.13674E+04, 0.13715E+04, 0.13755E+04, 0.13795E+04, 0.13835E+04, 0.13876E+04, 0.13916E+04, 0.13957E+04, 0.13997E+04 &
, 0.14038E+04, 0.14078E+04, 0.14119E+04, 0.14160E+04, 0.14200E+04, 0.14241E+04, 0.14282E+04, 0.14323E+04, 0.14364E+04, 0.14405E+04 &
, 0.14446E+04, 0.14487E+04, 0.14528E+04, 0.14570E+04, 0.14611E+04, 0.14652E+04, 0.14694E+04, 0.14735E+04, 0.14777E+04, 0.14818E+04 &
, 0.14860E+04, 0.14902E+04, 0.14943E+04, 0.14985E+04, 0.15027E+04, 0.15069E+04, 0.15111E+04, 0.15153E+04, 0.15195E+04, 0.15237E+04 &
, 0.15279E+04, 0.15321E+04, 0.15363E+04, 0.15406E+04, 0.15448E+04, 0.15490E+04, 0.15533E+04, 0.15575E+04, 0.15618E+04, 0.15661E+04 &
, 0.15703E+04, 0.15746E+04, 0.15789E+04, 0.15832E+04, 0.15874E+04, 0.15917E+04, 0.15960E+04, 0.16003E+04, 0.16047E+04, 0.16090E+04 &
, 0.16133E+04, 0.16176E+04, 0.16220E+04, 0.16263E+04, 0.16306E+04, 0.16350E+04, 0.16393E+04, 0.16437E+04, 0.16481E+04, 0.16524E+04 &
, 0.16568E+04, 0.16612E+04, 0.16656E+04, 0.16700E+04, 0.16744E+04, 0.16788E+04, 0.16832E+04, 0.16876E+04, 0.16920E+04, 0.16964E+04 &
, 0.17009E+04, 0.17053E+04, 0.17097E+04, 0.17142E+04, 0.17186E+04, 0.17231E+04, 0.17276E+04, 0.17320E+04, 0.17365E+04, 0.17410E+04 &
, 0.17455E+04, 0.17500E+04, 0.17545E+04, 0.17590E+04, 0.17635E+04, 0.17680E+04, 0.17725E+04, 0.17770E+04, 0.17816E+04, 0.17861E+04 /
 !            7           3
 DATA(QofT(          61 ,J),J=1,501)/   180.809997558594       &
, 0.19861E+03, 0.21644E+03, 0.23428E+03, 0.25213E+03, 0.27000E+03, 0.28787E+03, 0.30575E+03, 0.32364E+03, 0.34153E+03, 0.35942E+03 &
, 0.37732E+03, 0.39523E+03, 0.41313E+03, 0.43104E+03, 0.44896E+03, 0.46687E+03, 0.48479E+03, 0.50270E+03, 0.52062E+03, 0.53854E+03 &
, 0.55647E+03, 0.57439E+03, 0.59231E+03, 0.61024E+03, 0.62817E+03, 0.64610E+03, 0.66402E+03, 0.68195E+03, 0.69988E+03, 0.71782E+03 &
, 0.73575E+03, 0.75368E+03, 0.77161E+03, 0.78955E+03, 0.80748E+03, 0.82542E+03, 0.84335E+03, 0.86129E+03, 0.87922E+03, 0.89716E+03 &
, 0.91510E+03, 0.93304E+03, 0.95098E+03, 0.96891E+03, 0.98685E+03, 0.10048E+04, 0.10227E+04, 0.10407E+04, 0.10586E+04, 0.10766E+04 &
, 0.10945E+04, 0.11124E+04, 0.11304E+04, 0.11483E+04, 0.11663E+04, 0.11842E+04, 0.12022E+04, 0.12201E+04, 0.12380E+04, 0.12560E+04 &
, 0.12739E+04, 0.12919E+04, 0.13098E+04, 0.13278E+04, 0.13457E+04, 0.13637E+04, 0.13816E+04, 0.13996E+04, 0.14175E+04, 0.14355E+04 &
, 0.14534E+04, 0.14714E+04, 0.14893E+04, 0.15073E+04, 0.15252E+04, 0.15432E+04, 0.15611E+04, 0.15791E+04, 0.15970E+04, 0.16150E+04 &
, 0.16329E+04, 0.16509E+04, 0.16688E+04, 0.16868E+04, 0.17047E+04, 0.17227E+04, 0.17406E+04, 0.17586E+04, 0.17766E+04, 0.17945E+04 &
, 0.18125E+04, 0.18304E+04, 0.18484E+04, 0.18663E+04, 0.18843E+04, 0.19023E+04, 0.19202E+04, 0.19382E+04, 0.19562E+04, 0.19741E+04 &
, 0.19921E+04, 0.20101E+04, 0.20280E+04, 0.20460E+04, 0.20640E+04, 0.20819E+04, 0.20999E+04, 0.21179E+04, 0.21359E+04, 0.21538E+04 &
, 0.21718E+04, 0.21898E+04, 0.22078E+04, 0.22258E+04, 0.22437E+04, 0.22617E+04, 0.22797E+04, 0.22977E+04, 0.23157E+04, 0.23337E+04 &
, 0.23517E+04, 0.23697E+04, 0.23877E+04, 0.24057E+04, 0.24237E+04, 0.24417E+04, 0.24597E+04, 0.24778E+04, 0.24958E+04, 0.25138E+04 &
, 0.25318E+04, 0.25499E+04, 0.25679E+04, 0.25859E+04, 0.26040E+04, 0.26220E+04, 0.26401E+04, 0.26581E+04, 0.26762E+04, 0.26942E+04 &
, 0.27123E+04, 0.27304E+04, 0.27485E+04, 0.27665E+04, 0.27846E+04, 0.28027E+04, 0.28208E+04, 0.28389E+04, 0.28570E+04, 0.28751E+04 &
, 0.28932E+04, 0.29114E+04, 0.29295E+04, 0.29476E+04, 0.29658E+04, 0.29839E+04, 0.30021E+04, 0.30202E+04, 0.30384E+04, 0.30566E+04 &
, 0.30748E+04, 0.30929E+04, 0.31111E+04, 0.31293E+04, 0.31476E+04, 0.31658E+04, 0.31840E+04, 0.32022E+04, 0.32205E+04, 0.32387E+04 &
, 0.32570E+04, 0.32753E+04, 0.32935E+04, 0.33118E+04, 0.33301E+04, 0.33484E+04, 0.33668E+04, 0.33851E+04, 0.34034E+04, 0.34218E+04 &
, 0.34401E+04, 0.34585E+04, 0.34769E+04, 0.34952E+04, 0.35136E+04, 0.35320E+04, 0.35505E+04, 0.35689E+04, 0.35873E+04, 0.36058E+04 &
, 0.36243E+04, 0.36427E+04, 0.36612E+04, 0.36797E+04, 0.36982E+04, 0.37168E+04, 0.37353E+04, 0.37538E+04, 0.37724E+04, 0.37910E+04 &
, 0.38096E+04, 0.38282E+04, 0.38468E+04, 0.38654E+04, 0.38840E+04, 0.39027E+04, 0.39214E+04, 0.39401E+04, 0.39588E+04, 0.39775E+04 &
, 0.39962E+04, 0.40149E+04, 0.40337E+04, 0.40525E+04, 0.40712E+04, 0.40900E+04, 0.41089E+04, 0.41277E+04, 0.41465E+04, 0.41654E+04 &
, 0.41843E+04, 0.42032E+04, 0.42221E+04, 0.42410E+04, 0.42599E+04, 0.42789E+04, 0.42979E+04, 0.43169E+04, 0.43359E+04, 0.43549E+04 &
, 0.43739E+04, 0.43930E+04, 0.44121E+04, 0.44312E+04, 0.44503E+04, 0.44694E+04, 0.44886E+04, 0.45077E+04, 0.45269E+04, 0.45461E+04 &
, 0.45653E+04, 0.45846E+04, 0.46038E+04, 0.46231E+04, 0.46424E+04, 0.46617E+04, 0.46810E+04, 0.47004E+04, 0.47198E+04, 0.47392E+04 &
, 0.47586E+04, 0.47780E+04, 0.47974E+04, 0.48169E+04, 0.48364E+04, 0.48559E+04, 0.48754E+04, 0.48950E+04, 0.49146E+04, 0.49341E+04 &
, 0.49538E+04, 0.49734E+04, 0.49930E+04, 0.50127E+04, 0.50324E+04, 0.50521E+04, 0.50719E+04, 0.50916E+04, 0.51114E+04, 0.51312E+04 &
, 0.51510E+04, 0.51708E+04, 0.51907E+04, 0.52106E+04, 0.52305E+04, 0.52504E+04, 0.52704E+04, 0.52904E+04, 0.53104E+04, 0.53304E+04 &
, 0.53504E+04, 0.53705E+04, 0.53906E+04, 0.54107E+04, 0.54308E+04, 0.54510E+04, 0.54711E+04, 0.54913E+04, 0.55116E+04, 0.55318E+04 &
, 0.55521E+04, 0.55724E+04, 0.55927E+04, 0.56130E+04, 0.56334E+04, 0.56538E+04, 0.56742E+04, 0.56946E+04, 0.57151E+04, 0.57356E+04 &
, 0.57561E+04, 0.57766E+04, 0.57972E+04, 0.58178E+04, 0.58384E+04, 0.58590E+04, 0.58797E+04, 0.59003E+04, 0.59211E+04, 0.59418E+04 &
, 0.59625E+04, 0.59833E+04, 0.60041E+04, 0.60249E+04, 0.60458E+04, 0.60667E+04, 0.60876E+04, 0.61085E+04, 0.61295E+04, 0.61505E+04 &
, 0.61715E+04, 0.61925E+04, 0.62135E+04, 0.62346E+04, 0.62557E+04, 0.62769E+04, 0.62980E+04, 0.63192E+04, 0.63404E+04, 0.63617E+04 &
, 0.63829E+04, 0.64042E+04, 0.64255E+04, 0.64469E+04, 0.64683E+04, 0.64897E+04, 0.65111E+04, 0.65325E+04, 0.65540E+04, 0.65755E+04 &
, 0.65970E+04, 0.66186E+04, 0.66402E+04, 0.66618E+04, 0.66834E+04, 0.67051E+04, 0.67268E+04, 0.67485E+04, 0.67703E+04, 0.67920E+04 &
, 0.68138E+04, 0.68357E+04, 0.68575E+04, 0.68794E+04, 0.69013E+04, 0.69232E+04, 0.69452E+04, 0.69672E+04, 0.69892E+04, 0.70113E+04 &
, 0.70333E+04, 0.70554E+04, 0.70776E+04, 0.70997E+04, 0.71219E+04, 0.71441E+04, 0.71664E+04, 0.71887E+04, 0.72110E+04, 0.72333E+04 &
, 0.72556E+04, 0.72780E+04, 0.73004E+04, 0.73229E+04, 0.73453E+04, 0.73678E+04, 0.73904E+04, 0.74129E+04, 0.74355E+04, 0.74581E+04 &
, 0.74808E+04, 0.75034E+04, 0.75261E+04, 0.75489E+04, 0.75716E+04, 0.75944E+04, 0.76172E+04, 0.76400E+04, 0.76629E+04, 0.76858E+04 &
, 0.77087E+04, 0.77317E+04, 0.77547E+04, 0.77777E+04, 0.78008E+04, 0.78238E+04, 0.78469E+04, 0.78701E+04, 0.78932E+04, 0.79164E+04 &
, 0.79396E+04, 0.79629E+04, 0.79862E+04, 0.80095E+04, 0.80328E+04, 0.80562E+04, 0.80796E+04, 0.81030E+04, 0.81265E+04, 0.81500E+04 &
, 0.81735E+04, 0.81970E+04, 0.82206E+04, 0.82442E+04, 0.82678E+04, 0.82915E+04, 0.83152E+04, 0.83389E+04, 0.83627E+04, 0.83865E+04 &
, 0.84103E+04, 0.84341E+04, 0.84580E+04, 0.84819E+04, 0.85058E+04, 0.85298E+04, 0.85538E+04, 0.85778E+04, 0.86019E+04, 0.86260E+04 &
, 0.86501E+04, 0.86742E+04, 0.86984E+04, 0.87226E+04, 0.87469E+04, 0.87711E+04, 0.87954E+04, 0.88198E+04, 0.88441E+04, 0.88685E+04 &
, 0.88929E+04, 0.89174E+04, 0.89419E+04, 0.89664E+04, 0.89909E+04, 0.90155E+04, 0.90401E+04, 0.90648E+04, 0.90894E+04, 0.91141E+04 &
, 0.91388E+04, 0.91636E+04, 0.91884E+04, 0.92132E+04, 0.92381E+04, 0.92630E+04, 0.92879E+04, 0.93128E+04, 0.93378E+04, 0.93628E+04 &
, 0.93878E+04, 0.94129E+04, 0.94380E+04, 0.94631E+04, 0.94883E+04, 0.95135E+04, 0.95387E+04, 0.95640E+04, 0.95893E+04, 0.96146E+04 &
, 0.96399E+04, 0.96653E+04, 0.96907E+04, 0.97162E+04, 0.97416E+04, 0.97672E+04, 0.97927E+04, 0.98183E+04, 0.98439E+04, 0.98695E+04 &
, 0.98952E+04, 0.99209E+04, 0.99466E+04, 0.99723E+04, 0.99981E+04, 0.10024E+05, 0.10050E+05, 0.10076E+05, 0.10102E+05, 0.10128E+05 &
, 0.10154E+05, 0.10180E+05, 0.10206E+05, 0.10232E+05, 0.10258E+05, 0.10284E+05, 0.10310E+05, 0.10336E+05, 0.10363E+05, 0.10389E+05 /
 !            7           4
 DATA(QofT(          62 ,J),J=1,501)/   16.8780002593994       &
, 0.18496E+02, 0.20117E+02, 0.21739E+02, 0.23363E+02, 0.24988E+02, 0.26614E+02, 0.28240E+02, 0.29867E+02, 0.31495E+02, 0.33123E+02 &
, 0.34751E+02, 0.36380E+02, 0.38009E+02, 0.39638E+02, 0.41268E+02, 0.42897E+02, 0.44527E+02, 0.46157E+02, 0.47788E+02, 0.49418E+02 &
, 0.51049E+02, 0.52679E+02, 0.54310E+02, 0.55941E+02, 0.57572E+02, 0.59203E+02, 0.60834E+02, 0.62465E+02, 0.64097E+02, 0.65728E+02 &
, 0.67360E+02, 0.68991E+02, 0.70623E+02, 0.72254E+02, 0.73886E+02, 0.75518E+02, 0.77150E+02, 0.78781E+02, 0.80413E+02, 0.82045E+02 &
, 0.83677E+02, 0.85309E+02, 0.86941E+02, 0.88574E+02, 0.90206E+02, 0.91838E+02, 0.93470E+02, 0.95103E+02, 0.96735E+02, 0.98367E+02 &
, 0.10000E+03, 0.10163E+03, 0.10326E+03, 0.10490E+03, 0.10653E+03, 0.10816E+03, 0.10979E+03, 0.11143E+03, 0.11306E+03, 0.11469E+03 &
, 0.11633E+03, 0.11796E+03, 0.11959E+03, 0.12122E+03, 0.12286E+03, 0.12449E+03, 0.12612E+03, 0.12776E+03, 0.12939E+03, 0.13102E+03 &
, 0.13265E+03, 0.13429E+03, 0.13592E+03, 0.13755E+03, 0.13919E+03, 0.14082E+03, 0.14245E+03, 0.14409E+03, 0.14572E+03, 0.14735E+03 &
, 0.14899E+03, 0.15062E+03, 0.15225E+03, 0.15389E+03, 0.15552E+03, 0.15715E+03, 0.15879E+03, 0.16042E+03, 0.16206E+03, 0.16369E+03 &
, 0.16532E+03, 0.16696E+03, 0.16859E+03, 0.17023E+03, 0.17186E+03, 0.17350E+03, 0.17513E+03, 0.17676E+03, 0.17840E+03, 0.18003E+03 &
, 0.18167E+03, 0.18330E+03, 0.18494E+03, 0.18658E+03, 0.18821E+03, 0.18985E+03, 0.19148E+03, 0.19312E+03, 0.19475E+03, 0.19639E+03 &
, 0.19803E+03, 0.19966E+03, 0.20130E+03, 0.20294E+03, 0.20458E+03, 0.20621E+03, 0.20785E+03, 0.20949E+03, 0.21113E+03, 0.21277E+03 &
, 0.21440E+03, 0.21604E+03, 0.21768E+03, 0.21932E+03, 0.22096E+03, 0.22260E+03, 0.22424E+03, 0.22588E+03, 0.22753E+03, 0.22917E+03 &
, 0.23081E+03, 0.23245E+03, 0.23409E+03, 0.23574E+03, 0.23738E+03, 0.23903E+03, 0.24067E+03, 0.24232E+03, 0.24396E+03, 0.24561E+03 &
, 0.24725E+03, 0.24890E+03, 0.25055E+03, 0.25220E+03, 0.25385E+03, 0.25549E+03, 0.25714E+03, 0.25879E+03, 0.26045E+03, 0.26210E+03 &
, 0.26375E+03, 0.26540E+03, 0.26706E+03, 0.26871E+03, 0.27037E+03, 0.27202E+03, 0.27368E+03, 0.27533E+03, 0.27699E+03, 0.27865E+03 &
, 0.28031E+03, 0.28197E+03, 0.28363E+03, 0.28529E+03, 0.28696E+03, 0.28862E+03, 0.29028E+03, 0.29195E+03, 0.29362E+03, 0.29528E+03 &
, 0.29695E+03, 0.29862E+03, 0.30029E+03, 0.30196E+03, 0.30363E+03, 0.30531E+03, 0.30698E+03, 0.30866E+03, 0.31033E+03, 0.31201E+03 &
, 0.31369E+03, 0.31537E+03, 0.31705E+03, 0.31873E+03, 0.32041E+03, 0.32209E+03, 0.32378E+03, 0.32546E+03, 0.32715E+03, 0.32884E+03 &
, 0.33053E+03, 0.33222E+03, 0.33391E+03, 0.33561E+03, 0.33730E+03, 0.33900E+03, 0.34069E+03, 0.34239E+03, 0.34409E+03, 0.34579E+03 &
, 0.34750E+03, 0.34920E+03, 0.35091E+03, 0.35261E+03, 0.35432E+03, 0.35603E+03, 0.35774E+03, 0.35946E+03, 0.36117E+03, 0.36288E+03 &
, 0.36460E+03, 0.36632E+03, 0.36804E+03, 0.36976E+03, 0.37149E+03, 0.37321E+03, 0.37494E+03, 0.37666E+03, 0.37839E+03, 0.38012E+03 &
, 0.38186E+03, 0.38359E+03, 0.38533E+03, 0.38707E+03, 0.38880E+03, 0.39055E+03, 0.39229E+03, 0.39403E+03, 0.39578E+03, 0.39753E+03 &
, 0.39928E+03, 0.40103E+03, 0.40278E+03, 0.40454E+03, 0.40629E+03, 0.40805E+03, 0.40981E+03, 0.41157E+03, 0.41334E+03, 0.41510E+03 &
, 0.41687E+03, 0.41864E+03, 0.42041E+03, 0.42219E+03, 0.42396E+03, 0.42574E+03, 0.42752E+03, 0.42930E+03, 0.43108E+03, 0.43287E+03 &
, 0.43465E+03, 0.43644E+03, 0.43823E+03, 0.44003E+03, 0.44182E+03, 0.44362E+03, 0.44542E+03, 0.44722E+03, 0.44902E+03, 0.45082E+03 &
, 0.45263E+03, 0.45444E+03, 0.45625E+03, 0.45807E+03, 0.45988E+03, 0.46170E+03, 0.46352E+03, 0.46534E+03, 0.46716E+03, 0.46899E+03 &
, 0.47082E+03, 0.47265E+03, 0.47448E+03, 0.47632E+03, 0.47815E+03, 0.47999E+03, 0.48183E+03, 0.48368E+03, 0.48552E+03, 0.48737E+03 &
, 0.48922E+03, 0.49107E+03, 0.49293E+03, 0.49478E+03, 0.49664E+03, 0.49851E+03, 0.50037E+03, 0.50224E+03, 0.50410E+03, 0.50597E+03 &
, 0.50785E+03, 0.50972E+03, 0.51160E+03, 0.51348E+03, 0.51536E+03, 0.51725E+03, 0.51913E+03, 0.52102E+03, 0.52292E+03, 0.52481E+03 &
, 0.52671E+03, 0.52860E+03, 0.53051E+03, 0.53241E+03, 0.53432E+03, 0.53622E+03, 0.53814E+03, 0.54005E+03, 0.54196E+03, 0.54388E+03 &
, 0.54580E+03, 0.54773E+03, 0.54965E+03, 0.55158E+03, 0.55351E+03, 0.55544E+03, 0.55738E+03, 0.55932E+03, 0.56126E+03, 0.56320E+03 &
, 0.56515E+03, 0.56710E+03, 0.56905E+03, 0.57100E+03, 0.57295E+03, 0.57491E+03, 0.57687E+03, 0.57884E+03, 0.58080E+03, 0.58277E+03 &
, 0.58474E+03, 0.58672E+03, 0.58869E+03, 0.59067E+03, 0.59265E+03, 0.59464E+03, 0.59662E+03, 0.59861E+03, 0.60060E+03, 0.60260E+03 &
, 0.60459E+03, 0.60659E+03, 0.60860E+03, 0.61060E+03, 0.61261E+03, 0.61462E+03, 0.61663E+03, 0.61865E+03, 0.62066E+03, 0.62268E+03 &
, 0.62471E+03, 0.62673E+03, 0.62876E+03, 0.63079E+03, 0.63283E+03, 0.63486E+03, 0.63690E+03, 0.63895E+03, 0.64099E+03, 0.64304E+03 &
, 0.64509E+03, 0.64714E+03, 0.64920E+03, 0.65125E+03, 0.65332E+03, 0.65538E+03, 0.65745E+03, 0.65951E+03, 0.66159E+03, 0.66366E+03 &
, 0.66574E+03, 0.66782E+03, 0.66990E+03, 0.67199E+03, 0.67408E+03, 0.67617E+03, 0.67826E+03, 0.68036E+03, 0.68246E+03, 0.68456E+03 &
, 0.68666E+03, 0.68877E+03, 0.69088E+03, 0.69300E+03, 0.69511E+03, 0.69723E+03, 0.69935E+03, 0.70148E+03, 0.70360E+03, 0.70573E+03 &
, 0.70787E+03, 0.71000E+03, 0.71214E+03, 0.71428E+03, 0.71643E+03, 0.71857E+03, 0.72072E+03, 0.72288E+03, 0.72503E+03, 0.72719E+03 &
, 0.72935E+03, 0.73152E+03, 0.73368E+03, 0.73585E+03, 0.73803E+03, 0.74020E+03, 0.74238E+03, 0.74456E+03, 0.74674E+03, 0.74893E+03 &
, 0.75112E+03, 0.75331E+03, 0.75551E+03, 0.75771E+03, 0.75991E+03, 0.76211E+03, 0.76432E+03, 0.76653E+03, 0.76874E+03, 0.77096E+03 &
, 0.77318E+03, 0.77540E+03, 0.77762E+03, 0.77985E+03, 0.78208E+03, 0.78432E+03, 0.78655E+03, 0.78879E+03, 0.79103E+03, 0.79328E+03 &
, 0.79553E+03, 0.79778E+03, 0.80003E+03, 0.80229E+03, 0.80455E+03, 0.80681E+03, 0.80907E+03, 0.81134E+03, 0.81361E+03, 0.81589E+03 &
, 0.81816E+03, 0.82045E+03, 0.82273E+03, 0.82501E+03, 0.82730E+03, 0.82959E+03, 0.83189E+03, 0.83419E+03, 0.83649E+03, 0.83879E+03 &
, 0.84110E+03, 0.84341E+03, 0.84572E+03, 0.84804E+03, 0.85036E+03, 0.85268E+03, 0.85500E+03, 0.85733E+03, 0.85966E+03, 0.86199E+03 &
, 0.86433E+03, 0.86667E+03, 0.86901E+03, 0.87136E+03, 0.87371E+03, 0.87606E+03, 0.87841E+03, 0.88077E+03, 0.88313E+03, 0.88549E+03 &
, 0.88786E+03, 0.89023E+03, 0.89260E+03, 0.89498E+03, 0.89736E+03, 0.89974E+03, 0.90212E+03, 0.90451E+03, 0.90690E+03, 0.90929E+03 &
, 0.91169E+03, 0.91409E+03, 0.91649E+03, 0.91890E+03, 0.92131E+03, 0.92372E+03, 0.92613E+03, 0.92855E+03, 0.93097E+03, 0.93340E+03 &
, 0.93582E+03, 0.93825E+03, 0.94069E+03, 0.94312E+03, 0.94556E+03, 0.94800E+03, 0.95045E+03, 0.95290E+03, 0.95535E+03, 0.95780E+03 /
 !            7           5
 DATA(QofT(          63 ,J),J=1,501)/   191.020004272461       &
, 0.20990E+03, 0.22881E+03, 0.24773E+03, 0.26667E+03, 0.28561E+03, 0.30457E+03, 0.32354E+03, 0.34251E+03, 0.36149E+03, 0.38047E+03 &
, 0.39946E+03, 0.41845E+03, 0.43744E+03, 0.45644E+03, 0.47543E+03, 0.49444E+03, 0.51344E+03, 0.53244E+03, 0.55145E+03, 0.57046E+03 &
, 0.58947E+03, 0.60848E+03, 0.62749E+03, 0.64651E+03, 0.66552E+03, 0.68454E+03, 0.70355E+03, 0.72257E+03, 0.74159E+03, 0.76061E+03 &
, 0.77963E+03, 0.79865E+03, 0.81767E+03, 0.83669E+03, 0.85572E+03, 0.87474E+03, 0.89376E+03, 0.91279E+03, 0.93181E+03, 0.95084E+03 &
, 0.96987E+03, 0.98889E+03, 0.10079E+04, 0.10269E+04, 0.10460E+04, 0.10650E+04, 0.10840E+04, 0.11031E+04, 0.11221E+04, 0.11411E+04 &
, 0.11602E+04, 0.11792E+04, 0.11982E+04, 0.12172E+04, 0.12363E+04, 0.12553E+04, 0.12743E+04, 0.12934E+04, 0.13124E+04, 0.13314E+04 &
, 0.13505E+04, 0.13695E+04, 0.13886E+04, 0.14076E+04, 0.14266E+04, 0.14457E+04, 0.14647E+04, 0.14837E+04, 0.15028E+04, 0.15218E+04 &
, 0.15408E+04, 0.15599E+04, 0.15789E+04, 0.15980E+04, 0.16170E+04, 0.16360E+04, 0.16551E+04, 0.16741E+04, 0.16932E+04, 0.17122E+04 &
, 0.17313E+04, 0.17503E+04, 0.17693E+04, 0.17884E+04, 0.18074E+04, 0.18265E+04, 0.18455E+04, 0.18646E+04, 0.18836E+04, 0.19027E+04 &
, 0.19217E+04, 0.19408E+04, 0.19598E+04, 0.19789E+04, 0.19979E+04, 0.20170E+04, 0.20360E+04, 0.20551E+04, 0.20741E+04, 0.20932E+04 &
, 0.21122E+04, 0.21313E+04, 0.21504E+04, 0.21694E+04, 0.21885E+04, 0.22076E+04, 0.22266E+04, 0.22457E+04, 0.22648E+04, 0.22838E+04 &
, 0.23029E+04, 0.23220E+04, 0.23411E+04, 0.23602E+04, 0.23792E+04, 0.23983E+04, 0.24174E+04, 0.24365E+04, 0.24556E+04, 0.24747E+04 &
, 0.24938E+04, 0.25129E+04, 0.25320E+04, 0.25511E+04, 0.25702E+04, 0.25893E+04, 0.26085E+04, 0.26276E+04, 0.26467E+04, 0.26659E+04 &
, 0.26850E+04, 0.27041E+04, 0.27233E+04, 0.27424E+04, 0.27616E+04, 0.27807E+04, 0.27999E+04, 0.28191E+04, 0.28382E+04, 0.28574E+04 &
, 0.28766E+04, 0.28958E+04, 0.29150E+04, 0.29342E+04, 0.29534E+04, 0.29726E+04, 0.29918E+04, 0.30110E+04, 0.30303E+04, 0.30495E+04 &
, 0.30688E+04, 0.30880E+04, 0.31073E+04, 0.31265E+04, 0.31458E+04, 0.31651E+04, 0.31844E+04, 0.32037E+04, 0.32230E+04, 0.32423E+04 &
, 0.32616E+04, 0.32810E+04, 0.33003E+04, 0.33197E+04, 0.33390E+04, 0.33584E+04, 0.33778E+04, 0.33972E+04, 0.34166E+04, 0.34360E+04 &
, 0.34554E+04, 0.34748E+04, 0.34943E+04, 0.35137E+04, 0.35332E+04, 0.35526E+04, 0.35721E+04, 0.35916E+04, 0.36111E+04, 0.36307E+04 &
, 0.36502E+04, 0.36697E+04, 0.36893E+04, 0.37088E+04, 0.37284E+04, 0.37480E+04, 0.37676E+04, 0.37872E+04, 0.38069E+04, 0.38265E+04 &
, 0.38462E+04, 0.38659E+04, 0.38855E+04, 0.39052E+04, 0.39250E+04, 0.39447E+04, 0.39644E+04, 0.39842E+04, 0.40040E+04, 0.40237E+04 &
, 0.40435E+04, 0.40634E+04, 0.40832E+04, 0.41030E+04, 0.41229E+04, 0.41428E+04, 0.41627E+04, 0.41826E+04, 0.42025E+04, 0.42225E+04 &
, 0.42424E+04, 0.42624E+04, 0.42824E+04, 0.43024E+04, 0.43225E+04, 0.43425E+04, 0.43626E+04, 0.43827E+04, 0.44028E+04, 0.44229E+04 &
, 0.44430E+04, 0.44632E+04, 0.44834E+04, 0.45035E+04, 0.45238E+04, 0.45440E+04, 0.45642E+04, 0.45845E+04, 0.46048E+04, 0.46251E+04 &
, 0.46454E+04, 0.46658E+04, 0.46862E+04, 0.47065E+04, 0.47269E+04, 0.47474E+04, 0.47678E+04, 0.47883E+04, 0.48088E+04, 0.48293E+04 &
, 0.48498E+04, 0.48704E+04, 0.48910E+04, 0.49115E+04, 0.49322E+04, 0.49528E+04, 0.49735E+04, 0.49941E+04, 0.50148E+04, 0.50356E+04 &
, 0.50563E+04, 0.50771E+04, 0.50979E+04, 0.51187E+04, 0.51395E+04, 0.51604E+04, 0.51813E+04, 0.52022E+04, 0.52231E+04, 0.52441E+04 &
, 0.52650E+04, 0.52860E+04, 0.53071E+04, 0.53281E+04, 0.53492E+04, 0.53703E+04, 0.53914E+04, 0.54125E+04, 0.54337E+04, 0.54549E+04 &
, 0.54761E+04, 0.54973E+04, 0.55186E+04, 0.55399E+04, 0.55612E+04, 0.55825E+04, 0.56039E+04, 0.56253E+04, 0.56467E+04, 0.56681E+04 &
, 0.56896E+04, 0.57111E+04, 0.57326E+04, 0.57541E+04, 0.57757E+04, 0.57973E+04, 0.58189E+04, 0.58405E+04, 0.58622E+04, 0.58839E+04 &
, 0.59056E+04, 0.59274E+04, 0.59491E+04, 0.59709E+04, 0.59928E+04, 0.60146E+04, 0.60365E+04, 0.60584E+04, 0.60803E+04, 0.61023E+04 &
, 0.61243E+04, 0.61463E+04, 0.61683E+04, 0.61904E+04, 0.62125E+04, 0.62346E+04, 0.62568E+04, 0.62790E+04, 0.63012E+04, 0.63234E+04 &
, 0.63457E+04, 0.63680E+04, 0.63903E+04, 0.64126E+04, 0.64350E+04, 0.64574E+04, 0.64798E+04, 0.65023E+04, 0.65248E+04, 0.65473E+04 &
, 0.65698E+04, 0.65924E+04, 0.66150E+04, 0.66376E+04, 0.66603E+04, 0.66830E+04, 0.67057E+04, 0.67285E+04, 0.67512E+04, 0.67740E+04 &
, 0.67969E+04, 0.68197E+04, 0.68426E+04, 0.68655E+04, 0.68885E+04, 0.69115E+04, 0.69345E+04, 0.69575E+04, 0.69806E+04, 0.70037E+04 &
, 0.70268E+04, 0.70500E+04, 0.70732E+04, 0.70964E+04, 0.71196E+04, 0.71429E+04, 0.71662E+04, 0.71896E+04, 0.72129E+04, 0.72363E+04 &
, 0.72598E+04, 0.72832E+04, 0.73067E+04, 0.73302E+04, 0.73538E+04, 0.73774E+04, 0.74010E+04, 0.74246E+04, 0.74483E+04, 0.74720E+04 &
, 0.74957E+04, 0.75195E+04, 0.75433E+04, 0.75671E+04, 0.75910E+04, 0.76149E+04, 0.76388E+04, 0.76627E+04, 0.76867E+04, 0.77107E+04 &
, 0.77348E+04, 0.77589E+04, 0.77830E+04, 0.78071E+04, 0.78313E+04, 0.78555E+04, 0.78797E+04, 0.79040E+04, 0.79283E+04, 0.79526E+04 &
, 0.79770E+04, 0.80014E+04, 0.80258E+04, 0.80502E+04, 0.80747E+04, 0.80993E+04, 0.81238E+04, 0.81484E+04, 0.81730E+04, 0.81977E+04 &
, 0.82223E+04, 0.82470E+04, 0.82718E+04, 0.82966E+04, 0.83214E+04, 0.83462E+04, 0.83711E+04, 0.83960E+04, 0.84209E+04, 0.84459E+04 &
, 0.84709E+04, 0.84959E+04, 0.85210E+04, 0.85461E+04, 0.85713E+04, 0.85964E+04, 0.86216E+04, 0.86469E+04, 0.86721E+04, 0.86974E+04 &
, 0.87227E+04, 0.87481E+04, 0.87735E+04, 0.87989E+04, 0.88244E+04, 0.88499E+04, 0.88754E+04, 0.89010E+04, 0.89266E+04, 0.89522E+04 &
, 0.89779E+04, 0.90036E+04, 0.90293E+04, 0.90550E+04, 0.90808E+04, 0.91067E+04, 0.91325E+04, 0.91584E+04, 0.91843E+04, 0.92103E+04 &
, 0.92363E+04, 0.92623E+04, 0.92884E+04, 0.93145E+04, 0.93406E+04, 0.93668E+04, 0.93930E+04, 0.94192E+04, 0.94454E+04, 0.94717E+04 &
, 0.94981E+04, 0.95244E+04, 0.95508E+04, 0.95773E+04, 0.96037E+04, 0.96302E+04, 0.96567E+04, 0.96833E+04, 0.97099E+04, 0.97365E+04 &
, 0.97632E+04, 0.97899E+04, 0.98166E+04, 0.98434E+04, 0.98702E+04, 0.98971E+04, 0.99239E+04, 0.99508E+04, 0.99778E+04, 0.10005E+05 &
, 0.10032E+05, 0.10059E+05, 0.10086E+05, 0.10113E+05, 0.10140E+05, 0.10167E+05, 0.10195E+05, 0.10222E+05, 0.10249E+05, 0.10276E+05 &
, 0.10304E+05, 0.10331E+05, 0.10359E+05, 0.10386E+05, 0.10413E+05, 0.10441E+05, 0.10469E+05, 0.10496E+05, 0.10524E+05, 0.10551E+05 &
, 0.10579E+05, 0.10607E+05, 0.10635E+05, 0.10662E+05, 0.10690E+05, 0.10718E+05, 0.10746E+05, 0.10774E+05, 0.10802E+05, 0.10830E+05 &
, 0.10858E+05, 0.10886E+05, 0.10914E+05, 0.10942E+05, 0.10971E+05, 0.10999E+05, 0.11027E+05, 0.11055E+05, 0.11084E+05, 0.11112E+05 /
 !            7           6
 DATA(QofT(          64 ,J),J=1,501)/   558.619995117188       &
, 0.61370E+03, 0.66885E+03, 0.72404E+03, 0.77928E+03, 0.83455E+03, 0.88985E+03, 0.94517E+03, 0.10005E+04, 0.10559E+04, 0.11112E+04 &
, 0.11666E+04, 0.12220E+04, 0.12774E+04, 0.13328E+04, 0.13882E+04, 0.14437E+04, 0.14991E+04, 0.15545E+04, 0.16100E+04, 0.16654E+04 &
, 0.17209E+04, 0.17763E+04, 0.18318E+04, 0.18872E+04, 0.19427E+04, 0.19982E+04, 0.20536E+04, 0.21091E+04, 0.21646E+04, 0.22201E+04 &
, 0.22755E+04, 0.23310E+04, 0.23865E+04, 0.24420E+04, 0.24975E+04, 0.25530E+04, 0.26085E+04, 0.26639E+04, 0.27194E+04, 0.27749E+04 &
, 0.28304E+04, 0.28859E+04, 0.29414E+04, 0.29969E+04, 0.30524E+04, 0.31079E+04, 0.31634E+04, 0.32190E+04, 0.32745E+04, 0.33300E+04 &
, 0.33855E+04, 0.34410E+04, 0.34965E+04, 0.35520E+04, 0.36075E+04, 0.36631E+04, 0.37186E+04, 0.37741E+04, 0.38296E+04, 0.38851E+04 &
, 0.39406E+04, 0.39962E+04, 0.40517E+04, 0.41072E+04, 0.41627E+04, 0.42183E+04, 0.42738E+04, 0.43293E+04, 0.43849E+04, 0.44404E+04 &
, 0.44959E+04, 0.45515E+04, 0.46070E+04, 0.46625E+04, 0.47181E+04, 0.47736E+04, 0.48291E+04, 0.48847E+04, 0.49402E+04, 0.49958E+04 &
, 0.50513E+04, 0.51068E+04, 0.51624E+04, 0.52179E+04, 0.52735E+04, 0.53290E+04, 0.53846E+04, 0.54402E+04, 0.54957E+04, 0.55513E+04 &
, 0.56068E+04, 0.56624E+04, 0.57180E+04, 0.57735E+04, 0.58291E+04, 0.58847E+04, 0.59403E+04, 0.59958E+04, 0.60514E+04, 0.61070E+04 &
, 0.61626E+04, 0.62182E+04, 0.62738E+04, 0.63294E+04, 0.63850E+04, 0.64406E+04, 0.64962E+04, 0.65518E+04, 0.66074E+04, 0.66631E+04 &
, 0.67187E+04, 0.67743E+04, 0.68300E+04, 0.68856E+04, 0.69413E+04, 0.69969E+04, 0.70526E+04, 0.71083E+04, 0.71640E+04, 0.72197E+04 &
, 0.72754E+04, 0.73311E+04, 0.73868E+04, 0.74425E+04, 0.74982E+04, 0.75540E+04, 0.76097E+04, 0.76655E+04, 0.77213E+04, 0.77771E+04 &
, 0.78329E+04, 0.78887E+04, 0.79445E+04, 0.80003E+04, 0.80561E+04, 0.81120E+04, 0.81679E+04, 0.82238E+04, 0.82796E+04, 0.83356E+04 &
, 0.83915E+04, 0.84474E+04, 0.85034E+04, 0.85594E+04, 0.86153E+04, 0.86713E+04, 0.87274E+04, 0.87834E+04, 0.88395E+04, 0.88955E+04 &
, 0.89516E+04, 0.90078E+04, 0.90639E+04, 0.91201E+04, 0.91762E+04, 0.92324E+04, 0.92887E+04, 0.93449E+04, 0.94012E+04, 0.94574E+04 &
, 0.95138E+04, 0.95701E+04, 0.96265E+04, 0.96828E+04, 0.97392E+04, 0.97957E+04, 0.98521E+04, 0.99086E+04, 0.99652E+04, 0.10022E+05 &
, 0.10078E+05, 0.10135E+05, 0.10192E+05, 0.10248E+05, 0.10305E+05, 0.10362E+05, 0.10418E+05, 0.10475E+05, 0.10532E+05, 0.10589E+05 &
, 0.10646E+05, 0.10703E+05, 0.10760E+05, 0.10816E+05, 0.10874E+05, 0.10931E+05, 0.10988E+05, 0.11045E+05, 0.11102E+05, 0.11159E+05 &
, 0.11216E+05, 0.11274E+05, 0.11331E+05, 0.11388E+05, 0.11446E+05, 0.11503E+05, 0.11561E+05, 0.11618E+05, 0.11676E+05, 0.11733E+05 &
, 0.11791E+05, 0.11849E+05, 0.11906E+05, 0.11964E+05, 0.12022E+05, 0.12080E+05, 0.12138E+05, 0.12196E+05, 0.12254E+05, 0.12312E+05 &
, 0.12370E+05, 0.12428E+05, 0.12486E+05, 0.12544E+05, 0.12603E+05, 0.12661E+05, 0.12719E+05, 0.12778E+05, 0.12836E+05, 0.12895E+05 &
, 0.12953E+05, 0.13012E+05, 0.13071E+05, 0.13129E+05, 0.13188E+05, 0.13247E+05, 0.13306E+05, 0.13365E+05, 0.13424E+05, 0.13483E+05 &
, 0.13542E+05, 0.13601E+05, 0.13660E+05, 0.13720E+05, 0.13779E+05, 0.13838E+05, 0.13898E+05, 0.13957E+05, 0.14017E+05, 0.14076E+05 &
, 0.14136E+05, 0.14196E+05, 0.14256E+05, 0.14316E+05, 0.14375E+05, 0.14435E+05, 0.14495E+05, 0.14556E+05, 0.14616E+05, 0.14676E+05 &
, 0.14736E+05, 0.14797E+05, 0.14857E+05, 0.14918E+05, 0.14978E+05, 0.15039E+05, 0.15099E+05, 0.15160E+05, 0.15221E+05, 0.15282E+05 &
, 0.15343E+05, 0.15404E+05, 0.15465E+05, 0.15526E+05, 0.15587E+05, 0.15648E+05, 0.15710E+05, 0.15771E+05, 0.15833E+05, 0.15894E+05 &
, 0.15956E+05, 0.16017E+05, 0.16079E+05, 0.16141E+05, 0.16203E+05, 0.16265E+05, 0.16327E+05, 0.16389E+05, 0.16451E+05, 0.16513E+05 &
, 0.16576E+05, 0.16638E+05, 0.16701E+05, 0.16763E+05, 0.16826E+05, 0.16888E+05, 0.16951E+05, 0.17014E+05, 0.17077E+05, 0.17140E+05 &
, 0.17203E+05, 0.17266E+05, 0.17329E+05, 0.17392E+05, 0.17456E+05, 0.17519E+05, 0.17583E+05, 0.17646E+05, 0.17710E+05, 0.17774E+05 &
, 0.17837E+05, 0.17901E+05, 0.17965E+05, 0.18029E+05, 0.18093E+05, 0.18158E+05, 0.18222E+05, 0.18286E+05, 0.18351E+05, 0.18415E+05 &
, 0.18480E+05, 0.18544E+05, 0.18609E+05, 0.18674E+05, 0.18739E+05, 0.18804E+05, 0.18869E+05, 0.18934E+05, 0.18999E+05, 0.19065E+05 &
, 0.19130E+05, 0.19196E+05, 0.19261E+05, 0.19327E+05, 0.19392E+05, 0.19458E+05, 0.19524E+05, 0.19590E+05, 0.19656E+05, 0.19722E+05 &
, 0.19788E+05, 0.19855E+05, 0.19921E+05, 0.19988E+05, 0.20054E+05, 0.20121E+05, 0.20187E+05, 0.20254E+05, 0.20321E+05, 0.20388E+05 &
, 0.20455E+05, 0.20522E+05, 0.20590E+05, 0.20657E+05, 0.20724E+05, 0.20792E+05, 0.20859E+05, 0.20927E+05, 0.20995E+05, 0.21062E+05 &
, 0.21130E+05, 0.21198E+05, 0.21266E+05, 0.21335E+05, 0.21403E+05, 0.21471E+05, 0.21540E+05, 0.21608E+05, 0.21677E+05, 0.21745E+05 &
, 0.21814E+05, 0.21883E+05, 0.21952E+05, 0.22021E+05, 0.22090E+05, 0.22160E+05, 0.22229E+05, 0.22298E+05, 0.22368E+05, 0.22437E+05 &
, 0.22507E+05, 0.22577E+05, 0.22647E+05, 0.22717E+05, 0.22787E+05, 0.22857E+05, 0.22927E+05, 0.22997E+05, 0.23068E+05, 0.23138E+05 &
, 0.23209E+05, 0.23279E+05, 0.23350E+05, 0.23421E+05, 0.23492E+05, 0.23563E+05, 0.23634E+05, 0.23705E+05, 0.23776E+05, 0.23848E+05 &
, 0.23919E+05, 0.23991E+05, 0.24063E+05, 0.24134E+05, 0.24206E+05, 0.24278E+05, 0.24350E+05, 0.24422E+05, 0.24495E+05, 0.24567E+05 &
, 0.24639E+05, 0.24712E+05, 0.24784E+05, 0.24857E+05, 0.24930E+05, 0.25003E+05, 0.25076E+05, 0.25149E+05, 0.25222E+05, 0.25295E+05 &
, 0.25368E+05, 0.25442E+05, 0.25515E+05, 0.25589E+05, 0.25663E+05, 0.25737E+05, 0.25810E+05, 0.25884E+05, 0.25959E+05, 0.26033E+05 &
, 0.26107E+05, 0.26181E+05, 0.26256E+05, 0.26331E+05, 0.26405E+05, 0.26480E+05, 0.26555E+05, 0.26630E+05, 0.26705E+05, 0.26780E+05 &
, 0.26855E+05, 0.26931E+05, 0.27006E+05, 0.27081E+05, 0.27157E+05, 0.27233E+05, 0.27309E+05, 0.27385E+05, 0.27461E+05, 0.27537E+05 &
, 0.27613E+05, 0.27689E+05, 0.27766E+05, 0.27842E+05, 0.27919E+05, 0.27995E+05, 0.28072E+05, 0.28149E+05, 0.28226E+05, 0.28303E+05 &
, 0.28380E+05, 0.28457E+05, 0.28535E+05, 0.28612E+05, 0.28690E+05, 0.28768E+05, 0.28845E+05, 0.28923E+05, 0.29001E+05, 0.29079E+05 &
, 0.29157E+05, 0.29236E+05, 0.29314E+05, 0.29392E+05, 0.29471E+05, 0.29549E+05, 0.29628E+05, 0.29707E+05, 0.29786E+05, 0.29865E+05 &
, 0.29944E+05, 0.30023E+05, 0.30103E+05, 0.30182E+05, 0.30262E+05, 0.30341E+05, 0.30421E+05, 0.30501E+05, 0.30581E+05, 0.30661E+05 &
, 0.30741E+05, 0.30821E+05, 0.30901E+05, 0.30982E+05, 0.31062E+05, 0.31143E+05, 0.31224E+05, 0.31305E+05, 0.31385E+05, 0.31466E+05 &
, 0.31548E+05, 0.31629E+05, 0.31710E+05, 0.31792E+05, 0.31873E+05, 0.31955E+05, 0.32036E+05, 0.32118E+05, 0.32200E+05, 0.32282E+05 /
 !            8           1
 DATA(QofT(          65 ,J),J=1,501)/   54.3889999389648       &
, 0.59297E+02, 0.64217E+02, 0.69156E+02, 0.74117E+02, 0.79106E+02, 0.84130E+02, 0.89193E+02, 0.94302E+02, 0.99461E+02, 0.10468E+03 &
, 0.10995E+03, 0.11529E+03, 0.12069E+03, 0.12617E+03, 0.13171E+03, 0.13733E+03, 0.14302E+03, 0.14879E+03, 0.15463E+03, 0.16055E+03 &
, 0.16654E+03, 0.17261E+03, 0.17876E+03, 0.18498E+03, 0.19127E+03, 0.19763E+03, 0.20407E+03, 0.21057E+03, 0.21714E+03, 0.22378E+03 &
, 0.23049E+03, 0.23726E+03, 0.24410E+03, 0.25100E+03, 0.25796E+03, 0.26497E+03, 0.27205E+03, 0.27918E+03, 0.28637E+03, 0.29361E+03 &
, 0.30090E+03, 0.30825E+03, 0.31564E+03, 0.32308E+03, 0.33080E+03, 0.33843E+03, 0.34611E+03, 0.35382E+03, 0.36158E+03, 0.36938E+03 &
, 0.37722E+03, 0.38510E+03, 0.39301E+03, 0.40096E+03, 0.40895E+03, 0.41697E+03, 0.42503E+03, 0.43312E+03, 0.44124E+03, 0.44940E+03 &
, 0.45758E+03, 0.46580E+03, 0.47404E+03, 0.48231E+03, 0.49061E+03, 0.49894E+03, 0.50730E+03, 0.51568E+03, 0.52408E+03, 0.53251E+03 &
, 0.54096E+03, 0.54944E+03, 0.55794E+03, 0.56646E+03, 0.57500E+03, 0.58357E+03, 0.59215E+03, 0.60076E+03, 0.60938E+03, 0.61803E+03 &
, 0.62669E+03, 0.63537E+03, 0.64407E+03, 0.65279E+03, 0.66153E+03, 0.67028E+03, 0.67905E+03, 0.68783E+03, 0.69663E+03, 0.70545E+03 &
, 0.71428E+03, 0.72313E+03, 0.73199E+03, 0.74086E+03, 0.74975E+03, 0.75866E+03, 0.76757E+03, 0.77650E+03, 0.78544E+03, 0.79440E+03 &
, 0.80337E+03, 0.81235E+03, 0.82134E+03, 0.83035E+03, 0.83936E+03, 0.84839E+03, 0.85742E+03, 0.86647E+03, 0.87553E+03, 0.88460E+03 &
, 0.89367E+03, 0.90277E+03, 0.91187E+03, 0.92098E+03, 0.93009E+03, 0.93922E+03, 0.94835E+03, 0.95751E+03, 0.96666E+03, 0.97582E+03 &
, 0.98499E+03, 0.99418E+03, 0.10034E+04, 0.10126E+04, 0.10218E+04, 0.10310E+04, 0.10402E+04, 0.10494E+04, 0.10587E+04, 0.10679E+04 &
, 0.10772E+04, 0.10864E+04, 0.10957E+04, 0.11050E+04, 0.11142E+04, 0.11235E+04, 0.11328E+04, 0.11421E+04, 0.11514E+04, 0.11608E+04 &
, 0.11701E+04, 0.11794E+04, 0.11887E+04, 0.11981E+04, 0.12074E+04, 0.12168E+04, 0.12261E+04, 0.12355E+04, 0.12449E+04, 0.12543E+04 &
, 0.12636E+04, 0.12730E+04, 0.12824E+04, 0.12918E+04, 0.13012E+04, 0.13106E+04, 0.13201E+04, 0.13295E+04, 0.13389E+04, 0.13484E+04 &
, 0.13578E+04, 0.13672E+04, 0.13767E+04, 0.13861E+04, 0.13956E+04, 0.14051E+04, 0.14145E+04, 0.14240E+04, 0.14335E+04, 0.14430E+04 &
, 0.14525E+04, 0.14620E+04, 0.14715E+04, 0.14810E+04, 0.14905E+04, 0.15001E+04, 0.15096E+04, 0.15191E+04, 0.15286E+04, 0.15382E+04 &
, 0.15477E+04, 0.15573E+04, 0.15668E+04, 0.15764E+04, 0.15860E+04, 0.15955E+04, 0.16051E+04, 0.16147E+04, 0.16243E+04, 0.16339E+04 &
, 0.16435E+04, 0.16531E+04, 0.16627E+04, 0.16723E+04, 0.16820E+04, 0.16916E+04, 0.17012E+04, 0.17108E+04, 0.17205E+04, 0.17301E+04 &
, 0.17398E+04, 0.17494E+04, 0.17591E+04, 0.17688E+04, 0.17785E+04, 0.17881E+04, 0.17978E+04, 0.18075E+04, 0.18172E+04, 0.18269E+04 &
, 0.18366E+04, 0.18463E+04, 0.18561E+04, 0.18658E+04, 0.18755E+04, 0.18853E+04, 0.18950E+04, 0.19048E+04, 0.19145E+04, 0.19243E+04 &
, 0.19340E+04, 0.19438E+04, 0.19536E+04, 0.19634E+04, 0.19732E+04, 0.19830E+04, 0.19928E+04, 0.20026E+04, 0.20124E+04, 0.20222E+04 &
, 0.20320E+04, 0.20418E+04, 0.20517E+04, 0.20615E+04, 0.20714E+04, 0.20812E+04, 0.20911E+04, 0.21010E+04, 0.21109E+04, 0.21207E+04 &
, 0.21306E+04, 0.21405E+04, 0.21504E+04, 0.21603E+04, 0.21703E+04, 0.21802E+04, 0.21901E+04, 0.22000E+04, 0.22100E+04, 0.22199E+04 &
, 0.22299E+04, 0.22399E+04, 0.22498E+04, 0.22598E+04, 0.22698E+04, 0.22798E+04, 0.22898E+04, 0.22998E+04, 0.23098E+04, 0.23198E+04 &
, 0.23298E+04, 0.23399E+04, 0.23499E+04, 0.23600E+04, 0.23700E+04, 0.23801E+04, 0.23901E+04, 0.24002E+04, 0.24103E+04, 0.24204E+04 &
, 0.24305E+04, 0.24406E+04, 0.24507E+04, 0.24609E+04, 0.24710E+04, 0.24811E+04, 0.24913E+04, 0.25014E+04, 0.25116E+04, 0.25217E+04 &
, 0.25319E+04, 0.25421E+04, 0.25523E+04, 0.25625E+04, 0.25727E+04, 0.25830E+04, 0.25932E+04, 0.26034E+04, 0.26137E+04, 0.26239E+04 &
, 0.26342E+04, 0.26444E+04, 0.26547E+04, 0.26650E+04, 0.26753E+04, 0.26856E+04, 0.26959E+04, 0.27062E+04, 0.27165E+04, 0.27269E+04 &
, 0.27372E+04, 0.27476E+04, 0.27579E+04, 0.27683E+04, 0.27787E+04, 0.27891E+04, 0.27994E+04, 0.28099E+04, 0.28203E+04, 0.28307E+04 &
, 0.28411E+04, 0.28516E+04, 0.28620E+04, 0.28725E+04, 0.28829E+04, 0.28934E+04, 0.29039E+04, 0.29144E+04, 0.29249E+04, 0.29354E+04 &
, 0.29459E+04, 0.29564E+04, 0.29670E+04, 0.29775E+04, 0.29881E+04, 0.29987E+04, 0.30092E+04, 0.30199E+04, 0.30304E+04, 0.30410E+04 &
, 0.30517E+04, 0.30623E+04, 0.30729E+04, 0.30836E+04, 0.30942E+04, 0.31049E+04, 0.31156E+04, 0.31262E+04, 0.31369E+04, 0.31476E+04 &
, 0.31584E+04, 0.31691E+04, 0.31798E+04, 0.31906E+04, 0.32013E+04, 0.32121E+04, 0.32228E+04, 0.32336E+04, 0.32444E+04, 0.32552E+04 &
, 0.32660E+04, 0.32768E+04, 0.32877E+04, 0.32985E+04, 0.33094E+04, 0.33202E+04, 0.33311E+04, 0.33420E+04, 0.33529E+04, 0.33638E+04 &
, 0.33747E+04, 0.33856E+04, 0.33966E+04, 0.34075E+04, 0.34185E+04, 0.34295E+04, 0.34404E+04, 0.34514E+04, 0.34624E+04, 0.34734E+04 &
, 0.34845E+04, 0.34955E+04, 0.35065E+04, 0.35176E+04, 0.35286E+04, 0.35397E+04, 0.35508E+04, 0.35619E+04, 0.35730E+04, 0.35841E+04 &
, 0.35952E+04, 0.36064E+04, 0.36175E+04, 0.36287E+04, 0.36399E+04, 0.36511E+04, 0.36622E+04, 0.36734E+04, 0.36847E+04, 0.36959E+04 &
, 0.37071E+04, 0.37184E+04, 0.37296E+04, 0.37409E+04, 0.37522E+04, 0.37635E+04, 0.37748E+04, 0.37861E+04, 0.37974E+04, 0.38088E+04 &
, 0.38201E+04, 0.38315E+04, 0.38429E+04, 0.38543E+04, 0.38656E+04, 0.38770E+04, 0.38885E+04, 0.38999E+04, 0.39113E+04, 0.39228E+04 &
, 0.39343E+04, 0.39457E+04, 0.39572E+04, 0.39687E+04, 0.39802E+04, 0.39918E+04, 0.40033E+04, 0.40148E+04, 0.40264E+04, 0.40379E+04 &
, 0.40495E+04, 0.40611E+04, 0.40727E+04, 0.40843E+04, 0.40960E+04, 0.41076E+04, 0.41193E+04, 0.41309E+04, 0.41426E+04, 0.41543E+04 &
, 0.41660E+04, 0.41777E+04, 0.41894E+04, 0.42011E+04, 0.42129E+04, 0.42247E+04, 0.42364E+04, 0.42482E+04, 0.42600E+04, 0.42718E+04 &
, 0.42836E+04, 0.42955E+04, 0.43073E+04, 0.43191E+04, 0.43310E+04, 0.43429E+04, 0.43548E+04, 0.43667E+04, 0.43786E+04, 0.43905E+04 &
, 0.44024E+04, 0.44144E+04, 0.44264E+04, 0.44384E+04, 0.44503E+04, 0.44623E+04, 0.44743E+04, 0.44864E+04, 0.44984E+04, 0.45105E+04 &
, 0.45225E+04, 0.45346E+04, 0.45467E+04, 0.45588E+04, 0.45709E+04, 0.45830E+04, 0.45952E+04, 0.46073E+04, 0.46194E+04, 0.46316E+04 &
, 0.46438E+04, 0.46560E+04, 0.46682E+04, 0.46804E+04, 0.46927E+04, 0.47049E+04, 0.47172E+04, 0.47295E+04, 0.47418E+04, 0.47540E+04 &
, 0.47663E+04, 0.47787E+04, 0.47910E+04, 0.48033E+04, 0.48157E+04, 0.48281E+04, 0.48405E+04, 0.48529E+04, 0.48653E+04, 0.48777E+04 &
, 0.48901E+04, 0.49026E+04, 0.49151E+04, 0.49275E+04, 0.49400E+04, 0.49525E+04, 0.49651E+04, 0.49776E+04, 0.49901E+04, 0.50026E+04 /
 !            8           2
 DATA(QofT(          66 ,J),J=1,501)/   37.4589996337891       &
, 0.40852E+02, 0.44253E+02, 0.47667E+02, 0.51096E+02, 0.54544E+02, 0.58016E+02, 0.61516E+02, 0.65047E+02, 0.68612E+02, 0.72216E+02 &
, 0.75861E+02, 0.79550E+02, 0.83285E+02, 0.87067E+02, 0.90899E+02, 0.94782E+02, 0.98715E+02, 0.10270E+03, 0.10674E+03, 0.11083E+03 &
, 0.11497E+03, 0.11916E+03, 0.12341E+03, 0.12771E+03, 0.13205E+03, 0.13645E+03, 0.14090E+03, 0.14539E+03, 0.14994E+03, 0.15453E+03 &
, 0.15916E+03, 0.16384E+03, 0.16856E+03, 0.17333E+03, 0.17814E+03, 0.18299E+03, 0.18788E+03, 0.19280E+03, 0.19777E+03, 0.20277E+03 &
, 0.20781E+03, 0.21289E+03, 0.21811E+03, 0.22332E+03, 0.22856E+03, 0.23383E+03, 0.23913E+03, 0.24446E+03, 0.24983E+03, 0.25521E+03 &
, 0.26063E+03, 0.26608E+03, 0.27155E+03, 0.27704E+03, 0.28256E+03, 0.28811E+03, 0.29367E+03, 0.29927E+03, 0.30488E+03, 0.31051E+03 &
, 0.31617E+03, 0.32185E+03, 0.32754E+03, 0.33326E+03, 0.33900E+03, 0.34475E+03, 0.35052E+03, 0.35631E+03, 0.36212E+03, 0.36795E+03 &
, 0.37379E+03, 0.37965E+03, 0.38552E+03, 0.39141E+03, 0.39731E+03, 0.40323E+03, 0.40917E+03, 0.41511E+03, 0.42107E+03, 0.42705E+03 &
, 0.43304E+03, 0.43904E+03, 0.44505E+03, 0.45107E+03, 0.45711E+03, 0.46316E+03, 0.46922E+03, 0.47529E+03, 0.48137E+03, 0.48746E+03 &
, 0.49357E+03, 0.49968E+03, 0.50580E+03, 0.51194E+03, 0.51808E+03, 0.52423E+03, 0.53039E+03, 0.53657E+03, 0.54275E+03, 0.54894E+03 &
, 0.55514E+03, 0.56134E+03, 0.56756E+03, 0.57378E+03, 0.58001E+03, 0.58624E+03, 0.59249E+03, 0.59874E+03, 0.60500E+03, 0.61127E+03 &
, 0.61755E+03, 0.62383E+03, 0.63011E+03, 0.63641E+03, 0.64271E+03, 0.64902E+03, 0.65534E+03, 0.66165E+03, 0.66798E+03, 0.67432E+03 &
, 0.68066E+03, 0.68700E+03, 0.69335E+03, 0.69971E+03, 0.70607E+03, 0.71244E+03, 0.71881E+03, 0.72519E+03, 0.73157E+03, 0.73797E+03 &
, 0.74436E+03, 0.75076E+03, 0.75716E+03, 0.76357E+03, 0.76999E+03, 0.77640E+03, 0.78283E+03, 0.78926E+03, 0.79570E+03, 0.80214E+03 &
, 0.80857E+03, 0.81502E+03, 0.82147E+03, 0.82793E+03, 0.83439E+03, 0.84085E+03, 0.84733E+03, 0.85380E+03, 0.86028E+03, 0.86676E+03 &
, 0.87325E+03, 0.87974E+03, 0.88623E+03, 0.89273E+03, 0.89923E+03, 0.90574E+03, 0.91225E+03, 0.91876E+03, 0.92528E+03, 0.93181E+03 &
, 0.93833E+03, 0.94487E+03, 0.95140E+03, 0.95794E+03, 0.96448E+03, 0.97103E+03, 0.97758E+03, 0.98413E+03, 0.99069E+03, 0.99725E+03 &
, 0.10038E+04, 0.10104E+04, 0.10170E+04, 0.10235E+04, 0.10301E+04, 0.10367E+04, 0.10433E+04, 0.10499E+04, 0.10565E+04, 0.10631E+04 &
, 0.10697E+04, 0.10763E+04, 0.10829E+04, 0.10895E+04, 0.10961E+04, 0.11027E+04, 0.11094E+04, 0.11160E+04, 0.11226E+04, 0.11293E+04 &
, 0.11359E+04, 0.11425E+04, 0.11492E+04, 0.11559E+04, 0.11625E+04, 0.11692E+04, 0.11758E+04, 0.11825E+04, 0.11892E+04, 0.11958E+04 &
, 0.12025E+04, 0.12092E+04, 0.12159E+04, 0.12226E+04, 0.12293E+04, 0.12360E+04, 0.12427E+04, 0.12494E+04, 0.12561E+04, 0.12628E+04 &
, 0.12695E+04, 0.12763E+04, 0.12830E+04, 0.12897E+04, 0.12964E+04, 0.13032E+04, 0.13099E+04, 0.13167E+04, 0.13234E+04, 0.13302E+04 &
, 0.13369E+04, 0.13437E+04, 0.13505E+04, 0.13572E+04, 0.13640E+04, 0.13708E+04, 0.13776E+04, 0.13844E+04, 0.13912E+04, 0.13980E+04 &
, 0.14048E+04, 0.14116E+04, 0.14184E+04, 0.14252E+04, 0.14320E+04, 0.14389E+04, 0.14457E+04, 0.14525E+04, 0.14594E+04, 0.14662E+04 &
, 0.14731E+04, 0.14799E+04, 0.14868E+04, 0.14936E+04, 0.15005E+04, 0.15074E+04, 0.15143E+04, 0.15211E+04, 0.15280E+04, 0.15349E+04 &
, 0.15418E+04, 0.15487E+04, 0.15556E+04, 0.15625E+04, 0.15695E+04, 0.15764E+04, 0.15833E+04, 0.15902E+04, 0.15972E+04, 0.16041E+04 &
, 0.16111E+04, 0.16180E+04, 0.16250E+04, 0.16320E+04, 0.16389E+04, 0.16459E+04, 0.16529E+04, 0.16599E+04, 0.16669E+04, 0.16739E+04 &
, 0.16809E+04, 0.16879E+04, 0.16949E+04, 0.17019E+04, 0.17089E+04, 0.17159E+04, 0.17230E+04, 0.17300E+04, 0.17371E+04, 0.17441E+04 &
, 0.17512E+04, 0.17582E+04, 0.17653E+04, 0.17724E+04, 0.17795E+04, 0.17866E+04, 0.17937E+04, 0.18008E+04, 0.18079E+04, 0.18150E+04 &
, 0.18221E+04, 0.18292E+04, 0.18363E+04, 0.18435E+04, 0.18506E+04, 0.18578E+04, 0.18649E+04, 0.18721E+04, 0.18792E+04, 0.18864E+04 &
, 0.18936E+04, 0.19008E+04, 0.19079E+04, 0.19151E+04, 0.19223E+04, 0.19295E+04, 0.19368E+04, 0.19440E+04, 0.19512E+04, 0.19584E+04 &
, 0.19657E+04, 0.19729E+04, 0.19802E+04, 0.19874E+04, 0.19947E+04, 0.20020E+04, 0.20093E+04, 0.20166E+04, 0.20238E+04, 0.20311E+04 &
, 0.20385E+04, 0.20458E+04, 0.20531E+04, 0.20604E+04, 0.20677E+04, 0.20751E+04, 0.20824E+04, 0.20898E+04, 0.20971E+04, 0.21045E+04 &
, 0.21119E+04, 0.21192E+04, 0.21266E+04, 0.21340E+04, 0.21414E+04, 0.21488E+04, 0.21563E+04, 0.21637E+04, 0.21711E+04, 0.21785E+04 &
, 0.21860E+04, 0.21934E+04, 0.22009E+04, 0.22084E+04, 0.22158E+04, 0.22233E+04, 0.22308E+04, 0.22383E+04, 0.22458E+04, 0.22533E+04 &
, 0.22608E+04, 0.22683E+04, 0.22759E+04, 0.22834E+04, 0.22909E+04, 0.22985E+04, 0.23060E+04, 0.23136E+04, 0.23212E+04, 0.23288E+04 &
, 0.23363E+04, 0.23439E+04, 0.23515E+04, 0.23592E+04, 0.23668E+04, 0.23744E+04, 0.23820E+04, 0.23897E+04, 0.23973E+04, 0.24050E+04 &
, 0.24126E+04, 0.24203E+04, 0.24280E+04, 0.24357E+04, 0.24434E+04, 0.24511E+04, 0.24588E+04, 0.24665E+04, 0.24742E+04, 0.24819E+04 &
, 0.24897E+04, 0.24974E+04, 0.25052E+04, 0.25129E+04, 0.25207E+04, 0.25285E+04, 0.25363E+04, 0.25441E+04, 0.25519E+04, 0.25597E+04 &
, 0.25675E+04, 0.25753E+04, 0.25832E+04, 0.25910E+04, 0.25989E+04, 0.26067E+04, 0.26146E+04, 0.26225E+04, 0.26303E+04, 0.26382E+04 &
, 0.26461E+04, 0.26540E+04, 0.26619E+04, 0.26699E+04, 0.26778E+04, 0.26857E+04, 0.26937E+04, 0.27016E+04, 0.27096E+04, 0.27176E+04 &
, 0.27255E+04, 0.27335E+04, 0.27415E+04, 0.27495E+04, 0.27575E+04, 0.27656E+04, 0.27736E+04, 0.27816E+04, 0.27897E+04, 0.27977E+04 &
, 0.28058E+04, 0.28139E+04, 0.28219E+04, 0.28300E+04, 0.28381E+04, 0.28462E+04, 0.28543E+04, 0.28625E+04, 0.28706E+04, 0.28787E+04 &
, 0.28869E+04, 0.28950E+04, 0.29032E+04, 0.29114E+04, 0.29195E+04, 0.29277E+04, 0.29359E+04, 0.29441E+04, 0.29523E+04, 0.29606E+04 &
, 0.29688E+04, 0.29770E+04, 0.29853E+04, 0.29935E+04, 0.30018E+04, 0.30101E+04, 0.30184E+04, 0.30267E+04, 0.30349E+04, 0.30432E+04 &
, 0.30516E+04, 0.30599E+04, 0.30682E+04, 0.30766E+04, 0.30849E+04, 0.30933E+04, 0.31017E+04, 0.31100E+04, 0.31184E+04, 0.31268E+04 &
, 0.31352E+04, 0.31436E+04, 0.31520E+04, 0.31605E+04, 0.31689E+04, 0.31774E+04, 0.31858E+04, 0.31943E+04, 0.32027E+04, 0.32112E+04 &
, 0.32197E+04, 0.32282E+04, 0.32367E+04, 0.32452E+04, 0.32538E+04, 0.32623E+04, 0.32709E+04, 0.32794E+04, 0.32880E+04, 0.32966E+04 &
, 0.33051E+04, 0.33137E+04, 0.33223E+04, 0.33309E+04, 0.33395E+04, 0.33482E+04, 0.33568E+04, 0.33654E+04, 0.33741E+04, 0.33827E+04 &
, 0.33914E+04, 0.34001E+04, 0.34088E+04, 0.34175E+04, 0.34262E+04, 0.34349E+04, 0.34436E+04, 0.34524E+04, 0.34611E+04, 0.34699E+04 /
 !            8           3
 DATA(QofT(          67 ,J),J=1,501)/   57.0660018920898       &
, 0.62244E+02, 0.67435E+02, 0.72645E+02, 0.77879E+02, 0.83142E+02, 0.88441E+02, 0.93782E+02, 0.99170E+02, 0.10461E+03, 0.11011E+03 &
, 0.11567E+03, 0.12130E+03, 0.12700E+03, 0.13278E+03, 0.13862E+03, 0.14455E+03, 0.15055E+03, 0.15663E+03, 0.16279E+03, 0.16903E+03 &
, 0.17535E+03, 0.18175E+03, 0.18823E+03, 0.19479E+03, 0.20142E+03, 0.20813E+03, 0.21492E+03, 0.22178E+03, 0.22871E+03, 0.23571E+03 &
, 0.24279E+03, 0.24993E+03, 0.25713E+03, 0.26441E+03, 0.27175E+03, 0.27914E+03, 0.28661E+03, 0.29413E+03, 0.30170E+03, 0.30934E+03 &
, 0.31703E+03, 0.32477E+03, 0.33257E+03, 0.34076E+03, 0.34876E+03, 0.35680E+03, 0.36489E+03, 0.37303E+03, 0.38121E+03, 0.38944E+03 &
, 0.39771E+03, 0.40601E+03, 0.41436E+03, 0.42275E+03, 0.43117E+03, 0.43963E+03, 0.44813E+03, 0.45666E+03, 0.46523E+03, 0.47383E+03 &
, 0.48246E+03, 0.49112E+03, 0.49982E+03, 0.50854E+03, 0.51729E+03, 0.52608E+03, 0.53489E+03, 0.54372E+03, 0.55259E+03, 0.56148E+03 &
, 0.57039E+03, 0.57933E+03, 0.58829E+03, 0.59728E+03, 0.60629E+03, 0.61532E+03, 0.62438E+03, 0.63345E+03, 0.64255E+03, 0.65167E+03 &
, 0.66080E+03, 0.66996E+03, 0.67914E+03, 0.68833E+03, 0.69754E+03, 0.70677E+03, 0.71602E+03, 0.72529E+03, 0.73457E+03, 0.74387E+03 &
, 0.75318E+03, 0.76251E+03, 0.77185E+03, 0.78121E+03, 0.79059E+03, 0.79998E+03, 0.80938E+03, 0.81881E+03, 0.82824E+03, 0.83768E+03 &
, 0.84714E+03, 0.85661E+03, 0.86609E+03, 0.87559E+03, 0.88509E+03, 0.89461E+03, 0.90414E+03, 0.91368E+03, 0.92324E+03, 0.93281E+03 &
, 0.94238E+03, 0.95197E+03, 0.96156E+03, 0.97117E+03, 0.98079E+03, 0.99042E+03, 0.10000E+04, 0.10097E+04, 0.10194E+04, 0.10290E+04 &
, 0.10387E+04, 0.10484E+04, 0.10581E+04, 0.10678E+04, 0.10775E+04, 0.10872E+04, 0.10969E+04, 0.11067E+04, 0.11164E+04, 0.11262E+04 &
, 0.11359E+04, 0.11457E+04, 0.11555E+04, 0.11652E+04, 0.11750E+04, 0.11848E+04, 0.11946E+04, 0.12044E+04, 0.12143E+04, 0.12241E+04 &
, 0.12339E+04, 0.12438E+04, 0.12536E+04, 0.12635E+04, 0.12733E+04, 0.12832E+04, 0.12931E+04, 0.13029E+04, 0.13128E+04, 0.13227E+04 &
, 0.13326E+04, 0.13425E+04, 0.13524E+04, 0.13624E+04, 0.13723E+04, 0.13822E+04, 0.13922E+04, 0.14021E+04, 0.14121E+04, 0.14220E+04 &
, 0.14320E+04, 0.14419E+04, 0.14519E+04, 0.14619E+04, 0.14719E+04, 0.14819E+04, 0.14919E+04, 0.15019E+04, 0.15119E+04, 0.15219E+04 &
, 0.15319E+04, 0.15419E+04, 0.15520E+04, 0.15620E+04, 0.15721E+04, 0.15821E+04, 0.15922E+04, 0.16022E+04, 0.16123E+04, 0.16224E+04 &
, 0.16325E+04, 0.16426E+04, 0.16526E+04, 0.16627E+04, 0.16728E+04, 0.16830E+04, 0.16931E+04, 0.17032E+04, 0.17133E+04, 0.17234E+04 &
, 0.17336E+04, 0.17437E+04, 0.17539E+04, 0.17640E+04, 0.17742E+04, 0.17844E+04, 0.17945E+04, 0.18047E+04, 0.18149E+04, 0.18251E+04 &
, 0.18353E+04, 0.18455E+04, 0.18557E+04, 0.18659E+04, 0.18761E+04, 0.18864E+04, 0.18966E+04, 0.19068E+04, 0.19171E+04, 0.19273E+04 &
, 0.19376E+04, 0.19479E+04, 0.19581E+04, 0.19684E+04, 0.19787E+04, 0.19890E+04, 0.19993E+04, 0.20096E+04, 0.20199E+04, 0.20302E+04 &
, 0.20406E+04, 0.20509E+04, 0.20612E+04, 0.20716E+04, 0.20819E+04, 0.20923E+04, 0.21027E+04, 0.21130E+04, 0.21234E+04, 0.21338E+04 &
, 0.21442E+04, 0.21546E+04, 0.21650E+04, 0.21754E+04, 0.21858E+04, 0.21962E+04, 0.22067E+04, 0.22171E+04, 0.22276E+04, 0.22380E+04 &
, 0.22485E+04, 0.22590E+04, 0.22694E+04, 0.22799E+04, 0.22904E+04, 0.23009E+04, 0.23114E+04, 0.23220E+04, 0.23325E+04, 0.23430E+04 &
, 0.23535E+04, 0.23641E+04, 0.23747E+04, 0.23852E+04, 0.23958E+04, 0.24064E+04, 0.24169E+04, 0.24276E+04, 0.24382E+04, 0.24488E+04 &
, 0.24594E+04, 0.24700E+04, 0.24807E+04, 0.24913E+04, 0.25019E+04, 0.25126E+04, 0.25233E+04, 0.25339E+04, 0.25446E+04, 0.25553E+04 &
, 0.25660E+04, 0.25767E+04, 0.25875E+04, 0.25982E+04, 0.26089E+04, 0.26197E+04, 0.26304E+04, 0.26412E+04, 0.26519E+04, 0.26627E+04 &
, 0.26735E+04, 0.26843E+04, 0.26951E+04, 0.27059E+04, 0.27168E+04, 0.27276E+04, 0.27384E+04, 0.27493E+04, 0.27602E+04, 0.27710E+04 &
, 0.27819E+04, 0.27928E+04, 0.28037E+04, 0.28146E+04, 0.28255E+04, 0.28364E+04, 0.28474E+04, 0.28583E+04, 0.28693E+04, 0.28803E+04 &
, 0.28912E+04, 0.29022E+04, 0.29132E+04, 0.29242E+04, 0.29352E+04, 0.29463E+04, 0.29573E+04, 0.29683E+04, 0.29794E+04, 0.29904E+04 &
, 0.30015E+04, 0.30126E+04, 0.30237E+04, 0.30348E+04, 0.30459E+04, 0.30570E+04, 0.30682E+04, 0.30793E+04, 0.30904E+04, 0.31016E+04 &
, 0.31128E+04, 0.31240E+04, 0.31352E+04, 0.31464E+04, 0.31576E+04, 0.31688E+04, 0.31801E+04, 0.31913E+04, 0.32026E+04, 0.32138E+04 &
, 0.32251E+04, 0.32364E+04, 0.32477E+04, 0.32590E+04, 0.32703E+04, 0.32817E+04, 0.32930E+04, 0.33044E+04, 0.33157E+04, 0.33271E+04 &
, 0.33385E+04, 0.33499E+04, 0.33613E+04, 0.33727E+04, 0.33842E+04, 0.33956E+04, 0.34070E+04, 0.34185E+04, 0.34300E+04, 0.34415E+04 &
, 0.34530E+04, 0.34645E+04, 0.34760E+04, 0.34875E+04, 0.34991E+04, 0.35106E+04, 0.35222E+04, 0.35338E+04, 0.35454E+04, 0.35570E+04 &
, 0.35686E+04, 0.35802E+04, 0.35918E+04, 0.36035E+04, 0.36151E+04, 0.36268E+04, 0.36385E+04, 0.36502E+04, 0.36619E+04, 0.36736E+04 &
, 0.36853E+04, 0.36971E+04, 0.37088E+04, 0.37206E+04, 0.37324E+04, 0.37441E+04, 0.37560E+04, 0.37678E+04, 0.37796E+04, 0.37914E+04 &
, 0.38033E+04, 0.38151E+04, 0.38270E+04, 0.38389E+04, 0.38508E+04, 0.38627E+04, 0.38746E+04, 0.38866E+04, 0.38985E+04, 0.39105E+04 &
, 0.39224E+04, 0.39344E+04, 0.39464E+04, 0.39584E+04, 0.39704E+04, 0.39824E+04, 0.39945E+04, 0.40065E+04, 0.40186E+04, 0.40307E+04 &
, 0.40428E+04, 0.40549E+04, 0.40670E+04, 0.40792E+04, 0.40913E+04, 0.41035E+04, 0.41156E+04, 0.41278E+04, 0.41400E+04, 0.41522E+04 &
, 0.41644E+04, 0.41767E+04, 0.41889E+04, 0.42012E+04, 0.42134E+04, 0.42257E+04, 0.42380E+04, 0.42503E+04, 0.42626E+04, 0.42750E+04 &
, 0.42873E+04, 0.42997E+04, 0.43121E+04, 0.43244E+04, 0.43368E+04, 0.43492E+04, 0.43616E+04, 0.43741E+04, 0.43866E+04, 0.43990E+04 &
, 0.44115E+04, 0.44240E+04, 0.44365E+04, 0.44490E+04, 0.44615E+04, 0.44741E+04, 0.44866E+04, 0.44992E+04, 0.45118E+04, 0.45244E+04 &
, 0.45370E+04, 0.45496E+04, 0.45622E+04, 0.45749E+04, 0.45876E+04, 0.46002E+04, 0.46129E+04, 0.46256E+04, 0.46383E+04, 0.46510E+04 &
, 0.46638E+04, 0.46766E+04, 0.46893E+04, 0.47021E+04, 0.47149E+04, 0.47277E+04, 0.47405E+04, 0.47534E+04, 0.47662E+04, 0.47791E+04 &
, 0.47919E+04, 0.48048E+04, 0.48177E+04, 0.48306E+04, 0.48436E+04, 0.48565E+04, 0.48695E+04, 0.48825E+04, 0.48954E+04, 0.49084E+04 &
, 0.49214E+04, 0.49344E+04, 0.49475E+04, 0.49605E+04, 0.49736E+04, 0.49867E+04, 0.49998E+04, 0.50129E+04, 0.50260E+04, 0.50391E+04 &
, 0.50523E+04, 0.50655E+04, 0.50786E+04, 0.50918E+04, 0.51050E+04, 0.51183E+04, 0.51315E+04, 0.51447E+04, 0.51580E+04, 0.51712E+04 &
, 0.51846E+04, 0.51979E+04, 0.52112E+04, 0.52245E+04, 0.52378E+04, 0.52512E+04, 0.52645E+04, 0.52779E+04, 0.52914E+04, 0.53048E+04 /
 !            9           1
 DATA(QofT(          68 ,J),J=1,501)/   102.230003356934       &
, 0.11787E+03, 0.13424E+03, 0.15129E+03, 0.16902E+03, 0.18739E+03, 0.20638E+03, 0.22598E+03, 0.24616E+03, 0.26690E+03, 0.28820E+03 &
, 0.31004E+03, 0.33240E+03, 0.35528E+03, 0.37866E+03, 0.40253E+03, 0.42689E+03, 0.45172E+03, 0.47701E+03, 0.50276E+03, 0.52896E+03 &
, 0.55560E+03, 0.58267E+03, 0.61017E+03, 0.63812E+03, 0.66646E+03, 0.69520E+03, 0.72436E+03, 0.75391E+03, 0.78385E+03, 0.81419E+03 &
, 0.84491E+03, 0.87601E+03, 0.90750E+03, 0.93936E+03, 0.97159E+03, 0.10042E+04, 0.10372E+04, 0.10705E+04, 0.11042E+04, 0.11382E+04 &
, 0.11727E+04, 0.12074E+04, 0.12426E+04, 0.12780E+04, 0.13139E+04, 0.13501E+04, 0.13866E+04, 0.14235E+04, 0.14608E+04, 0.14984E+04 &
, 0.15363E+04, 0.15746E+04, 0.16133E+04, 0.16523E+04, 0.16916E+04, 0.17313E+04, 0.17714E+04, 0.18118E+04, 0.18526E+04, 0.18937E+04 &
, 0.19351E+04, 0.19770E+04, 0.20192E+04, 0.20617E+04, 0.21046E+04, 0.21479E+04, 0.21915E+04, 0.22355E+04, 0.22799E+04, 0.23246E+04 &
, 0.23698E+04, 0.24152E+04, 0.24611E+04, 0.25073E+04, 0.25539E+04, 0.26008E+04, 0.26482E+04, 0.26959E+04, 0.27440E+04, 0.27925E+04 &
, 0.28414E+04, 0.28907E+04, 0.29403E+04, 0.29904E+04, 0.30408E+04, 0.30917E+04, 0.31429E+04, 0.31946E+04, 0.32466E+04, 0.32991E+04 &
, 0.33519E+04, 0.34052E+04, 0.34589E+04, 0.35130E+04, 0.35675E+04, 0.36224E+04, 0.36777E+04, 0.37335E+04, 0.37897E+04, 0.38463E+04 &
, 0.39034E+04, 0.39608E+04, 0.40188E+04, 0.40771E+04, 0.41359E+04, 0.41951E+04, 0.42548E+04, 0.43149E+04, 0.43755E+04, 0.44365E+04 &
, 0.44981E+04, 0.45600E+04, 0.46224E+04, 0.46853E+04, 0.47486E+04, 0.48124E+04, 0.48767E+04, 0.49414E+04, 0.50066E+04, 0.50723E+04 &
, 0.51385E+04, 0.52052E+04, 0.52724E+04, 0.53400E+04, 0.54081E+04, 0.54768E+04, 0.55459E+04, 0.56155E+04, 0.56857E+04, 0.57563E+04 &
, 0.58275E+04, 0.58992E+04, 0.59714E+04, 0.60441E+04, 0.61174E+04, 0.61912E+04, 0.62654E+04, 0.63403E+04, 0.64157E+04, 0.64915E+04 &
, 0.65680E+04, 0.66450E+04, 0.67226E+04, 0.68007E+04, 0.68794E+04, 0.69585E+04, 0.70383E+04, 0.71187E+04, 0.71996E+04, 0.72811E+04 &
, 0.73632E+04, 0.74458E+04, 0.75290E+04, 0.76129E+04, 0.76972E+04, 0.77823E+04, 0.78678E+04, 0.79540E+04, 0.80408E+04, 0.81283E+04 &
, 0.82162E+04, 0.83049E+04, 0.83942E+04, 0.84840E+04, 0.85746E+04, 0.86657E+04, 0.87574E+04, 0.88498E+04, 0.89428E+04, 0.90366E+04 &
, 0.91309E+04, 0.92258E+04, 0.93215E+04, 0.94178E+04, 0.95147E+04, 0.96124E+04, 0.97107E+04, 0.98096E+04, 0.99093E+04, 0.10010E+05 &
, 0.10111E+05, 0.10212E+05, 0.10315E+05, 0.10418E+05, 0.10522E+05, 0.10626E+05, 0.10731E+05, 0.10837E+05, 0.10944E+05, 0.11051E+05 &
, 0.11159E+05, 0.11268E+05, 0.11378E+05, 0.11488E+05, 0.11599E+05, 0.11711E+05, 0.11824E+05, 0.11937E+05, 0.12051E+05, 0.12166E+05 &
, 0.12281E+05, 0.12398E+05, 0.12515E+05, 0.12633E+05, 0.12751E+05, 0.12871E+05, 0.12991E+05, 0.13112E+05, 0.13234E+05, 0.13357E+05 &
, 0.13480E+05, 0.13605E+05, 0.13730E+05, 0.13856E+05, 0.13982E+05, 0.14110E+05, 0.14239E+05, 0.14368E+05, 0.14498E+05, 0.14629E+05 &
, 0.14761E+05, 0.14894E+05, 0.15027E+05, 0.15162E+05, 0.15297E+05, 0.15433E+05, 0.15570E+05, 0.15708E+05, 0.15847E+05, 0.15987E+05 &
, 0.16127E+05, 0.16269E+05, 0.16411E+05, 0.16555E+05, 0.16699E+05, 0.16844E+05, 0.16990E+05, 0.17137E+05, 0.17285E+05, 0.17434E+05 &
, 0.17584E+05, 0.17735E+05, 0.17887E+05, 0.18040E+05, 0.18193E+05, 0.18348E+05, 0.18504E+05, 0.18661E+05, 0.18818E+05, 0.18977E+05 &
, 0.19137E+05, 0.19297E+05, 0.19459E+05, 0.19622E+05, 0.19785E+05, 0.19950E+05, 0.20116E+05, 0.20282E+05, 0.20450E+05, 0.20619E+05 &
, 0.20789E+05, 0.20960E+05, 0.21132E+05, 0.21305E+05, 0.21479E+05, 0.21654E+05, 0.21830E+05, 0.22008E+05, 0.22186E+05, 0.22366E+05 &
, 0.22546E+05, 0.22728E+05, 0.22911E+05, 0.23095E+05, 0.23280E+05, 0.23466E+05, 0.23654E+05, 0.23842E+05, 0.24032E+05, 0.24222E+05 &
, 0.24414E+05, 0.24607E+05, 0.24802E+05, 0.24997E+05, 0.25194E+05, 0.25391E+05, 0.25590E+05, 0.25790E+05, 0.25992E+05, 0.26194E+05 &
, 0.26398E+05, 0.26603E+05, 0.26809E+05, 0.27016E+05, 0.27225E+05, 0.27435E+05, 0.27646E+05, 0.27858E+05, 0.28072E+05, 0.28286E+05 &
, 0.28502E+05, 0.28720E+05, 0.28938E+05, 0.29158E+05, 0.29379E+05, 0.29602E+05, 0.29826E+05, 0.30051E+05, 0.30277E+05, 0.30505E+05 &
, 0.30734E+05, 0.30964E+05, 0.31196E+05, 0.31429E+05, 0.31663E+05, 0.31899E+05, 0.32136E+05, 0.32374E+05, 0.32614E+05, 0.32855E+05 &
, 0.33097E+05, 0.33341E+05, 0.33586E+05, 0.33833E+05, 0.34081E+05, 0.34331E+05, 0.34582E+05, 0.34834E+05, 0.35087E+05, 0.35343E+05 &
, 0.35599E+05, 0.35857E+05, 0.36117E+05, 0.36378E+05, 0.36640E+05, 0.36904E+05, 0.37170E+05, 0.37436E+05, 0.37705E+05, 0.37974E+05 &
, 0.38246E+05, 0.38519E+05, 0.38793E+05, 0.39069E+05, 0.39346E+05, 0.39625E+05, 0.39905E+05, 0.40187E+05, 0.40471E+05, 0.40756E+05 &
, 0.41043E+05, 0.41331E+05, 0.41621E+05, 0.41912E+05, 0.42205E+05, 0.42499E+05, 0.42796E+05, 0.43093E+05, 0.43393E+05, 0.43694E+05 &
, 0.43996E+05, 0.44301E+05, 0.44607E+05, 0.44914E+05, 0.45223E+05, 0.45534E+05, 0.45847E+05, 0.46161E+05, 0.46477E+05, 0.46794E+05 &
, 0.47113E+05, 0.47435E+05, 0.47757E+05, 0.48082E+05, 0.48408E+05, 0.48736E+05, 0.49065E+05, 0.49396E+05, 0.49729E+05, 0.50064E+05 &
, 0.50401E+05, 0.50739E+05, 0.51079E+05, 0.51421E+05, 0.51765E+05, 0.52110E+05, 0.52458E+05, 0.52807E+05, 0.53158E+05, 0.53510E+05 &
, 0.53865E+05, 0.54221E+05, 0.54580E+05, 0.54940E+05, 0.55302E+05, 0.55665E+05, 0.56031E+05, 0.56399E+05, 0.56768E+05, 0.57140E+05 &
, 0.57513E+05, 0.57888E+05, 0.58265E+05, 0.58644E+05, 0.59025E+05, 0.59408E+05, 0.59793E+05, 0.60180E+05, 0.60569E+05, 0.60959E+05 &
, 0.61352E+05, 0.61747E+05, 0.62143E+05, 0.62542E+05, 0.62943E+05, 0.63345E+05, 0.63750E+05, 0.64157E+05, 0.64566E+05, 0.64977E+05 &
, 0.65390E+05, 0.65805E+05, 0.66222E+05, 0.66641E+05, 0.67062E+05, 0.67485E+05, 0.67910E+05, 0.68338E+05, 0.68767E+05, 0.69199E+05 &
, 0.69633E+05, 0.70069E+05, 0.70507E+05, 0.70947E+05, 0.71390E+05, 0.71834E+05, 0.72281E+05, 0.72730E+05, 0.73181E+05, 0.73634E+05 &
, 0.74090E+05, 0.74547E+05, 0.75007E+05, 0.75469E+05, 0.75934E+05, 0.76400E+05, 0.76869E+05, 0.77341E+05, 0.77814E+05, 0.78290E+05 &
, 0.78768E+05, 0.79248E+05, 0.79731E+05, 0.80215E+05, 0.80702E+05, 0.81192E+05, 0.81684E+05, 0.82178E+05, 0.82675E+05, 0.83174E+05 &
, 0.83675E+05, 0.84178E+05, 0.84684E+05, 0.85193E+05, 0.85704E+05, 0.86217E+05, 0.86733E+05, 0.87251E+05, 0.87771E+05, 0.88294E+05 &
, 0.88820E+05, 0.89347E+05, 0.89878E+05, 0.90410E+05, 0.90946E+05, 0.91483E+05, 0.92024E+05, 0.92567E+05, 0.93112E+05, 0.93660E+05 &
, 0.94210E+05, 0.94763E+05, 0.95318E+05, 0.95876E+05, 0.96437E+05, 0.97000E+05, 0.97566E+05, 0.98134E+05, 0.98705E+05, 0.99279E+05 &
, 0.99855E+05, 0.10043E+06, 0.10101E+06, 0.10160E+06, 0.10219E+06, 0.10278E+06, 0.10337E+06, 0.10396E+06, 0.10456E+06, 0.10516E+06 /
 !            9           2
 DATA(QofT(          69 ,J),J=1,501)/   102.699996948242       &
, 0.11842E+03, 0.13486E+03, 0.15199E+03, 0.16980E+03, 0.18826E+03, 0.20734E+03, 0.22702E+03, 0.24730E+03, 0.26814E+03, 0.28954E+03 &
, 0.31147E+03, 0.33394E+03, 0.35693E+03, 0.38041E+03, 0.40440E+03, 0.42886E+03, 0.45381E+03, 0.47922E+03, 0.50508E+03, 0.53140E+03 &
, 0.55816E+03, 0.58536E+03, 0.61298E+03, 0.64104E+03, 0.66950E+03, 0.69838E+03, 0.72766E+03, 0.75735E+03, 0.78743E+03, 0.81790E+03 &
, 0.84876E+03, 0.88001E+03, 0.91164E+03, 0.94365E+03, 0.97602E+03, 0.10088E+04, 0.10419E+04, 0.10754E+04, 0.11092E+04, 0.11434E+04 &
, 0.11780E+04, 0.12129E+04, 0.12482E+04, 0.12839E+04, 0.13199E+04, 0.13562E+04, 0.13930E+04, 0.14300E+04, 0.14674E+04, 0.15052E+04 &
, 0.15433E+04, 0.15818E+04, 0.16206E+04, 0.16598E+04, 0.16993E+04, 0.17392E+04, 0.17795E+04, 0.18201E+04, 0.18610E+04, 0.19023E+04 &
, 0.19440E+04, 0.19860E+04, 0.20284E+04, 0.20712E+04, 0.21143E+04, 0.21577E+04, 0.22016E+04, 0.22457E+04, 0.22903E+04, 0.23353E+04 &
, 0.23806E+04, 0.24263E+04, 0.24723E+04, 0.25187E+04, 0.25655E+04, 0.26127E+04, 0.26603E+04, 0.27082E+04, 0.27566E+04, 0.28053E+04 &
, 0.28544E+04, 0.29039E+04, 0.29538E+04, 0.30040E+04, 0.30547E+04, 0.31058E+04, 0.31573E+04, 0.32092E+04, 0.32614E+04, 0.33141E+04 &
, 0.33672E+04, 0.34207E+04, 0.34746E+04, 0.35289E+04, 0.35837E+04, 0.36389E+04, 0.36945E+04, 0.37505E+04, 0.38069E+04, 0.38638E+04 &
, 0.39211E+04, 0.39789E+04, 0.40370E+04, 0.40957E+04, 0.41547E+04, 0.42142E+04, 0.42742E+04, 0.43346E+04, 0.43954E+04, 0.44567E+04 &
, 0.45185E+04, 0.45807E+04, 0.46434E+04, 0.47065E+04, 0.47702E+04, 0.48342E+04, 0.48988E+04, 0.49638E+04, 0.50294E+04, 0.50953E+04 &
, 0.51618E+04, 0.52288E+04, 0.52962E+04, 0.53642E+04, 0.54327E+04, 0.55016E+04, 0.55710E+04, 0.56410E+04, 0.57115E+04, 0.57824E+04 &
, 0.58540E+04, 0.59259E+04, 0.59984E+04, 0.60715E+04, 0.61450E+04, 0.62192E+04, 0.62938E+04, 0.63690E+04, 0.64447E+04, 0.65210E+04 &
, 0.65977E+04, 0.66751E+04, 0.67530E+04, 0.68314E+04, 0.69104E+04, 0.69900E+04, 0.70701E+04, 0.71509E+04, 0.72321E+04, 0.73140E+04 &
, 0.73964E+04, 0.74794E+04, 0.75630E+04, 0.76472E+04, 0.77320E+04, 0.78174E+04, 0.79034E+04, 0.79899E+04, 0.80771E+04, 0.81649E+04 &
, 0.82534E+04, 0.83424E+04, 0.84320E+04, 0.85223E+04, 0.86132E+04, 0.87047E+04, 0.87969E+04, 0.88897E+04, 0.89832E+04, 0.90772E+04 &
, 0.91720E+04, 0.92674E+04, 0.93635E+04, 0.94602E+04, 0.95576E+04, 0.96556E+04, 0.97544E+04, 0.98538E+04, 0.99539E+04, 0.10055E+05 &
, 0.10156E+05, 0.10258E+05, 0.10361E+05, 0.10465E+05, 0.10569E+05, 0.10674E+05, 0.10780E+05, 0.10886E+05, 0.10993E+05, 0.11101E+05 &
, 0.11210E+05, 0.11319E+05, 0.11429E+05, 0.11540E+05, 0.11651E+05, 0.11764E+05, 0.11877E+05, 0.11990E+05, 0.12105E+05, 0.12220E+05 &
, 0.12336E+05, 0.12453E+05, 0.12571E+05, 0.12689E+05, 0.12808E+05, 0.12929E+05, 0.13049E+05, 0.13171E+05, 0.13293E+05, 0.13417E+05 &
, 0.13541E+05, 0.13666E+05, 0.13791E+05, 0.13918E+05, 0.14045E+05, 0.14173E+05, 0.14302E+05, 0.14432E+05, 0.14563E+05, 0.14694E+05 &
, 0.14827E+05, 0.14960E+05, 0.15094E+05, 0.15229E+05, 0.15365E+05, 0.15502E+05, 0.15640E+05, 0.15778E+05, 0.15918E+05, 0.16058E+05 &
, 0.16199E+05, 0.16341E+05, 0.16484E+05, 0.16629E+05, 0.16773E+05, 0.16919E+05, 0.17066E+05, 0.17214E+05, 0.17362E+05, 0.17512E+05 &
, 0.17663E+05, 0.17814E+05, 0.17967E+05, 0.18120E+05, 0.18274E+05, 0.18430E+05, 0.18586E+05, 0.18744E+05, 0.18902E+05, 0.19061E+05 &
, 0.19222E+05, 0.19383E+05, 0.19545E+05, 0.19709E+05, 0.19873E+05, 0.20039E+05, 0.20205E+05, 0.20373E+05, 0.20541E+05, 0.20711E+05 &
, 0.20881E+05, 0.21053E+05, 0.21226E+05, 0.21399E+05, 0.21574E+05, 0.21750E+05, 0.21927E+05, 0.22106E+05, 0.22285E+05, 0.22465E+05 &
, 0.22646E+05, 0.22829E+05, 0.23013E+05, 0.23197E+05, 0.23383E+05, 0.23570E+05, 0.23758E+05, 0.23948E+05, 0.24138E+05, 0.24330E+05 &
, 0.24523E+05, 0.24717E+05, 0.24912E+05, 0.25108E+05, 0.25305E+05, 0.25504E+05, 0.25704E+05, 0.25905E+05, 0.26107E+05, 0.26310E+05 &
, 0.26515E+05, 0.26720E+05, 0.26928E+05, 0.27136E+05, 0.27345E+05, 0.27556E+05, 0.27768E+05, 0.27981E+05, 0.28196E+05, 0.28411E+05 &
, 0.28628E+05, 0.28847E+05, 0.29066E+05, 0.29287E+05, 0.29509E+05, 0.29733E+05, 0.29957E+05, 0.30183E+05, 0.30411E+05, 0.30639E+05 &
, 0.30869E+05, 0.31101E+05, 0.31333E+05, 0.31567E+05, 0.31803E+05, 0.32039E+05, 0.32277E+05, 0.32517E+05, 0.32758E+05, 0.33000E+05 &
, 0.33243E+05, 0.33488E+05, 0.33735E+05, 0.33982E+05, 0.34231E+05, 0.34482E+05, 0.34734E+05, 0.34987E+05, 0.35242E+05, 0.35498E+05 &
, 0.35756E+05, 0.36015E+05, 0.36276E+05, 0.36538E+05, 0.36801E+05, 0.37067E+05, 0.37333E+05, 0.37601E+05, 0.37871E+05, 0.38141E+05 &
, 0.38414E+05, 0.38688E+05, 0.38963E+05, 0.39241E+05, 0.39519E+05, 0.39799E+05, 0.40081E+05, 0.40364E+05, 0.40649E+05, 0.40935E+05 &
, 0.41223E+05, 0.41512E+05, 0.41803E+05, 0.42096E+05, 0.42390E+05, 0.42686E+05, 0.42984E+05, 0.43283E+05, 0.43583E+05, 0.43886E+05 &
, 0.44189E+05, 0.44495E+05, 0.44802E+05, 0.45111E+05, 0.45422E+05, 0.45734E+05, 0.46048E+05, 0.46363E+05, 0.46681E+05, 0.46999E+05 &
, 0.47320E+05, 0.47643E+05, 0.47967E+05, 0.48292E+05, 0.48620E+05, 0.48949E+05, 0.49280E+05, 0.49613E+05, 0.49947E+05, 0.50284E+05 &
, 0.50622E+05, 0.50962E+05, 0.51303E+05, 0.51646E+05, 0.51992E+05, 0.52339E+05, 0.52687E+05, 0.53038E+05, 0.53390E+05, 0.53745E+05 &
, 0.54101E+05, 0.54459E+05, 0.54819E+05, 0.55180E+05, 0.55544E+05, 0.55909E+05, 0.56276E+05, 0.56646E+05, 0.57017E+05, 0.57390E+05 &
, 0.57764E+05, 0.58141E+05, 0.58520E+05, 0.58901E+05, 0.59283E+05, 0.59668E+05, 0.60054E+05, 0.60443E+05, 0.60833E+05, 0.61226E+05 &
, 0.61620E+05, 0.62017E+05, 0.62415E+05, 0.62815E+05, 0.63218E+05, 0.63622E+05, 0.64029E+05, 0.64437E+05, 0.64848E+05, 0.65260E+05 &
, 0.65675E+05, 0.66092E+05, 0.66511E+05, 0.66931E+05, 0.67354E+05, 0.67779E+05, 0.68207E+05, 0.68636E+05, 0.69067E+05, 0.69501E+05 &
, 0.69936E+05, 0.70374E+05, 0.70814E+05, 0.71256E+05, 0.71700E+05, 0.72147E+05, 0.72596E+05, 0.73047E+05, 0.73500E+05, 0.73955E+05 &
, 0.74412E+05, 0.74872E+05, 0.75334E+05, 0.75798E+05, 0.76264E+05, 0.76733E+05, 0.77204E+05, 0.77677E+05, 0.78153E+05, 0.78630E+05 &
, 0.79110E+05, 0.79593E+05, 0.80077E+05, 0.80564E+05, 0.81054E+05, 0.81545E+05, 0.82039E+05, 0.82535E+05, 0.83034E+05, 0.83535E+05 &
, 0.84039E+05, 0.84544E+05, 0.85053E+05, 0.85563E+05, 0.86076E+05, 0.86592E+05, 0.87110E+05, 0.87630E+05, 0.88153E+05, 0.88678E+05 &
, 0.89206E+05, 0.89736E+05, 0.90268E+05, 0.90803E+05, 0.91341E+05, 0.91881E+05, 0.92423E+05, 0.92969E+05, 0.93516E+05, 0.94067E+05 &
, 0.94619E+05, 0.95174E+05, 0.95732E+05, 0.96293E+05, 0.96855E+05, 0.97421E+05, 0.97989E+05, 0.98560E+05, 0.99133E+05, 0.99710E+05 &
, 0.10029E+06, 0.10087E+06, 0.10145E+06, 0.10204E+06, 0.10263E+06, 0.10322E+06, 0.10382E+06, 0.10441E+06, 0.10501E+06, 0.10562E+06 /
 !           10           1
 DATA(QofT(          70 ,J),J=1,501)/   233.270004272461       &
, 0.26888E+03, 0.30614E+03, 0.34499E+03, 0.38535E+03, 0.42717E+03, 0.47041E+03, 0.51501E+03, 0.56095E+03, 0.60817E+03, 0.65665E+03 &
, 0.70636E+03, 0.75727E+03, 0.80934E+03, 0.86255E+03, 0.91689E+03, 0.97232E+03, 0.10288E+04, 0.10864E+04, 0.11450E+04, 0.12046E+04 &
, 0.12652E+04, 0.13268E+04, 0.13894E+04, 0.14529E+04, 0.15174E+04, 0.15828E+04, 0.16491E+04, 0.17164E+04, 0.17845E+04, 0.18535E+04 &
, 0.19233E+04, 0.19940E+04, 0.20656E+04, 0.21380E+04, 0.22113E+04, 0.22853E+04, 0.23602E+04, 0.24359E+04, 0.25123E+04, 0.25896E+04 &
, 0.26676E+04, 0.27464E+04, 0.28260E+04, 0.29064E+04, 0.29875E+04, 0.30693E+04, 0.31519E+04, 0.32352E+04, 0.33192E+04, 0.34040E+04 &
, 0.34895E+04, 0.35757E+04, 0.36627E+04, 0.37503E+04, 0.38386E+04, 0.39277E+04, 0.40174E+04, 0.41079E+04, 0.41990E+04, 0.42908E+04 &
, 0.43833E+04, 0.44765E+04, 0.45704E+04, 0.46649E+04, 0.47601E+04, 0.48560E+04, 0.49526E+04, 0.50498E+04, 0.51477E+04, 0.52463E+04 &
, 0.53456E+04, 0.54455E+04, 0.55461E+04, 0.56474E+04, 0.57493E+04, 0.58519E+04, 0.59551E+04, 0.60591E+04, 0.61637E+04, 0.62689E+04 &
, 0.63748E+04, 0.64815E+04, 0.65888E+04, 0.66966E+04, 0.68053E+04, 0.69146E+04, 0.70245E+04, 0.71351E+04, 0.72465E+04, 0.73585E+04 &
, 0.74711E+04, 0.75845E+04, 0.76985E+04, 0.78131E+04, 0.79286E+04, 0.80447E+04, 0.81615E+04, 0.82789E+04, 0.83970E+04, 0.85159E+04 &
, 0.86355E+04, 0.87557E+04, 0.88767E+04, 0.89983E+04, 0.91207E+04, 0.92438E+04, 0.93676E+04, 0.94921E+04, 0.96174E+04, 0.97433E+04 &
, 0.98700E+04, 0.99974E+04, 0.10126E+05, 0.10254E+05, 0.10384E+05, 0.10514E+05, 0.10646E+05, 0.10777E+05, 0.10910E+05, 0.11043E+05 &
, 0.11177E+05, 0.11312E+05, 0.11448E+05, 0.11584E+05, 0.11722E+05, 0.11860E+05, 0.11998E+05, 0.12138E+05, 0.12278E+05, 0.12419E+05 &
, 0.12561E+05, 0.12704E+05, 0.12848E+05, 0.12992E+05, 0.13137E+05, 0.13283E+05, 0.13430E+05, 0.13577E+05, 0.13726E+05, 0.13875E+05 &
, 0.14025E+05, 0.14176E+05, 0.14328E+05, 0.14481E+05, 0.14634E+05, 0.14789E+05, 0.14944E+05, 0.15100E+05, 0.15257E+05, 0.15415E+05 &
, 0.15574E+05, 0.15733E+05, 0.15894E+05, 0.16055E+05, 0.16218E+05, 0.16381E+05, 0.16545E+05, 0.16710E+05, 0.16876E+05, 0.17043E+05 &
, 0.17211E+05, 0.17380E+05, 0.17550E+05, 0.17720E+05, 0.17892E+05, 0.18065E+05, 0.18238E+05, 0.18413E+05, 0.18588E+05, 0.18765E+05 &
, 0.18942E+05, 0.19121E+05, 0.19300E+05, 0.19481E+05, 0.19662E+05, 0.19845E+05, 0.20028E+05, 0.20213E+05, 0.20398E+05, 0.20585E+05 &
, 0.20772E+05, 0.20961E+05, 0.21151E+05, 0.21342E+05, 0.21533E+05, 0.21726E+05, 0.21920E+05, 0.22115E+05, 0.22312E+05, 0.22509E+05 &
, 0.22707E+05, 0.22907E+05, 0.23107E+05, 0.23309E+05, 0.23511E+05, 0.23715E+05, 0.23920E+05, 0.24126E+05, 0.24334E+05, 0.24542E+05 &
, 0.24752E+05, 0.24963E+05, 0.25175E+05, 0.25388E+05, 0.25602E+05, 0.25817E+05, 0.26034E+05, 0.26252E+05, 0.26471E+05, 0.26691E+05 &
, 0.26913E+05, 0.27135E+05, 0.27359E+05, 0.27584E+05, 0.27811E+05, 0.28038E+05, 0.28267E+05, 0.28497E+05, 0.28729E+05, 0.28961E+05 &
, 0.29195E+05, 0.29430E+05, 0.29667E+05, 0.29905E+05, 0.30144E+05, 0.30384E+05, 0.30626E+05, 0.30869E+05, 0.31113E+05, 0.31359E+05 &
, 0.31606E+05, 0.31854E+05, 0.32104E+05, 0.32355E+05, 0.32608E+05, 0.32862E+05, 0.33117E+05, 0.33373E+05, 0.33631E+05, 0.33890E+05 &
, 0.34151E+05, 0.34413E+05, 0.34677E+05, 0.34942E+05, 0.35209E+05, 0.35477E+05, 0.35746E+05, 0.36017E+05, 0.36289E+05, 0.36563E+05 &
, 0.36838E+05, 0.37114E+05, 0.37393E+05, 0.37672E+05, 0.37954E+05, 0.38236E+05, 0.38520E+05, 0.38806E+05, 0.39093E+05, 0.39382E+05 &
, 0.39672E+05, 0.39964E+05, 0.40258E+05, 0.40553E+05, 0.40849E+05, 0.41147E+05, 0.41447E+05, 0.41748E+05, 0.42051E+05, 0.42355E+05 &
, 0.42661E+05, 0.42969E+05, 0.43279E+05, 0.43589E+05, 0.43902E+05, 0.44217E+05, 0.44532E+05, 0.44850E+05, 0.45169E+05, 0.45490E+05 &
, 0.45812E+05, 0.46137E+05, 0.46463E+05, 0.46791E+05, 0.47120E+05, 0.47451E+05, 0.47784E+05, 0.48119E+05, 0.48455E+05, 0.48793E+05 &
, 0.49133E+05, 0.49474E+05, 0.49818E+05, 0.50163E+05, 0.50510E+05, 0.50858E+05, 0.51209E+05, 0.51562E+05, 0.51916E+05, 0.52272E+05 &
, 0.52629E+05, 0.52989E+05, 0.53351E+05, 0.53714E+05, 0.54079E+05, 0.54446E+05, 0.54815E+05, 0.55186E+05, 0.55559E+05, 0.55934E+05 &
, 0.56310E+05, 0.56689E+05, 0.57069E+05, 0.57451E+05, 0.57836E+05, 0.58222E+05, 0.58610E+05, 0.59000E+05, 0.59393E+05, 0.59787E+05 &
, 0.60183E+05, 0.60581E+05, 0.60981E+05, 0.61383E+05, 0.61787E+05, 0.62194E+05, 0.62602E+05, 0.63012E+05, 0.63424E+05, 0.63839E+05 &
, 0.64255E+05, 0.64674E+05, 0.65094E+05, 0.65517E+05, 0.65942E+05, 0.66369E+05, 0.66798E+05, 0.67229E+05, 0.67663E+05, 0.68098E+05 &
, 0.68535E+05, 0.68976E+05, 0.69417E+05, 0.69861E+05, 0.70308E+05, 0.70756E+05, 0.71207E+05, 0.71660E+05, 0.72115E+05, 0.72572E+05 &
, 0.73032E+05, 0.73494E+05, 0.73958E+05, 0.74424E+05, 0.74893E+05, 0.75364E+05, 0.75837E+05, 0.76312E+05, 0.76790E+05, 0.77270E+05 &
, 0.77753E+05, 0.78237E+05, 0.78725E+05, 0.79215E+05, 0.79706E+05, 0.80200E+05, 0.80697E+05, 0.81196E+05, 0.81697E+05, 0.82201E+05 &
, 0.82707E+05, 0.83215E+05, 0.83727E+05, 0.84240E+05, 0.84756E+05, 0.85274E+05, 0.85795E+05, 0.86319E+05, 0.86844E+05, 0.87372E+05 &
, 0.87903E+05, 0.88437E+05, 0.88973E+05, 0.89510E+05, 0.90051E+05, 0.90595E+05, 0.91141E+05, 0.91690E+05, 0.92241E+05, 0.92794E+05 &
, 0.93350E+05, 0.93909E+05, 0.94471E+05, 0.95035E+05, 0.95602E+05, 0.96171E+05, 0.96743E+05, 0.97318E+05, 0.97895E+05, 0.98475E+05 &
, 0.99058E+05, 0.99643E+05, 0.10023E+06, 0.10082E+06, 0.10142E+06, 0.10201E+06, 0.10261E+06, 0.10321E+06, 0.10382E+06, 0.10442E+06 &
, 0.10504E+06, 0.10565E+06, 0.10626E+06, 0.10688E+06, 0.10750E+06, 0.10813E+06, 0.10875E+06, 0.10938E+06, 0.11002E+06, 0.11065E+06 &
, 0.11129E+06, 0.11193E+06, 0.11258E+06, 0.11322E+06, 0.11387E+06, 0.11453E+06, 0.11518E+06, 0.11584E+06, 0.11650E+06, 0.11717E+06 &
, 0.11784E+06, 0.11851E+06, 0.11918E+06, 0.11986E+06, 0.12054E+06, 0.12122E+06, 0.12191E+06, 0.12259E+06, 0.12329E+06, 0.12398E+06 &
, 0.12468E+06, 0.12538E+06, 0.12609E+06, 0.12679E+06, 0.12750E+06, 0.12822E+06, 0.12893E+06, 0.12965E+06, 0.13038E+06, 0.13110E+06 &
, 0.13183E+06, 0.13257E+06, 0.13330E+06, 0.13404E+06, 0.13478E+06, 0.13553E+06, 0.13628E+06, 0.13703E+06, 0.13778E+06, 0.13854E+06 &
, 0.13930E+06, 0.14007E+06, 0.14084E+06, 0.14161E+06, 0.14239E+06, 0.14316E+06, 0.14395E+06, 0.14473E+06, 0.14552E+06, 0.14631E+06 &
, 0.14711E+06, 0.14791E+06, 0.14871E+06, 0.14951E+06, 0.15032E+06, 0.15114E+06, 0.15195E+06, 0.15277E+06, 0.15359E+06, 0.15442E+06 &
, 0.15525E+06, 0.15608E+06, 0.15692E+06, 0.15776E+06, 0.15861E+06, 0.15945E+06, 0.16030E+06, 0.16116E+06, 0.16202E+06, 0.16288E+06 &
, 0.16375E+06, 0.16461E+06, 0.16549E+06, 0.16636E+06, 0.16724E+06, 0.16813E+06, 0.16902E+06, 0.16991E+06, 0.17080E+06, 0.17170E+06 /
 !           11           1
 DATA(QofT(          71 ,J),J=1,501)/   35.4500007629395       &
, 0.39949E+02, 0.44653E+02, 0.49562E+02, 0.54672E+02, 0.59980E+02, 0.65482E+02, 0.71174E+02, 0.77049E+02, 0.83103E+02, 0.89331E+02 &
, 0.95726E+02, 0.10229E+03, 0.10900E+03, 0.11587E+03, 0.12290E+03, 0.13006E+03, 0.13737E+03, 0.14482E+03, 0.15241E+03, 0.16012E+03 &
, 0.16797E+03, 0.17595E+03, 0.18405E+03, 0.19227E+03, 0.20062E+03, 0.20909E+03, 0.21767E+03, 0.22637E+03, 0.23518E+03, 0.24411E+03 &
, 0.25314E+03, 0.26229E+03, 0.27155E+03, 0.28091E+03, 0.29038E+03, 0.29995E+03, 0.30963E+03, 0.31941E+03, 0.32929E+03, 0.33927E+03 &
, 0.34935E+03, 0.35953E+03, 0.36980E+03, 0.38018E+03, 0.39064E+03, 0.40121E+03, 0.41187E+03, 0.42262E+03, 0.43346E+03, 0.44440E+03 &
, 0.45542E+03, 0.46654E+03, 0.47775E+03, 0.48905E+03, 0.50043E+03, 0.51190E+03, 0.52346E+03, 0.53511E+03, 0.54684E+03, 0.55866E+03 &
, 0.57056E+03, 0.58255E+03, 0.59462E+03, 0.60678E+03, 0.61902E+03, 0.63134E+03, 0.64374E+03, 0.65623E+03, 0.66880E+03, 0.68145E+03 &
, 0.69418E+03, 0.70699E+03, 0.71988E+03, 0.73285E+03, 0.74590E+03, 0.75902E+03, 0.77223E+03, 0.78552E+03, 0.79888E+03, 0.81233E+03 &
, 0.82585E+03, 0.83945E+03, 0.85313E+03, 0.86688E+03, 0.88072E+03, 0.89463E+03, 0.90862E+03, 0.92268E+03, 0.93682E+03, 0.95104E+03 &
, 0.96534E+03, 0.97971E+03, 0.99416E+03, 0.10087E+04, 0.10233E+04, 0.10380E+04, 0.10527E+04, 0.10676E+04, 0.10825E+04, 0.10975E+04 &
, 0.11125E+04, 0.11277E+04, 0.11429E+04, 0.11582E+04, 0.11735E+04, 0.11890E+04, 0.12045E+04, 0.12201E+04, 0.12358E+04, 0.12516E+04 &
, 0.12674E+04, 0.12833E+04, 0.12993E+04, 0.13154E+04, 0.13315E+04, 0.13478E+04, 0.13641E+04, 0.13805E+04, 0.13969E+04, 0.14135E+04 &
, 0.14301E+04, 0.14468E+04, 0.14636E+04, 0.14805E+04, 0.14974E+04, 0.15144E+04, 0.15315E+04, 0.15487E+04, 0.15660E+04, 0.15834E+04 &
, 0.16008E+04, 0.16183E+04, 0.16359E+04, 0.16536E+04, 0.16714E+04, 0.16893E+04, 0.17072E+04, 0.17252E+04, 0.17433E+04, 0.17615E+04 &
, 0.17798E+04, 0.17982E+04, 0.18166E+04, 0.18352E+04, 0.18538E+04, 0.18725E+04, 0.18913E+04, 0.19102E+04, 0.19292E+04, 0.19483E+04 &
, 0.19675E+04, 0.19867E+04, 0.20061E+04, 0.20255E+04, 0.20450E+04, 0.20646E+04, 0.20843E+04, 0.21042E+04, 0.21241E+04, 0.21440E+04 &
, 0.21641E+04, 0.21843E+04, 0.22046E+04, 0.22250E+04, 0.22454E+04, 0.22660E+04, 0.22867E+04, 0.23074E+04, 0.23283E+04, 0.23493E+04 &
, 0.23703E+04, 0.23915E+04, 0.24127E+04, 0.24341E+04, 0.24556E+04, 0.24771E+04, 0.24988E+04, 0.25206E+04, 0.25424E+04, 0.25644E+04 &
, 0.25865E+04, 0.26087E+04, 0.26310E+04, 0.26534E+04, 0.26759E+04, 0.26985E+04, 0.27212E+04, 0.27441E+04, 0.27670E+04, 0.27900E+04 &
, 0.28132E+04, 0.28365E+04, 0.28599E+04, 0.28834E+04, 0.29070E+04, 0.29307E+04, 0.29545E+04, 0.29785E+04, 0.30025E+04, 0.30267E+04 &
, 0.30510E+04, 0.30754E+04, 0.31000E+04, 0.31246E+04, 0.31494E+04, 0.31743E+04, 0.31993E+04, 0.32244E+04, 0.32497E+04, 0.32750E+04 &
, 0.33005E+04, 0.33262E+04, 0.33519E+04, 0.33778E+04, 0.34038E+04, 0.34299E+04, 0.34561E+04, 0.34825E+04, 0.35090E+04, 0.35356E+04 &
, 0.35624E+04, 0.35893E+04, 0.36163E+04, 0.36435E+04, 0.36707E+04, 0.36982E+04, 0.37257E+04, 0.37534E+04, 0.37812E+04, 0.38092E+04 &
, 0.38373E+04, 0.38655E+04, 0.38938E+04, 0.39223E+04, 0.39510E+04, 0.39798E+04, 0.40087E+04, 0.40377E+04, 0.40670E+04, 0.40963E+04 &
, 0.41258E+04, 0.41554E+04, 0.41852E+04, 0.42151E+04, 0.42452E+04, 0.42754E+04, 0.43058E+04, 0.43363E+04, 0.43669E+04, 0.43977E+04 &
, 0.44287E+04, 0.44598E+04, 0.44911E+04, 0.45225E+04, 0.45541E+04, 0.45858E+04, 0.46177E+04, 0.46497E+04, 0.46819E+04, 0.47142E+04 &
, 0.47467E+04, 0.47794E+04, 0.48122E+04, 0.48452E+04, 0.48783E+04, 0.49116E+04, 0.49451E+04, 0.49787E+04, 0.50125E+04, 0.50465E+04 &
, 0.50806E+04, 0.51149E+04, 0.51494E+04, 0.51840E+04, 0.52188E+04, 0.52538E+04, 0.52889E+04, 0.53242E+04, 0.53597E+04, 0.53954E+04 &
, 0.54312E+04, 0.54672E+04, 0.55034E+04, 0.55397E+04, 0.55763E+04, 0.56130E+04, 0.56499E+04, 0.56869E+04, 0.57242E+04, 0.57616E+04 &
, 0.57992E+04, 0.58370E+04, 0.58750E+04, 0.59132E+04, 0.59515E+04, 0.59900E+04, 0.60288E+04, 0.60677E+04, 0.61068E+04, 0.61461E+04 &
, 0.61856E+04, 0.62252E+04, 0.62651E+04, 0.63052E+04, 0.63454E+04, 0.63859E+04, 0.64265E+04, 0.64674E+04, 0.65084E+04, 0.65496E+04 &
, 0.65911E+04, 0.66327E+04, 0.66746E+04, 0.67166E+04, 0.67589E+04, 0.68013E+04, 0.68440E+04, 0.68868E+04, 0.69299E+04, 0.69732E+04 &
, 0.70167E+04, 0.70604E+04, 0.71043E+04, 0.71484E+04, 0.71928E+04, 0.72373E+04, 0.72821E+04, 0.73270E+04, 0.73722E+04, 0.74177E+04 &
, 0.74633E+04, 0.75091E+04, 0.75552E+04, 0.76015E+04, 0.76480E+04, 0.76948E+04, 0.77417E+04, 0.77889E+04, 0.78364E+04, 0.78840E+04 &
, 0.79319E+04, 0.79800E+04, 0.80283E+04, 0.80769E+04, 0.81257E+04, 0.81747E+04, 0.82240E+04, 0.82735E+04, 0.83232E+04, 0.83732E+04 &
, 0.84234E+04, 0.84739E+04, 0.85246E+04, 0.85755E+04, 0.86267E+04, 0.86782E+04, 0.87298E+04, 0.87818E+04, 0.88339E+04, 0.88863E+04 &
, 0.89390E+04, 0.89919E+04, 0.90451E+04, 0.90985E+04, 0.91522E+04, 0.92061E+04, 0.92603E+04, 0.93148E+04, 0.93695E+04, 0.94244E+04 &
, 0.94797E+04, 0.95352E+04, 0.95909E+04, 0.96469E+04, 0.97032E+04, 0.97597E+04, 0.98166E+04, 0.98736E+04, 0.99310E+04, 0.99886E+04 &
, 0.10047E+05, 0.10105E+05, 0.10163E+05, 0.10222E+05, 0.10281E+05, 0.10340E+05, 0.10400E+05, 0.10460E+05, 0.10520E+05, 0.10580E+05 &
, 0.10641E+05, 0.10702E+05, 0.10763E+05, 0.10824E+05, 0.10886E+05, 0.10948E+05, 0.11011E+05, 0.11074E+05, 0.11137E+05, 0.11200E+05 &
, 0.11263E+05, 0.11327E+05, 0.11392E+05, 0.11456E+05, 0.11521E+05, 0.11586E+05, 0.11651E+05, 0.11717E+05, 0.11783E+05, 0.11849E+05 &
, 0.11916E+05, 0.11983E+05, 0.12050E+05, 0.12118E+05, 0.12186E+05, 0.12254E+05, 0.12322E+05, 0.12391E+05, 0.12460E+05, 0.12530E+05 &
, 0.12599E+05, 0.12670E+05, 0.12740E+05, 0.12811E+05, 0.12882E+05, 0.12953E+05, 0.13025E+05, 0.13097E+05, 0.13170E+05, 0.13242E+05 &
, 0.13315E+05, 0.13389E+05, 0.13463E+05, 0.13537E+05, 0.13611E+05, 0.13686E+05, 0.13761E+05, 0.13837E+05, 0.13912E+05, 0.13989E+05 &
, 0.14065E+05, 0.14142E+05, 0.14219E+05, 0.14297E+05, 0.14375E+05, 0.14453E+05, 0.14532E+05, 0.14611E+05, 0.14690E+05, 0.14770E+05 &
, 0.14850E+05, 0.14930E+05, 0.15011E+05, 0.15092E+05, 0.15174E+05, 0.15256E+05, 0.15338E+05, 0.15421E+05, 0.15504E+05, 0.15587E+05 &
, 0.15671E+05, 0.15756E+05, 0.15840E+05, 0.15925E+05, 0.16010E+05, 0.16096E+05, 0.16182E+05, 0.16269E+05, 0.16356E+05, 0.16443E+05 &
, 0.16531E+05, 0.16619E+05, 0.16708E+05, 0.16796E+05, 0.16886E+05, 0.16976E+05, 0.17066E+05, 0.17156E+05, 0.17247E+05, 0.17338E+05 &
, 0.17430E+05, 0.17522E+05, 0.17615E+05, 0.17708E+05, 0.17801E+05, 0.17895E+05, 0.17990E+05, 0.18084E+05, 0.18179E+05, 0.18275E+05 &
, 0.18371E+05, 0.18467E+05, 0.18564E+05, 0.18661E+05, 0.18759E+05, 0.18857E+05, 0.18956E+05, 0.19055E+05, 0.19154E+05, 0.19254E+05 /
 !           11           2
 DATA(QofT(          72 ,J),J=1,501)/   23.7169990539551       &
, 0.26724E+02, 0.29868E+02, 0.33149E+02, 0.36564E+02, 0.40112E+02, 0.43790E+02, 0.47594E+02, 0.51522E+02, 0.55568E+02, 0.59731E+02 &
, 0.64006E+02, 0.68390E+02, 0.72880E+02, 0.77472E+02, 0.82165E+02, 0.86956E+02, 0.91841E+02, 0.96820E+02, 0.10189E+03, 0.10705E+03 &
, 0.11229E+03, 0.11762E+03, 0.12304E+03, 0.12854E+03, 0.13411E+03, 0.13977E+03, 0.14551E+03, 0.15132E+03, 0.15721E+03, 0.16318E+03 &
, 0.16922E+03, 0.17533E+03, 0.18152E+03, 0.18777E+03, 0.19410E+03, 0.20050E+03, 0.20697E+03, 0.21350E+03, 0.22011E+03, 0.22678E+03 &
, 0.23351E+03, 0.24032E+03, 0.24719E+03, 0.25412E+03, 0.26112E+03, 0.26818E+03, 0.27530E+03, 0.28248E+03, 0.28973E+03, 0.29704E+03 &
, 0.30441E+03, 0.31184E+03, 0.31933E+03, 0.32688E+03, 0.33449E+03, 0.34216E+03, 0.34988E+03, 0.35767E+03, 0.36551E+03, 0.37341E+03 &
, 0.38136E+03, 0.38937E+03, 0.39744E+03, 0.40557E+03, 0.41375E+03, 0.42198E+03, 0.43027E+03, 0.43862E+03, 0.44702E+03, 0.45547E+03 &
, 0.46398E+03, 0.47254E+03, 0.48115E+03, 0.48982E+03, 0.49854E+03, 0.50732E+03, 0.51615E+03, 0.52503E+03, 0.53396E+03, 0.54295E+03 &
, 0.55199E+03, 0.56108E+03, 0.57022E+03, 0.57941E+03, 0.58866E+03, 0.59796E+03, 0.60731E+03, 0.61671E+03, 0.62616E+03, 0.63567E+03 &
, 0.64522E+03, 0.65483E+03, 0.66449E+03, 0.67420E+03, 0.68396E+03, 0.69377E+03, 0.70364E+03, 0.71355E+03, 0.72352E+03, 0.73354E+03 &
, 0.74361E+03, 0.75373E+03, 0.76391E+03, 0.77413E+03, 0.78441E+03, 0.79474E+03, 0.80512E+03, 0.81555E+03, 0.82603E+03, 0.83657E+03 &
, 0.84716E+03, 0.85780E+03, 0.86849E+03, 0.87924E+03, 0.89004E+03, 0.90089E+03, 0.91179E+03, 0.92275E+03, 0.93376E+03, 0.94483E+03 &
, 0.95594E+03, 0.96711E+03, 0.97834E+03, 0.98962E+03, 0.10010E+04, 0.10123E+04, 0.10238E+04, 0.10353E+04, 0.10468E+04, 0.10584E+04 &
, 0.10701E+04, 0.10818E+04, 0.10936E+04, 0.11054E+04, 0.11173E+04, 0.11292E+04, 0.11412E+04, 0.11533E+04, 0.11654E+04, 0.11776E+04 &
, 0.11898E+04, 0.12021E+04, 0.12144E+04, 0.12268E+04, 0.12393E+04, 0.12518E+04, 0.12644E+04, 0.12770E+04, 0.12897E+04, 0.13025E+04 &
, 0.13153E+04, 0.13282E+04, 0.13411E+04, 0.13541E+04, 0.13672E+04, 0.13803E+04, 0.13935E+04, 0.14068E+04, 0.14201E+04, 0.14335E+04 &
, 0.14469E+04, 0.14604E+04, 0.14740E+04, 0.14876E+04, 0.15013E+04, 0.15150E+04, 0.15289E+04, 0.15427E+04, 0.15567E+04, 0.15707E+04 &
, 0.15848E+04, 0.15990E+04, 0.16132E+04, 0.16275E+04, 0.16418E+04, 0.16562E+04, 0.16707E+04, 0.16853E+04, 0.16999E+04, 0.17146E+04 &
, 0.17294E+04, 0.17442E+04, 0.17592E+04, 0.17741E+04, 0.17892E+04, 0.18043E+04, 0.18195E+04, 0.18348E+04, 0.18501E+04, 0.18655E+04 &
, 0.18810E+04, 0.18966E+04, 0.19122E+04, 0.19279E+04, 0.19437E+04, 0.19596E+04, 0.19755E+04, 0.19915E+04, 0.20076E+04, 0.20238E+04 &
, 0.20400E+04, 0.20563E+04, 0.20727E+04, 0.20892E+04, 0.21058E+04, 0.21224E+04, 0.21391E+04, 0.21559E+04, 0.21728E+04, 0.21897E+04 &
, 0.22068E+04, 0.22239E+04, 0.22411E+04, 0.22584E+04, 0.22757E+04, 0.22932E+04, 0.23107E+04, 0.23283E+04, 0.23460E+04, 0.23638E+04 &
, 0.23816E+04, 0.23996E+04, 0.24176E+04, 0.24358E+04, 0.24540E+04, 0.24723E+04, 0.24906E+04, 0.25091E+04, 0.25277E+04, 0.25463E+04 &
, 0.25651E+04, 0.25839E+04, 0.26028E+04, 0.26218E+04, 0.26409E+04, 0.26601E+04, 0.26794E+04, 0.26987E+04, 0.27182E+04, 0.27378E+04 &
, 0.27574E+04, 0.27771E+04, 0.27970E+04, 0.28169E+04, 0.28369E+04, 0.28571E+04, 0.28773E+04, 0.28976E+04, 0.29180E+04, 0.29385E+04 &
, 0.29591E+04, 0.29798E+04, 0.30006E+04, 0.30215E+04, 0.30425E+04, 0.30636E+04, 0.30848E+04, 0.31061E+04, 0.31275E+04, 0.31490E+04 &
, 0.31706E+04, 0.31923E+04, 0.32141E+04, 0.32360E+04, 0.32580E+04, 0.32801E+04, 0.33023E+04, 0.33246E+04, 0.33471E+04, 0.33696E+04 &
, 0.33922E+04, 0.34150E+04, 0.34378E+04, 0.34608E+04, 0.34838E+04, 0.35070E+04, 0.35303E+04, 0.35537E+04, 0.35772E+04, 0.36008E+04 &
, 0.36245E+04, 0.36483E+04, 0.36723E+04, 0.36963E+04, 0.37205E+04, 0.37448E+04, 0.37691E+04, 0.37936E+04, 0.38183E+04, 0.38430E+04 &
, 0.38678E+04, 0.38928E+04, 0.39179E+04, 0.39431E+04, 0.39684E+04, 0.39938E+04, 0.40193E+04, 0.40450E+04, 0.40708E+04, 0.40967E+04 &
, 0.41227E+04, 0.41488E+04, 0.41751E+04, 0.42014E+04, 0.42279E+04, 0.42546E+04, 0.42813E+04, 0.43082E+04, 0.43351E+04, 0.43623E+04 &
, 0.43895E+04, 0.44168E+04, 0.44443E+04, 0.44719E+04, 0.44997E+04, 0.45275E+04, 0.45555E+04, 0.45836E+04, 0.46119E+04, 0.46402E+04 &
, 0.46687E+04, 0.46974E+04, 0.47261E+04, 0.47550E+04, 0.47840E+04, 0.48132E+04, 0.48424E+04, 0.48718E+04, 0.49014E+04, 0.49311E+04 &
, 0.49609E+04, 0.49908E+04, 0.50209E+04, 0.50511E+04, 0.50814E+04, 0.51119E+04, 0.51425E+04, 0.51733E+04, 0.52042E+04, 0.52352E+04 &
, 0.52664E+04, 0.52977E+04, 0.53291E+04, 0.53607E+04, 0.53924E+04, 0.54243E+04, 0.54563E+04, 0.54885E+04, 0.55208E+04, 0.55532E+04 &
, 0.55858E+04, 0.56185E+04, 0.56514E+04, 0.56844E+04, 0.57175E+04, 0.57508E+04, 0.57843E+04, 0.58179E+04, 0.58516E+04, 0.58855E+04 &
, 0.59195E+04, 0.59537E+04, 0.59881E+04, 0.60226E+04, 0.60572E+04, 0.60920E+04, 0.61269E+04, 0.61620E+04, 0.61973E+04, 0.62327E+04 &
, 0.62682E+04, 0.63039E+04, 0.63398E+04, 0.63758E+04, 0.64120E+04, 0.64483E+04, 0.64848E+04, 0.65214E+04, 0.65582E+04, 0.65952E+04 &
, 0.66323E+04, 0.66696E+04, 0.67070E+04, 0.67446E+04, 0.67824E+04, 0.68203E+04, 0.68584E+04, 0.68966E+04, 0.69350E+04, 0.69736E+04 &
, 0.70123E+04, 0.70512E+04, 0.70903E+04, 0.71295E+04, 0.71689E+04, 0.72085E+04, 0.72482E+04, 0.72881E+04, 0.73282E+04, 0.73684E+04 &
, 0.74089E+04, 0.74494E+04, 0.74902E+04, 0.75311E+04, 0.75722E+04, 0.76135E+04, 0.76549E+04, 0.76966E+04, 0.77384E+04, 0.77803E+04 &
, 0.78225E+04, 0.78648E+04, 0.79073E+04, 0.79500E+04, 0.79928E+04, 0.80359E+04, 0.80791E+04, 0.81225E+04, 0.81660E+04, 0.82098E+04 &
, 0.82537E+04, 0.82978E+04, 0.83422E+04, 0.83866E+04, 0.84313E+04, 0.84762E+04, 0.85212E+04, 0.85664E+04, 0.86118E+04, 0.86574E+04 &
, 0.87032E+04, 0.87492E+04, 0.87954E+04, 0.88417E+04, 0.88883E+04, 0.89350E+04, 0.89819E+04, 0.90291E+04, 0.90764E+04, 0.91239E+04 &
, 0.91716E+04, 0.92195E+04, 0.92676E+04, 0.93159E+04, 0.93643E+04, 0.94130E+04, 0.94619E+04, 0.95110E+04, 0.95603E+04, 0.96097E+04 &
, 0.96594E+04, 0.97093E+04, 0.97594E+04, 0.98097E+04, 0.98601E+04, 0.99108E+04, 0.99617E+04, 0.10013E+05, 0.10064E+05, 0.10116E+05 &
, 0.10167E+05, 0.10219E+05, 0.10271E+05, 0.10324E+05, 0.10376E+05, 0.10429E+05, 0.10482E+05, 0.10535E+05, 0.10589E+05, 0.10642E+05 &
, 0.10696E+05, 0.10750E+05, 0.10804E+05, 0.10859E+05, 0.10914E+05, 0.10968E+05, 0.11024E+05, 0.11079E+05, 0.11134E+05, 0.11190E+05 &
, 0.11246E+05, 0.11302E+05, 0.11359E+05, 0.11416E+05, 0.11472E+05, 0.11530E+05, 0.11587E+05, 0.11644E+05, 0.11702E+05, 0.11760E+05 &
, 0.11818E+05, 0.11877E+05, 0.11936E+05, 0.11995E+05, 0.12054E+05, 0.12113E+05, 0.12173E+05, 0.12233E+05, 0.12293E+05, 0.12353E+05 /
 !           12           1
 DATA(QofT(          73 ,J),J=1,501)/   2896.80004882812       &
, 0.33406E+04, 0.38049E+04, 0.42890E+04, 0.47921E+04, 0.53134E+04, 0.58524E+04, 0.64085E+04, 0.69812E+04, 0.75700E+04, 0.81745E+04 &
, 0.87944E+04, 0.94291E+04, 0.10078E+05, 0.10742E+05, 0.11420E+05, 0.12111E+05, 0.12816E+05, 0.13534E+05, 0.14265E+05, 0.15008E+05 &
, 0.15765E+05, 0.16533E+05, 0.17314E+05, 0.18107E+05, 0.18911E+05, 0.19728E+05, 0.20556E+05, 0.21395E+05, 0.22246E+05, 0.23108E+05 &
, 0.23982E+05, 0.24866E+05, 0.25762E+05, 0.26669E+05, 0.27586E+05, 0.28515E+05, 0.29454E+05, 0.30405E+05, 0.31367E+05, 0.32339E+05 &
, 0.33323E+05, 0.34317E+05, 0.35323E+05, 0.36339E+05, 0.37367E+05, 0.38407E+05, 0.39457E+05, 0.40520E+05, 0.41593E+05, 0.42679E+05 &
, 0.43776E+05, 0.44885E+05, 0.46007E+05, 0.47140E+05, 0.48286E+05, 0.49445E+05, 0.50616E+05, 0.51800E+05, 0.52997E+05, 0.54207E+05 &
, 0.55431E+05, 0.56669E+05, 0.57920E+05, 0.59186E+05, 0.60466E+05, 0.61760E+05, 0.63069E+05, 0.64393E+05, 0.65733E+05, 0.67087E+05 &
, 0.68458E+05, 0.69845E+05, 0.71248E+05, 0.72667E+05, 0.74104E+05, 0.75557E+05, 0.77028E+05, 0.78517E+05, 0.80023E+05, 0.81548E+05 &
, 0.83092E+05, 0.84654E+05, 0.86236E+05, 0.87837E+05, 0.89458E+05, 0.91099E+05, 0.92760E+05, 0.94443E+05, 0.96147E+05, 0.97872E+05 &
, 0.99619E+05, 0.10139E+06, 0.10318E+06, 0.10500E+06, 0.10683E+06, 0.10870E+06, 0.11058E+06, 0.11249E+06, 0.11443E+06, 0.11639E+06 &
, 0.11837E+06, 0.12039E+06, 0.12242E+06, 0.12449E+06, 0.12658E+06, 0.12870E+06, 0.13085E+06, 0.13303E+06, 0.13523E+06, 0.13747E+06 &
, 0.13973E+06, 0.14203E+06, 0.14435E+06, 0.14671E+06, 0.14910E+06, 0.15152E+06, 0.15397E+06, 0.15646E+06, 0.15898E+06, 0.16153E+06 &
, 0.16412E+06, 0.16675E+06, 0.16941E+06, 0.17210E+06, 0.17483E+06, 0.17760E+06, 0.18041E+06, 0.18325E+06, 0.18614E+06, 0.18906E+06 &
, 0.19202E+06, 0.19503E+06, 0.19807E+06, 0.20116E+06, 0.20428E+06, 0.20746E+06, 0.21067E+06, 0.21393E+06, 0.21723E+06, 0.22058E+06 &
, 0.22397E+06, 0.22741E+06, 0.23090E+06, 0.23443E+06, 0.23801E+06, 0.24165E+06, 0.24533E+06, 0.24906E+06, 0.25284E+06, 0.25668E+06 &
, 0.26056E+06, 0.26450E+06, 0.26850E+06, 0.27255E+06, 0.27665E+06, 0.28081E+06, 0.28503E+06, 0.28930E+06, 0.29364E+06, 0.29803E+06 &
, 0.30248E+06, 0.30699E+06, 0.31157E+06, 0.31620E+06, 0.32091E+06, 0.32567E+06, 0.33050E+06, 0.33539E+06, 0.34035E+06, 0.34538E+06 &
, 0.35048E+06, 0.35565E+06, 0.36088E+06, 0.36619E+06, 0.37157E+06, 0.37702E+06, 0.38255E+06, 0.38815E+06, 0.39383E+06, 0.39958E+06 &
, 0.40541E+06, 0.41132E+06, 0.41731E+06, 0.42338E+06, 0.42954E+06, 0.43577E+06, 0.44209E+06, 0.44850E+06, 0.45499E+06, 0.46157E+06 &
, 0.46823E+06, 0.47499E+06, 0.48184E+06, 0.48878E+06, 0.49581E+06, 0.50293E+06, 0.51015E+06, 0.51747E+06, 0.52488E+06, 0.53240E+06 &
, 0.54001E+06, 0.54773E+06, 0.55555E+06, 0.56347E+06, 0.57149E+06, 0.57963E+06, 0.58787E+06, 0.59622E+06, 0.60468E+06, 0.61326E+06 &
, 0.62194E+06, 0.63074E+06, 0.63966E+06, 0.64870E+06, 0.65785E+06, 0.66713E+06, 0.67652E+06, 0.68605E+06, 0.69569E+06, 0.70546E+06 &
, 0.71536E+06, 0.72539E+06, 0.73556E+06, 0.74585E+06, 0.75628E+06, 0.76684E+06, 0.77755E+06, 0.78839E+06, 0.79937E+06, 0.81050E+06 &
, 0.82177E+06, 0.83318E+06, 0.84475E+06, 0.85646E+06, 0.86833E+06, 0.88035E+06, 0.89253E+06, 0.90486E+06, 0.91735E+06, 0.93000E+06 &
, 0.94282E+06, 0.95580E+06, 0.96895E+06, 0.98226E+06, 0.99575E+06, 0.10094E+07, 0.10232E+07, 0.10373E+07, 0.10514E+07, 0.10658E+07 &
, 0.10804E+07, 0.10951E+07, 0.11100E+07, 0.11252E+07, 0.11405E+07, 0.11560E+07, 0.11717E+07, 0.11876E+07, 0.12037E+07, 0.12200E+07 &
, 0.12365E+07, 0.12532E+07, 0.12701E+07, 0.12873E+07, 0.13046E+07, 0.13222E+07, 0.13400E+07, 0.13580E+07, 0.13763E+07, 0.13947E+07 &
, 0.14134E+07, 0.14324E+07, 0.14515E+07, 0.14710E+07, 0.14906E+07, 0.15105E+07, 0.15307E+07, 0.15511E+07, 0.15717E+07, 0.15926E+07 &
, 0.16138E+07, 0.16352E+07, 0.16569E+07, 0.16789E+07, 0.17011E+07, 0.17236E+07, 0.17464E+07, 0.17694E+07, 0.17928E+07, 0.18164E+07 &
, 0.18403E+07, 0.18645E+07, 0.18890E+07, 0.19139E+07, 0.19390E+07, 0.19644E+07, 0.19901E+07, 0.20161E+07, 0.20425E+07, 0.20692E+07 &
, 0.20962E+07, 0.21235E+07, 0.21512E+07, 0.21792E+07, 0.22075E+07, 0.22362E+07, 0.22652E+07, 0.22945E+07, 0.23243E+07, 0.23543E+07 &
, 0.23848E+07, 0.24156E+07, 0.24468E+07, 0.24783E+07, 0.25102E+07, 0.25425E+07, 0.25752E+07, 0.26083E+07, 0.26418E+07, 0.26756E+07 &
, 0.27099E+07, 0.27446E+07, 0.27797E+07, 0.28152E+07, 0.28511E+07, 0.28875E+07, 0.29243E+07, 0.29615E+07, 0.29991E+07, 0.30372E+07 &
, 0.30758E+07, 0.31148E+07, 0.31543E+07, 0.31942E+07, 0.32346E+07, 0.32754E+07, 0.33168E+07, 0.33586E+07, 0.34009E+07, 0.34437E+07 &
, 0.34870E+07, 0.35308E+07, 0.35751E+07, 0.36200E+07, 0.36653E+07, 0.37112E+07, 0.37576E+07, 0.38045E+07, 0.38520E+07, 0.39000E+07 &
, 0.39486E+07, 0.39978E+07, 0.40475E+07, 0.40978E+07, 0.41486E+07, 0.42001E+07, 0.42521E+07, 0.43047E+07, 0.43580E+07, 0.44118E+07 &
, 0.44662E+07, 0.45213E+07, 0.45770E+07, 0.46333E+07, 0.46903E+07, 0.47479E+07, 0.48062E+07, 0.48651E+07, 0.49247E+07, 0.49850E+07 &
, 0.50459E+07, 0.51076E+07, 0.51699E+07, 0.52330E+07, 0.52967E+07, 0.53612E+07, 0.54263E+07, 0.54923E+07, 0.55589E+07, 0.56263E+07 &
, 0.56945E+07, 0.57634E+07, 0.58331E+07, 0.59035E+07, 0.59748E+07, 0.60468E+07, 0.61196E+07, 0.61933E+07, 0.62677E+07, 0.63430E+07 &
, 0.64191E+07, 0.64961E+07, 0.65739E+07, 0.66526E+07, 0.67321E+07, 0.68125E+07, 0.68938E+07, 0.69760E+07, 0.70591E+07, 0.71431E+07 &
, 0.72280E+07, 0.73139E+07, 0.74007E+07, 0.74884E+07, 0.75771E+07, 0.76668E+07, 0.77574E+07, 0.78490E+07, 0.79416E+07, 0.80353E+07 &
, 0.81299E+07, 0.82256E+07, 0.83223E+07, 0.84200E+07, 0.85188E+07, 0.86187E+07, 0.87196E+07, 0.88216E+07, 0.89248E+07, 0.90290E+07 &
, 0.91344E+07, 0.92409E+07, 0.93485E+07, 0.94573E+07, 0.95672E+07, 0.96783E+07, 0.97906E+07, 0.99041E+07, 0.10019E+08, 0.10135E+08 &
, 0.10252E+08, 0.10370E+08, 0.10490E+08, 0.10611E+08, 0.10733E+08, 0.10857E+08, 0.10981E+08, 0.11108E+08, 0.11235E+08, 0.11364E+08 &
, 0.11494E+08, 0.11626E+08, 0.11758E+08, 0.11893E+08, 0.12028E+08, 0.12166E+08, 0.12304E+08, 0.12444E+08, 0.12586E+08, 0.12729E+08 &
, 0.12873E+08, 0.13019E+08, 0.13166E+08, 0.13315E+08, 0.13466E+08, 0.13618E+08, 0.13772E+08, 0.13927E+08, 0.14084E+08, 0.14242E+08 &
, 0.14402E+08, 0.14564E+08, 0.14728E+08, 0.14893E+08, 0.15060E+08, 0.15228E+08, 0.15399E+08, 0.15571E+08, 0.15744E+08, 0.15920E+08 &
, 0.16097E+08, 0.16276E+08, 0.16457E+08, 0.16640E+08, 0.16825E+08, 0.17011E+08, 0.17200E+08, 0.17390E+08, 0.17583E+08, 0.17777E+08 &
, 0.17973E+08, 0.18171E+08, 0.18372E+08, 0.18574E+08, 0.18778E+08, 0.18984E+08, 0.19193E+08, 0.19403E+08, 0.19616E+08, 0.19831E+08 &
, 0.20048E+08, 0.20267E+08, 0.20488E+08, 0.20712E+08, 0.20937E+08, 0.21165E+08, 0.21395E+08, 0.21628E+08, 0.21863E+08, 0.22100E+08 /
 !           12           2
 DATA(QofT(          74 ,J),J=1,501)/   1931.40002441406       &
, 0.22273E+04, 0.25369E+04, 0.28597E+04, 0.31951E+04, 0.35427E+04, 0.39021E+04, 0.42728E+04, 0.46547E+04, 0.50473E+04, 0.54503E+04 &
, 0.58636E+04, 0.62868E+04, 0.67197E+04, 0.71622E+04, 0.76140E+04, 0.80749E+04, 0.85448E+04, 0.90235E+04, 0.95108E+04, 0.10007E+05 &
, 0.10511E+05, 0.11023E+05, 0.11544E+05, 0.12073E+05, 0.12609E+05, 0.13153E+05, 0.13705E+05, 0.14265E+05, 0.14832E+05, 0.15407E+05 &
, 0.15990E+05, 0.16579E+05, 0.17177E+05, 0.17781E+05, 0.18393E+05, 0.19012E+05, 0.19639E+05, 0.20272E+05, 0.20914E+05, 0.21562E+05 &
, 0.22218E+05, 0.22881E+05, 0.23552E+05, 0.24230E+05, 0.24915E+05, 0.25608E+05, 0.26309E+05, 0.27017E+05, 0.27733E+05, 0.28457E+05 &
, 0.29189E+05, 0.29929E+05, 0.30677E+05, 0.31433E+05, 0.32197E+05, 0.32970E+05, 0.33751E+05, 0.34541E+05, 0.35340E+05, 0.36147E+05 &
, 0.36964E+05, 0.37790E+05, 0.38625E+05, 0.39469E+05, 0.40323E+05, 0.41187E+05, 0.42061E+05, 0.42945E+05, 0.43839E+05, 0.44744E+05 &
, 0.45659E+05, 0.46585E+05, 0.47522E+05, 0.48470E+05, 0.49429E+05, 0.50400E+05, 0.51382E+05, 0.52377E+05, 0.53383E+05, 0.54402E+05 &
, 0.55434E+05, 0.56478E+05, 0.57535E+05, 0.58605E+05, 0.59688E+05, 0.60785E+05, 0.61896E+05, 0.63021E+05, 0.64161E+05, 0.65314E+05 &
, 0.66483E+05, 0.67666E+05, 0.68865E+05, 0.70079E+05, 0.71309E+05, 0.72555E+05, 0.73817E+05, 0.75096E+05, 0.76391E+05, 0.77704E+05 &
, 0.79033E+05, 0.80380E+05, 0.81745E+05, 0.83128E+05, 0.84530E+05, 0.85950E+05, 0.87389E+05, 0.88848E+05, 0.90325E+05, 0.91823E+05 &
, 0.93341E+05, 0.94879E+05, 0.96438E+05, 0.98018E+05, 0.99620E+05, 0.10124E+06, 0.10289E+06, 0.10456E+06, 0.10625E+06, 0.10796E+06 &
, 0.10969E+06, 0.11145E+06, 0.11324E+06, 0.11505E+06, 0.11688E+06, 0.11874E+06, 0.12062E+06, 0.12253E+06, 0.12447E+06, 0.12643E+06 &
, 0.12842E+06, 0.13044E+06, 0.13248E+06, 0.13455E+06, 0.13665E+06, 0.13878E+06, 0.14094E+06, 0.14313E+06, 0.14535E+06, 0.14760E+06 &
, 0.14988E+06, 0.15219E+06, 0.15453E+06, 0.15691E+06, 0.15932E+06, 0.16176E+06, 0.16423E+06, 0.16674E+06, 0.16928E+06, 0.17186E+06 &
, 0.17448E+06, 0.17713E+06, 0.17981E+06, 0.18254E+06, 0.18530E+06, 0.18810E+06, 0.19093E+06, 0.19381E+06, 0.19673E+06, 0.19968E+06 &
, 0.20268E+06, 0.20571E+06, 0.20879E+06, 0.21192E+06, 0.21508E+06, 0.21829E+06, 0.22154E+06, 0.22483E+06, 0.22817E+06, 0.23156E+06 &
, 0.23499E+06, 0.23847E+06, 0.24200E+06, 0.24558E+06, 0.24920E+06, 0.25287E+06, 0.25660E+06, 0.26037E+06, 0.26420E+06, 0.26807E+06 &
, 0.27200E+06, 0.27599E+06, 0.28003E+06, 0.28412E+06, 0.28827E+06, 0.29247E+06, 0.29673E+06, 0.30105E+06, 0.30543E+06, 0.30986E+06 &
, 0.31436E+06, 0.31891E+06, 0.32353E+06, 0.32821E+06, 0.33296E+06, 0.33776E+06, 0.34263E+06, 0.34757E+06, 0.35257E+06, 0.35764E+06 &
, 0.36278E+06, 0.36799E+06, 0.37326E+06, 0.37861E+06, 0.38403E+06, 0.38952E+06, 0.39508E+06, 0.40072E+06, 0.40643E+06, 0.41222E+06 &
, 0.41809E+06, 0.42403E+06, 0.43005E+06, 0.43616E+06, 0.44234E+06, 0.44860E+06, 0.45495E+06, 0.46138E+06, 0.46790E+06, 0.47450E+06 &
, 0.48119E+06, 0.48797E+06, 0.49483E+06, 0.50179E+06, 0.50884E+06, 0.51598E+06, 0.52321E+06, 0.53054E+06, 0.53796E+06, 0.54548E+06 &
, 0.55310E+06, 0.56082E+06, 0.56864E+06, 0.57656E+06, 0.58458E+06, 0.59271E+06, 0.60095E+06, 0.60929E+06, 0.61773E+06, 0.62629E+06 &
, 0.63496E+06, 0.64374E+06, 0.65263E+06, 0.66164E+06, 0.67076E+06, 0.68001E+06, 0.68937E+06, 0.69885E+06, 0.70845E+06, 0.71817E+06 &
, 0.72802E+06, 0.73800E+06, 0.74810E+06, 0.75834E+06, 0.76870E+06, 0.77919E+06, 0.78982E+06, 0.80058E+06, 0.81148E+06, 0.82252E+06 &
, 0.83370E+06, 0.84502E+06, 0.85649E+06, 0.86809E+06, 0.87985E+06, 0.89175E+06, 0.90381E+06, 0.91601E+06, 0.92837E+06, 0.94088E+06 &
, 0.95355E+06, 0.96638E+06, 0.97937E+06, 0.99253E+06, 0.10058E+07, 0.10193E+07, 0.10330E+07, 0.10468E+07, 0.10608E+07, 0.10750E+07 &
, 0.10893E+07, 0.11038E+07, 0.11185E+07, 0.11334E+07, 0.11485E+07, 0.11637E+07, 0.11792E+07, 0.11948E+07, 0.12107E+07, 0.12267E+07 &
, 0.12429E+07, 0.12593E+07, 0.12759E+07, 0.12928E+07, 0.13098E+07, 0.13270E+07, 0.13445E+07, 0.13621E+07, 0.13800E+07, 0.13981E+07 &
, 0.14164E+07, 0.14350E+07, 0.14537E+07, 0.14727E+07, 0.14920E+07, 0.15114E+07, 0.15311E+07, 0.15510E+07, 0.15712E+07, 0.15916E+07 &
, 0.16123E+07, 0.16332E+07, 0.16543E+07, 0.16758E+07, 0.16974E+07, 0.17194E+07, 0.17416E+07, 0.17640E+07, 0.17867E+07, 0.18097E+07 &
, 0.18330E+07, 0.18566E+07, 0.18804E+07, 0.19045E+07, 0.19289E+07, 0.19536E+07, 0.19786E+07, 0.20039E+07, 0.20294E+07, 0.20553E+07 &
, 0.20815E+07, 0.21080E+07, 0.21348E+07, 0.21619E+07, 0.21894E+07, 0.22171E+07, 0.22452E+07, 0.22737E+07, 0.23024E+07, 0.23315E+07 &
, 0.23609E+07, 0.23907E+07, 0.24208E+07, 0.24513E+07, 0.24821E+07, 0.25133E+07, 0.25448E+07, 0.25768E+07, 0.26090E+07, 0.26417E+07 &
, 0.26747E+07, 0.27081E+07, 0.27419E+07, 0.27761E+07, 0.28107E+07, 0.28457E+07, 0.28811E+07, 0.29169E+07, 0.29531E+07, 0.29897E+07 &
, 0.30267E+07, 0.30642E+07, 0.31021E+07, 0.31404E+07, 0.31791E+07, 0.32183E+07, 0.32580E+07, 0.32981E+07, 0.33386E+07, 0.33796E+07 &
, 0.34211E+07, 0.34630E+07, 0.35054E+07, 0.35483E+07, 0.35917E+07, 0.36356E+07, 0.36800E+07, 0.37248E+07, 0.37702E+07, 0.38161E+07 &
, 0.38624E+07, 0.39094E+07, 0.39568E+07, 0.40048E+07, 0.40533E+07, 0.41023E+07, 0.41519E+07, 0.42020E+07, 0.42527E+07, 0.43040E+07 &
, 0.43558E+07, 0.44082E+07, 0.44612E+07, 0.45148E+07, 0.45689E+07, 0.46237E+07, 0.46791E+07, 0.47351E+07, 0.47916E+07, 0.48489E+07 &
, 0.49067E+07, 0.49652E+07, 0.50243E+07, 0.50841E+07, 0.51445E+07, 0.52056E+07, 0.52673E+07, 0.53298E+07, 0.53929E+07, 0.54567E+07 &
, 0.55212E+07, 0.55863E+07, 0.56522E+07, 0.57189E+07, 0.57862E+07, 0.58542E+07, 0.59230E+07, 0.59926E+07, 0.60629E+07, 0.61339E+07 &
, 0.62057E+07, 0.62783E+07, 0.63517E+07, 0.64259E+07, 0.65008E+07, 0.65766E+07, 0.66531E+07, 0.67305E+07, 0.68087E+07, 0.68878E+07 &
, 0.69677E+07, 0.70484E+07, 0.71300E+07, 0.72125E+07, 0.72958E+07, 0.73800E+07, 0.74652E+07, 0.75512E+07, 0.76381E+07, 0.77260E+07 &
, 0.78147E+07, 0.79044E+07, 0.79951E+07, 0.80867E+07, 0.81793E+07, 0.82728E+07, 0.83674E+07, 0.84629E+07, 0.85594E+07, 0.86570E+07 &
, 0.87555E+07, 0.88551E+07, 0.89557E+07, 0.90574E+07, 0.91601E+07, 0.92639E+07, 0.93688E+07, 0.94748E+07, 0.95818E+07, 0.96900E+07 &
, 0.97993E+07, 0.99097E+07, 0.10021E+08, 0.10134E+08, 0.10248E+08, 0.10363E+08, 0.10479E+08, 0.10597E+08, 0.10715E+08, 0.10835E+08 &
, 0.10956E+08, 0.11078E+08, 0.11202E+08, 0.11327E+08, 0.11453E+08, 0.11580E+08, 0.11709E+08, 0.11839E+08, 0.11970E+08, 0.12103E+08 &
, 0.12237E+08, 0.12372E+08, 0.12509E+08, 0.12647E+08, 0.12787E+08, 0.12928E+08, 0.13070E+08, 0.13214E+08, 0.13359E+08, 0.13506E+08 &
, 0.13654E+08, 0.13804E+08, 0.13955E+08, 0.14108E+08, 0.14262E+08, 0.14418E+08, 0.14575E+08, 0.14734E+08, 0.14894E+08, 0.15056E+08 /

! -------------------------------------------------------------------------------
END

!********************************************************************************
!* Program for Voigt function calculations *									*
!* <<<Kuntz,M. A new implementaition of the humlicek algorithm for the		*
!* calculation of the voigt profile function.//JQSRT,v57,n6,pp.819-824,1997.>>> *
!* ---------------------------------------------------------------------------- *
	FUNCTION VOIGT(XXX) ! X=(V-Vi)/Ad ; Y=Al/Ad
	PARAMETER (SPI=1.772454)
		common/SHAPE/ SL,AL,ADD,ALAD,VI,MOTYPE,TT,CTF2,SLL,ALAL
	SAVE A1,A2,A3,A4,A5,A6,A7,A8,B1,B2,B3,B4,B5,B6,B7,B8,	&
	C3,C4,C5,C6,C7,C8,D3,D4,D5,D6,D7,D8,E5,E6,E7,E8,	&
	F7,F8,G7,G8,H7,H8,O7,O8,P7,P8,Q7,Q8,R7,R8,S7,S8,T7,T8
	DATA Y1,Y2,Y3,Y4/4*0/
	Y=ALAD
	X=ABS(XXX/ADD)
!****		IF(X > 12.5)THEN ! Fomin
		IF(X > 15.)THEN	! Kuntz
		VOIGT=VV_LOR(XXX)
		RETURN
		END IF
	X2=X**2
	IF(X+Y >= 15.0)THEN ! Region 1
	IF(Y /= Y1)THEN
	Y1=Y
	Y_2=Y1**2
	A1=(0.2820948+0.5641896*Y_2)*Y1
	B1=0.5641896*Y1
	A2=0.25+Y_2+Y_2**2
	B2=Y_2+Y_2-1.
	END IF
		VOIGT=(A1+B1*X2)/(A2+B2*X2+X2**2)*sl/add*spi
	ELSE
		IF(X+Y >= 5.5) THEN ! Region 2
	IF(Y /= Y2)THEN
	Y2=Y
	Y_2=Y2**2
	A3=Y2*(((0.56419*Y_2+3.10304)*Y_2+4.65456)*Y_2+1.05786)
	B3=Y2*((1.69257*Y_2+0.56419)*Y_2+2.962)
	C3=Y2*(1.69257*Y_2-2.53885)
	D3=Y2*0.56419
	A4=(((Y_2+6.0)*Y_2+10.5)*Y_2+4.5)*Y_2+0.5625
	B4=((4.0*Y_2+6.0)*Y_2+9.0)*Y_2-4.5
	C4=10.5+6.0*(Y_2-1.0)*Y_2
	D4=4.0*Y_2-6.0
	END IF

	VOIGT=(((D3*X2+C3)*X2+B3)*X2+A3)/((((X2+D4)*X2+C4)*X2+B4)*X2+A4)	&
	*sl/add*spi

		ELSE
!CCCC		IF(Y >= 0.195*X-0.176)THEN ! Region 3 in accordance with Kuntz
			IF(X <= 1.0 .OR. Y >= 0.02)THEN ! Region 3 - my suggestion
	IF(Y /= Y3)THEN
	Y3=Y
	A5=((((((((0.564224*Y3+7.55895)*Y3+49.5213)*Y3+204.510)*Y3+	&
	581.746)*Y3+1174.8)*Y3+1678.33)*Y3+1629.76)*Y3+973.778)*Y3+272.102
	B5=((((((2.25689*Y3+22.6778)*Y3+100.705)*Y3+247.198)*Y3+336.364)*	&
	Y3+220.843)*Y3-2.34403)*Y3-60.5644
	C5=((((3.38534*Y3+22.6798)*Y3+52.8454)*Y3+42.5683)*Y3+18.546)*Y3+	&
	4.58029
	D5=((2.25689*Y3+7.56186)*Y3+1.66203)*Y3-0.128922
	E5=0.971457E-3+0.564224*Y3
	A6=(((((((((Y3+13.3988)*Y3+88.2674)*Y3+369.199)*Y3+1074.41)*Y3+	&
	2256.98)*Y3+3447.63)*Y3+3764.97)*Y3+2802.87)*Y3+1280.83)*Y3+	&
	272.102
	B6=(((((((5.*Y3+53.5952)*Y3+266.299)*Y3+793.427)*Y3+1549.68)*Y3+	&
	2037.31)*Y3+1758.34)*Y3+902.306)*Y3+211.678
	C6=(((((10.*Y3+80.3928)*Y3+269.292)*Y3+479.258)*Y3+497.302)*Y3+	&
	308.186)*Y3+78.866
	D6=(((10.*Y3+53.5952)*Y3+92.7586)*Y3+55.0293)*Y3+22.0353
	E6=(5.0*Y3+13.3988)*Y3+1.49645
	END IF
!*			WRITE(10,*)' REGION 3 '

	VOIGT=((((E5*X2+D5)*X2+C5)*X2+B5)*X2+A5)/	&
	(((((X2+E6)*X2+D6)*X2+C6)*X2+B6)*X2+A6)*sl/add*spi

		ELSE					! Region 4

	VOIGT=DOPLER(XXX)
		END IF
		END IF
	END IF
	END
!****************************************************************************
				FUNCTION VV_LOR(X)
!*---------------------------------------------*
!* Lorentz or Van Vleck + Huber line shape. *
!* If line cut off = 25 and line position less *
!* than 125 (in cm**-1) the VVH is recommended *
!*---------------------------------------------*
	IMPLICIT INTEGER*4 (I-N)
	REAL*8 VI
	common/SHAPE/ SL,AL,ADD,ALAD,VI,MOTYPE,TT,CTF2,SLL,ALAL
!* --------------------------------------------------------------------*
!*		For VV_LOR are used :										*
!* X	- distance from line center VI ( X = V - VI cm**-1 ),	*
!* ALAL - (Lorentz half width)**2,								*
!* SLL	- (intensity*density*half_width*VVH_factor)/pi (see LBL93), *
!* MOTYPE - to define type of the molecule : 2 -CO2, 1 -H2O, 0 -other,*
!* CTF2 - (line cut off)**2 = (25 cm**-1)**2,					*
!* --------------------------------------------------------------------*
	DATA VVH/0./
!* Line cut off *
	VV_LOR=0.
	XX=X*X
!#######			IF(XX > 100..AND.MOTYPE /= 1)RETURN
!#######			IF(XX >= 624.9.AND.MOTYPE == 1)RETURN
!*
!* Lorentz *
	VV_LOR=SLL/(XX+ALAL)
	IF(MOTYPE == 0)RETURN ! LORENTZ for H2O
!*
!* Far wing correction *
	IF(MOTYPE == 2) THEN
!* CO2 - (Kunde et al.) *
	IF(VI > 1000. .OR. XX <= 12.25)RETURN
	IF(VI < 500.)RETURN
	Y=ABS(X)
	VV_LOR=VV_LOR   !#### *EXP(-1.4*SQRT(Y-3.5))
	RETURN
									ELSE

!* H2O - (Clough et al.) *
		VV_LOR=(VV_LOR -SLL/625.)
!*------------------------------------------------------------*
!* Correction by Clough et al. for H2O continuum (foreign) *
!*------------------------------------------------------------*
		IF(XX > 9.)VV_LOR=VV_LOR*(1.-(1.-6.65*EXP(-XX/5625.))*XX/625.)
			RETURN
	END IF
	END
	FUNCTION DOPLER(X)
!*-----------------------*
!* Doppler line shape. *
!*-----------------------*
	PARAMETER (SPI=1.772454,BOUND=12.5*12.5)
	REAL*8 VI
!*	IMPLICIT INTEGER*4 (I-N)
	common/SHAPE/ SL,AL,ADD,ALAD,VI,MOTYPE,TT,CTF2,SLL,ALAL
		dimension u(9),w(9)
		DATA U/1.,1.5,2.,2.5,3.,3.5,4.,4.5,5./,	&
	W/-0.688,0.2667,0.6338,0.4405,0.2529,0.1601,0.1131,0.0853,0.068/
!* ---------------------------------------------------------------*
!*		For DOPLER are used :								*
!* X - distance from line center VI ( X = V - VI cm**-1 )	*
!* ADD - (Doppler half width)/(Ln(2))**0.5					*
!* SL - (intensity*density*VVH_factor)/pi (see LBL93.f lbf93.f) *
!* ---------------------------------------------------------------*
			XX=(X/ADD)**2
		IF(XX < BOUND)THEN
!* Lorentz correction *
			IF(XX >= 25.)THEN
				DOPLER=VV_LOR(X)*(1.+1.5/XX)
			RETURN
			END IF
!* Pure Doppler *
			DOPLER=SL*EXP(-XX)*SPI/ADD
		IF(XX <= 1.4)RETURN
!* Doppler shape begins to transform to Lorentz shape *
		XI=ABS(X/ADD)
		I=XI/0.5-1.00001
		F=2.*(W(I)*(U(I+1)-XI)+W(I+1)*(XI-U(I)))
		DOPLER=DOPLER+SLL/X**2*(1.+F)
		ELSE
!* Doppler must be changed to Lorentz *
			DOPLER=VV_LOR(X)
			RETURN
		END IF
			END

				FUNCTION VAN_VLE(T,V)
!* the radiation factor
	PARAMETER (PLANCK=6.626075E-27,BOLTZ=1.380658E-16,	&
	CLIGHT=2.99792458E10,RADCN2=PLANCK*CLIGHT/BOLTZ)
!*
			XKT=T/RADCN2
			XV=V/XKT
			IF(XV <= 0.01)THEN
			VAN_VLE=0.5*XV*V
			ELSE
			IF(XV <= 10.)THEN
			EX=EXP(-XV)
			VAN_VLE=V*(1.-EX)/(1.+EX)
				ELSE
			VAN_VLE=V
			END IF
			END IF
			END


