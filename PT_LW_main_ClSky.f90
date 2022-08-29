!******************************************************************
!* This program calculates FLUXES & Cooling Rates              *
!*  PT-tables are used                                                                     *
!******************************************************************
PROGRAM PT_LW_main_ClSky  !  05-08-2022 ; 31 Dec.,2003.
USE INITIAL_PT
  REAL*8 VSTART,VFINISH,DIAP,S1,S2,V1,V2,VSOLD 
   PARAMETER(DIAP=10.D0,iout=47) 
    CHARACTER FIA*20,FIW*20,MET*2
!*
COMMON/R_A/VSTART       
LINE_PATH='/home/taskerbliss/fortran_projects/K-1_terms_Creation_WIND/PT_tables/'
!********** SETTINGS **********!
OPEN(2001,FILE='PT-CONTR_ClSky.LW')
READ(2001,*)VSTART,VFINISH
READ(2001,2003)FIA
READ(2001,2002)MET
2002 FORMAT(A2)
READ(2001,2003)FIW
2003 FORMAT(A20)
CLOSE(2001)
!*********** END SETTINGS **************************************************

!* ---------------------------------------------------------------- *
       V1=VSTART ; V2=VFINISH
! --------------------------- !
  CALL ATM_PROF_READING(FIA)
FUP=0.D0 ; FDO=0.D0 ;  XQ=0.D0
!* ---------------------------------------------------------------- *
Nc = (VFINISH - VSTART + 1.0) / 10.D0
write(*,*) 'Nc=', Nc
DO nnc=1, Nc 
 VFINISH=VSTART+10.D0
 WRITE(*,*)NNC,NC
FLUXUP=0.D0   ;FLUXDO=0.D0
!*
CALL PT_K_COEF(MET) ! 12 April,2002.
CALL FLUX_PT 
 FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO
VSOLD=VSTART
VSTART=VSTART+DIAP
		END DO
! ************* Finish *************** !
 WRITE(*,*)' *** OK! ***'
! ************************** !
!*			Writing			*
! ==================================================!

!###  CP_CONST=8.442 ; WRITE(*,*)' *** ',CP_CONST,' ->  Earth ***'
!###   CP_CONST=3.963 ; WRITE(*,*)' *** ',CP_CONST,' ->  Mars ***'
CP_CONST=3.963 ; WRITE(*,*)' *** ',CP_CONST,' ->  Venus ***'

 OPEN(567,FILE=FIW)
		DO J=JMAX,1,-1 !  100,99,98...km
        QT=0.
        IF(J>1) THEN
        S1=FDO(J-1)-FUP(J-1)
        S2=FDO(J)-FUP(J)
        QT=-(S1-S2)/(P1(J)-P1(J-1))*CP_CONST/1013.25 ! ***
		XQ(J)=QT
        END IF
!	WRITE(*,*)Z(J),QT,FUP(J),FDO(J) ! Here ZI = ZJ is used
	WRITE(567,172)Z(J),P1(J),T1(J),FUP(J),FDO(J),QT ! Here ZI = ZJ is used
 172		FORMAT(F7.2,E14.6,F8.2,2E14.7,F7.3)
			END DO
			XQ(1)=XQ(2)
  CLOSE(567)

!-------------- Plotting a graph  --------------- !
 ! Upward and Downward fluxes
                   open(98765,file='plot')  
  write(98765,*)' 2 ',jmax 
  do j=1,jmax 
   write(98765,*)z(j),fup(j),fdo(j)
  end do
                   close(98765)
 
 ! Cooling rates
                   open(98765,file='plot')  
  write(98765,*)' 1 ',jmax 
  do j=1,jmax 
   write(98765,*)z(j),xq(j)
  end do
                   close(98765)

	END
