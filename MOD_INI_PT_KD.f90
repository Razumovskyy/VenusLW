 MODULE INITIAL_PT_KD ! Information from atmospheric profile
       PARAMETER (NTH=20481, NCOMP = 49, JIM=200)  ! #### DLT8
	CHARACTER LINE_PATH*15  
	CHARACTER*80 TITLE
	CHARACTER*5  MOLECULE(NCOMP)
   INTEGER*4 :: NGAS,JMAX 
   REAL*4  Z(JIM),P1(JIM),T1(JIM),RO1(10,JIM),RK(NTH),RABMA(NTH,JIM)
!### REAL*8 FUP(JIM),FDO(JIM),FLUXUP(JIM),FLUXDO(JIM),ZI(JIM),XQ(JIM)
!                                                                                                                   *** 5.08.2022 ***  
CONTAINS
  SUBROUTINE ATM_PROF_READING(NATT)
  CHARACTER NATT*20
        OPEN(66,FILE='k_coef.chk')
	OPEN(55,FILE='./Atmospheres/'//NATT)
!*
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
	DO J=1,JMAX
	READ(55,*)Z(J),P1(J),T1(J),( RO1(I,J), I = 1, NGAS )
	WRITE(66,*)Z(J),P1(J),T1(J),( RO1(I,J), I = 1, NGAS )
	END DO
		CLOSE(55)
             CLOSE(66)
  END SUBROUTINE ATM_PROF_READING
  END MODULE INITIAL_PT_KD 
