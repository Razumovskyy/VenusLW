 PROGRAM LW_TOT_Terms  ! 06-07-2022 
 USE INITIAL_PT_KD 
REAL*8 VSTART,VFINISH,DIAP,FUP,FDO,S1,S2,FLUXUP,FLUXDO,V1,V2,VSOLD
 PARAMETER(DIAP=10.D0,IOUT=47) 
		CHARACTER FI*20,KISHKANAME*15                                                     
		CHARACTER NAWW*20,NATT*20,MET*2
       CHARACTER*3 N_K,NEW_OLD,BAN
!*
    DIMENSION LEV_SEP(0:20)
    COMMON/R_A/VSTART,SEP(0:20),KISH(10250)
    COMMON/KSHK/KISHKANAME
	ALLOCATABLE FUP(:,:),FDO(:,:),FLUXUP(:,:),FLUXDO(:,:),ZI(:)
save FUP, FDO, FLUXUP, FLUXDO, ZI
LINE_PATH='/srv/PT_TABLES/'

!********** BAND SETTINGS **********!
OPEN(2001,FILE='band_settings.txt')
READ(2001,2002)N_K
READ(2001,2002)BAN
2002 FORMAT(A3)
READ(2001,*)VSTART,VFINISH
READ(2001,2022)MET
2022 FORMAT(A2)
READ(2001,2003)FI
2003 FORMAT(A20)
CLOSE(2001)
NATT=FI      ! Atmospheric models- READING
!*  ------------------------------ *

! -------------- Threshoulds definition ------------------------------- 
OPEN(2001,FILE='channels_settings.txt')
READ(2001,*)N_SEP   ! <***  NUMBER OF  K-terms!!!
      IF(N_SEP>15)THEN
      WRITE(*,*)' *** Number of Channels should be < 16 *** ', N_SEP, ' !!!'
      PAUSE
      END IF 
DO J=N_SEP,0,-1
READ(2001,*)SEP(J),LEV_SEP(J)
END DO
 KISHKANAME='LENTA.'//N_K      
 CLOSE(2001)
!*********** END SETTINGS **************************************************
!****************************************!

         NEW_OLD='NEW'   ! "new" separation will be used
!* ---------------------------------------------------------------- *
       V1=VSTART ; V2=VFINISH
OPEN (89,FILE='Ch-Fd_Fu(LBL-10cm).'//BAN)

	OPEN(99,FILE='k_coef.in')
WRITE(99,'(A)')'./HITRAN16/'
    CLOSE(99)
! --------------------------- !
    CALL ATM_PROF_READING(NATT)
    ALLOCATE (FUP(JMAX,N_SEP),FDO(JMAX,N_SEP), &
     FLUXUP(JMAX,N_SEP),FLUXDO(JMAX,N_SEP),ZI(JMAX),STAT=IERR)
         IF(IERR/=0)THEN
          WRITE(*,*)' Allocaion is wrong (main) !!!'
          STOP
         END IF
		OPEN(IOUT,FILE='Fl_N-'//FI)
		WRITE(IOUT,*)N_SEP
!*	Last line definition *
			DO J=1,JMAX
                        DO NN=1,N_SEP 
			FUP(J,NN)=0.D0
	        FDO(J,NN)=0.D0
                        END DO
			END DO
	WRITE(*,*)'------------------------------------------------'
	WRITE(*,*)' *** START since ',VSTART,' cm**(-1) ***'
!* ---------------------------------------------------------------- *
		DO WHILE (VSTART<VFINISH) ! Loop over DIAP=10 cm**(-1)
VSOLD=VSTART
			DO J=1,JMAX
                        DO NN=1,N_SEP
			FLUXUP(J,NN)=0.D0
			FLUXDO(J,NN)=0.D0
			END DO
                        END DO
!*
        CALL PT_K_COEF(MET,N_SEP,LEV_SEP,NEW_OLD) !
        CALL FLUX_H(FLUXUP,FLUXDO,N_SEP) 
!*
			WRITE(89,*)VSTART,VSTART+DIAP
			DO J=1,JMAX
                                                !a1 = FLUXDO(J, 1)
                                                !a2 = FLUXDO(J, 2)
                                                !a3 = FLUXUP(J, 1)
                                                !a4 = FLUXUP(J, 2)
                                                !write(89, 300) Z(J), a1, a3, a2, a4
                                                write(89, 300) Z(J), (FLUXDO(J, N1), FLUXUP(J, N1), N1=1,N_SEP)
300 FORMAT(F7.1,15(2E12.4),2X)
            AA=0.D0
            BB=0.D0
               DO NS=1,N_SEP
               AA=AA+FLUXUP(J,NS)
               BB=BB+FLUXDO(J,NS)
               END DO
                       DO NS=1,N_SEP
			S1=FLUXUP(J,NS)
			S2=FLUXDO(J,NS)
			FUP(J,NS)=FUP(J,NS) +S1
			FDO(J,NS)=FDO(J,NS) +S2
	         	END DO
                       END DO
			VSTART=VSTART+DIAP
		END DO
!* ---------------------------------- *
!*			Writing			*
			DO J=1,JMAX
                                write(IOUT, 300)Z(J), (FDO(J, N1), FUP(J, N1), N1=1, N_SEP)
			END DO
!* ---------------------------------- *
        CLOSE(89)
		CLOSE(97) ; CLOSE(198)

! ************* Finish *************** !
500 CALL SLON_AVT_SEPAR(JMAX,NATT,N_SEP,V1,V2,BAN)
  write(*,*) ' slon_avt_separ => OK! '

   if(NEW_OLD.eq.'NEW')  THEN
   WRITE(*,*)' *** "NEW" -> planck_separate_H(V1,V2)' 
          CALL PLANCK_SEPARATE(V1,V2)
!###              CALL PLANCK_SEPARATE_H(V1,V2,N_K)
 !###   END IF
 
	
	END IF
            IF(NEW_OLD=='OLD')WRITE(*,*)' Without planck_separate '
!*
	END
