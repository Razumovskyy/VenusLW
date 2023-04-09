!******************************************************************
!* This program creates only 1 K-terms  using PT-tables*
!******************************************************************
PROGRAM LW_1_Term
USE atmosphere
REAL*8 VSTART,VFINISH,DIAP,FUP,FDO,S1,S2,FLUXUP,FLUXDO,V1,V2,VSOLD
PARAMETER(DIAP=10.D0,iout=47) 
CHARACTER FI*20,KISHKANAME*50
CHARACTER NAWW*20,NATT*20,MET*2
CHARACTER*3 N_K,NEW_OLD,BAN
DIMENSION LEV_SEP(0:2)     
COMMON/R_A/VSTART,SEP(0:2),KISH(10250)       
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
!****************************************!

!*  -----   GENERAL SETTING -----  *
  NATT=FI      ! Atmospheric models- READING
!*  ------------------------------ *
NEW_OLD='NEW'   ! "new" separation will be used
! -------------- Threshoulds definition ------------------------------- 
SEP(0)=0.0 ; LEV_SEP(0)=1   ! Must be defined in ANY case            
N_SEP=2   
KISHKANAME='./Control/LENTA.'//N_K      
OPEN(2011,FILE='channel_settings.txt')                              !wwwwwwwwwwwwwwwwwwWWWWW
READ(2011,*)SEP(1),LEV_SEP(1)
READ(2011,*)SEP(2),LEV_SEP(2)
CLOSE(2011)
!*********** END SETTINGS **************************************************

!* ---------------------------------------------------------------- *
V1=VSTART ; V2=VFINISH
OPEN (89,file='./Control/speH__separ.'//N_K)
OPEN(99,FILE='./Control/k_coef.in')
WRITE(99,'(A)')'./Atmospheres/'
WRITE(99,'(A)')'./HITRAN16/'
CLOSE(99)
! ----------------------------------------------------------------- !
CALL ATM_PROF_READING(NATT)
        ALLOCATE (FUP(JMAX,N_SEP),FDO(JMAX,N_SEP), &
        FLUXUP(JMAX,N_SEP),FLUXDO(JMAX,N_SEP),ZI(JMAX),STAT=IERR)
        IF(IERR/=0)THEN
                WRITE(*,*)' Allocation is wrong (main) !!!'
                STOP
        END IF
	OPEN(IOUT,FILE='Flux-Term'//N_K)
!* -------------Last line definition --------------------------------*
			DO J=1,JMAX
                        DO NN=1,N_SEP 
			FUP(J,NN)=0.D0
	        FDO(J,NN)=0.D0
                        END DO
			END DO
	WRITE(*,*)'------------------------------------------------'
	WRITE(*,*)' *** START since ',VSTART,' cm**(-1) ***'
!* ---------------------------------------------------------------- *
Nc = (VFINISH - VSTART + 1.0) / 10.D0
write(*,*) 'Nc=', Nc
DO nnc=1, Nc 
WRITE(*,*)NNC,NC
 VFINISH=VSTART+10.D0
	DO J=1,JMAX
                DO NN=1,N_SEP
	        	FLUXUP(J,NN)=0.D0
			FLUXDO(J,NN)=0.D0
		END DO
        END DO
!*
        CALL PT_K_COEF(MET,N_SEP,LEV_SEP,NEW_OLD) ! 12 April,2002.: calc spectrum in current 10 cm-1 subinterval; spectroscopy
       CALL FLUX_H(FLUXUP, FLUXDO, N_SEP)
	WRITE(89,*) VSTART, VSTART + 10.D0
	DO J=1,JMAX
                a1 = FLUXDO(J, 1)
                a2 = FLUXDO(J, 2)
                a3 = FLUXUP(J, 1)
                a4 = FLUXUP(J, 2)
                write(89, 300) Z(J), a1, a3, a2, a4
       300 FORMAT(F7.1,2(2E12.4),3X)
                AA=0.D0
                BB=0.D0
                DO NS=1,N_SEP
                        AA=AA+FLUXUP(J,NS)
                        BB=BB+FLUXDO(J,NS)
                END DO
                DO NS=1,N_SEP
			S1 = FLUXUP(J,NS)
			S2 = FLUXDO(J,NS)
			FUP(J,NS) = FUP(J,NS) +S1
			FDO(J,NS) = FDO(J,NS) +S2
	        END DO
        END DO
        VSOLD = VSTART
        VSTART = VSTART + DIAP
        VFINISH = VSTART + DIAP
END DO
!* ---------------------------------- *
!*			Writing			*
			DO J=1,JMAX
write(IOUT, 300)Z(J), (FDO(J, N1), FUP(J, N1), N1=1, N_SEP)      
			END DO
!* ---------------------------------- *
CLOSE(89) ;  CLOSE(97) ; CLOSE(198)

! ************* Finish *************** !

500 CONTINUE
 !### CALL SLON_AVT_SEPAR(JMAX,nar_nar,naww,natt,N_SEP,V1,V2,BAN)
 !###                  WRITE(*,*) ' SLON_AVT_SEPAR => OK! '

   IF(NEW_OLD.eq.'NEW')  THEN
   WRITE(*,*)' *** "NEW" -> PLANCK_SEPARATE_H(V1,V2)' 
           CALL PLANCK_SEPARATE_H(V1,V2,N_K)
    END IF
            IF(NEW_OLD.EQ.'OLD')WRITE(*,*)' Without planck_separate '
    write(*,*) 'Congratulations!'
	END
