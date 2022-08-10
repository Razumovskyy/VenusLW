!******************************************************************
!* This program creates only 1 K-terms  *
!******************************************************************
program LW_1_Term !  28-07-2022 ; 31 Dec.,2003.
use INITIAL

real*8 VSTART, VFINISH, DIAP, FUP, FDO, S1, S2, FLUXUP, FLUXDO, V1, V2, VSOLD
parameter (DIAP=10.D0, iout=47) 
character FI*20, KISHKANAME*15
character NAWW*20, NATT*20
character*3 N_K, NEW_OLD, BAN

dimension LEV_SEP(0:2)     
common/R_A/VSTART, SEP(0:2), KISH(10250)       
common/KSHK/KISHKANAME
allocatable FUP(:,:),FDO(:,:),FLUXUP(:,:),FLUXDO(:,:),ZI(:)
save FUP, FDO, FLUXUP, FLUXDO, ZI

!********** BAND SETTINGS **********!
open(2001, FILE='Band_V1-V2.txt')
read(2001,2002) N_K
read(2001,2002) BAN
2002 format(A3)
read(2001,*) VSTART, VFINISH
read(2001,2003) FI
2003 format(A20)
close(2001)
!****************************************!

!*  -----   GENERAL SETTING -----  *
!### NAWW= BAN//'TRO_CO2.'//N_K  !  WRITING
NATT = FI      ! Atmospheric models- READING

!*  ------------------------------ *

NEW_OLD = 'NEW'   ! "new" separation will be used
! -------------- Threshoulds definition ------------------------------- 
SEP(0) = 0.0 ; LEV_SEP(0) = 1   ! Must be defined in ANY case            
N_SEP = 2   
KISHKANAME = 'LENTA.'//N_K      
open(2010, FILE = 'Chanel_Setting.txt')                                    !wwwwwwwwwwwwwwwwwwWWWWW
read(2010,*) SEP(1), LEV_SEP(1)
read(2010,*) SEP(2), LEV_SEP(2)
close(2010)
!*********** END SETTINGS **************************************************

!* ---------------------------------------------------------------- *
V1 = VSTART ; V2 = VFINISH
write(*,*) VSTART, VFINISH
open (89,file=' speH__separ.'//N_K)
!OPEN (89, file=' speH_separ1.'//N_K)
!OPEN (889, file=' speH_separ2.'//N_K)
open(99, FILE='k_coef.in')
write(99, '(A)')'./Atmospheres/'
write(99,'(A)')'/srv/HITRAN16/'
CLOSE(99)
! ----------------------------------------------------------------- !
CALL ATM_PROF_READING(NATT)
        ALLOCATE (FUP(JMAX,N_SEP),FDO(JMAX,N_SEP), &
        FLUXUP(JMAX,N_SEP),FLUXDO(JMAX,N_SEP),ZI(JMAX),STAT=IERR)
        IF(IERR/=0)THEN
                WRITE(*,*)' Allocaion is wrong (main) !!!'
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
!write(*,*) 'Nc=', Nc
DO nnc=1, Nc 
        VFINISH = VSTART + 10.D0
        !DO WHILE (VSTART < VFINISH) ! Loop over DIAP=10 cm**(-1)
	do J=1, JMAX
                do NN=1,N_SEP
	        	FLUXUP(J,NN)=0.D0
			FLUXDO(J,NN)=0.D0
		end do
        end do
!*
        call K_COEF00B(N_SEP,LEV_SEP,NEW_OLD) ! 12 April,2002.: calc spectrum in current 10 cm-1 subinterval; spectroscopy
        call FLUX_H(FLUXUP, FLUXDO, N_SEP)
        
	write(89,*) VSTART, VSTART + 10.D0
	do J=1, JMAX
                a1 = FLUXDO(J, 1)
                a2 = FLUXDO(J, 2)
                a3 = FLUXUP(J, 1)
                a4 = FLUXUP(J, 2)
                write(89, 300) Z(J), a1, a3, a2, a4
                !WRITE(889, 300) Z(J), (FLUXUP(J, N1), N1=1, N_SEP)
                300 format(F7.1, 2(2E12.4), 3X)
                
                AA = 0.D0
                BB = 0.D0
                DO NS=1, N_SEP
                        AA = AA + FLUXUP(J, NS)
                        BB = BB + FLUXDO(J, NS)
                END DO
                DO NS=1, N_SEP
			S1 = FLUXUP(J,NS)
			S2 = FLUXDO(J,NS)
			FUP(J,NS) = FUP(J,NS) + S1
			FDO(J,NS) = FDO(J,NS) + S2
	        END DO
        END DO
        VSOLD = VSTART
        VSTART = VSTART + DIAP
END DO
!* ---------------------------------- *
!*			Writing			*
!			WRITE(IOUT,*)VSOLD,' - ',VFINISH,' cm**(-1)'
			DO J=1,JMAX
write(IOUT, 300)Z(J), (FDO(J, N1), FUP(J, N1), N1=1, N_SEP)      
                                !WRITE(IOUT,300)Z(J),((FDO(J,N1), FUP(J,N1)),N1=1,N_SEP)
			END DO
!* ---------------------------------- *
CLOSE(89) ;  CLOSE(97) ; CLOSE(198)

! ************* Finish *************** !

500 CONTINUE
 !#### CALL SLON_AVT_SEPAR(JMAX,nar_nar,naww,natt,N_SEP,V1,V2,BAN)
 !####                  WRITE(*,*) ' SLON_AVT_SEPAR => OK! '

   IF(NEW_OLD.eq.'NEW')  THEN
   WRITE(*,*)' *** "NEW" -> PLANCK_SEPARATE_H(V1,V2)' 
           CALL PLANCK_SEPARATE_H(V1,V2)
    END IF
            IF(NEW_OLD.EQ.'OLD')WRITE(*,*)' Without planck_separate '
    write(*,*) 'Congratulations!'
	END
