! *****  SLON_AVT_SEPAR.for ****
         SUBROUTINE SLON_AVT_SEPAR(JMAX,NATT,N_SEP,V1,V2,BAN)
         IMPLICIT REAL*8 (A-H,O-Z)
          PARAMETER (NCASE=1,NL=200) !               <=    e.g. NCASE=6  !!!!
CHARACTER*50  NAW(NCASE),NAR(NCASE),NAT(NCASE)
CHARACTER BAN*3,SEPAR(20)*2,NATT*20
        CHARACTER TITLE*60, PATH*14
         REAL*4 XUPJ,XDOJ,QT
DIMENSION FU(NL),FD(NL),P(NL),RO(9),XEF(NL),QT(NL),SEU(20),SED(20)
      DATA SEPAR/'1_','2_','3_','4_','5_','6_','7_','8_','9_' ,'10','11','12','13','14','15','16','17','18','19','20'/
  PATH='./Atmospheres/'

! *** Attention: here only ONE Atmospheric Model is supposed ***
         NAR(1)='Ch-Fd_Fu(LBL-10cm).'//BAN  
         NAW(1)=BAN//NATT      
          NAT(1)=PATH//NATT 

WRITE(*,62)NAR(1) ; WRITE(*,62)NAW(1)  ; WRITE(*,62)NAT(1) 
62 FORMAT(A50) ! ; ### 30 ### pause

 DO NSE=1,N_SEP
         DO NC=1,NCASE
         OPEN(10,FILE=NAR(NC))
       OPEN (15,FILE=SEPAR(NSE)//NAW(NC))
!*-------------------------------------------*
         FU=0. ; FD=0.
         DO K=1,100000
         READ(10,*)V
!         write(*,*)v
         DO L=1,JMAX
  READ(10,*)ZZZ,(SED(I),SEU(I),I=1,N_SEP)
          U=SEU(NSE)
          D=SED(NSE)
          IF(V.GE.V1-0.1)THEN
! *** Fluxes in EACH CHANNEL *** !
         FU(L)=FU(L)+U ;  FD(L)=FD(L)+D
           END IF
         END DO
          IF(V.GT.V2-11.)GO TO 1
         END DO
1    CONTINUE
         CLOSE(10) 

!* ----- Pressure and Temperature ------- *
      OPEN(55,FILE=NAT(NC))
           READ(55,'(A)')TITLE
      READ(55,*)NGAS,JZ
      IF(JZ.GT.JMAX) WRITE(*,*)' Number of levels',JZ,'>'  ,JMAX,' in DIMENSION'
      IF(JZ.GT.JMAX) STOP
         DO I=1,NGAS
         READ(55,'(A)')ATM_PATH
         END DO
      DO J=1,JMAX
      READ(55,*)Z,P(J),T,( RO(I), I = 1, NGAS )
      P(J)=P(J)*1013.25  !   <====   mBAR
      END DO
      CLOSE(55) 

!*  -----------  QOOLING RATES, PRINT, etc. -------------- *
            DO J=1,JMAX
            XEF(J)=FU(J)-FD(J)
            END DO
! **********************************************************************
!*   'LAYER'      Derivatives                 *
         DO J=2,JMAX
      QT(J)=(XEF(J)-XEF(J-1))/(P(J)-P(J-1))*8.442
         END DO
            DO J=JMAX,2,-1
	XUPJ=FU(J)
	XDOJ=FD(J)
      WRITE(15,*)J-1,XDOJ,XUPJ,-QT(J)
            END DO
            J=1
	XUPJ=FU(J)
	XDOJ=FD(J)
      WRITE(15,*)J-1,XDOJ,XUPJ,-QT(2) ! <= In Each Channel !!!
!* ---------------------------------- *
          CLOSE(15)
         END DO
           END DO
         END

