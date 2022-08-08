! *****  SLON_AVT_SEPAR.for ****
subroutine slon_avt_separ(JMAX, nar_nar, naww, natt, N_SEP, V1, V2, BAN)
  implicit real*8 (a-h,o-z)
  parameter (NCASE=1, NL=200)
  character*20  NAW(NCASE),NAR(NCASE),NAT(NCASE), naww, natt
  character*3 BAN, nar_nar, SEPAR(20)*2
	character title*60, path*55
  real*4 XUPJ,XDOJ,QT
  dimension fu(nl), fd(nl), p(nl), ro(9), xef(nl), qt(nl), SEU(20), SED(20)
  DATA SEPAR/'1_','2_','3_','4_','5_','6_','7_','8_','9_'    &
     ,'10','11','12','13','14','15','16','17','18','19','20'/
  
  !####     path='D:\Moscow-2k2\K-Distribution\LW-SEPARATE\LBL_H\'

  NAR(1) = 'speH__separ.'//nar_nar
  NAW(1) = naww !BAN//'TRP_H2O.'//nar_nar
  NAT(1) = natt !'TRP_1_MOD.H2O'

	write(*,*) nar(1), naw(1), nat(1)

  do NSE=1, N_SEP
    do NC=1, NCASE
      open(10, file=NAR(NC))
      open (15, file=SEPAR(NSE)//NAW(NC))
!*-------------------------------------------*
      do l=1, JMAX  ! ###  Analog JMAX
        fu(l) = 0.
        fd(l) = 0.
      end do
      
      do k=1, 100000
        read(10,*)v
        do l=1, JMAX
      read(10,*)(SEU(I),I=1,N_SEP),(SED(I),I=1,N_SEP)
!          read(10,*) ZZZ, SED(I),SEU(I), I=1, N_SEP ! compilation error with ZZZ
          U = SEU(NSE)
          D = SED(NSE)
          if (v .ge. v1-0.1) then
! *** Fluxes in EACH CHANNEL *** !
            fu(l)=fu(l)+u
            fd(l)=fd(l)+d
          end if
         end do
          if (v .gt. v2-11.) go to 1
         end do

 1       continue
         close(10) 

!* ----- Pressure and Temperature ------- *
      OPEN(55,FILE=NAT(NC))

           READ(55,'(A)')TITLE
      READ(55,*)NGAS,JZ
      IF(JZ.GT.JMAX) WRITE(*,*)' Number of levels',JZ,'>' &
     ,JMAX,' in DIMENSION'
      IF(JZ.GT.JMAX) STOP
         DO I=1,NGAS
         READ(55,'(A)')ATM_PATH
         END DO
      DO J=1,JMAX
      READ(55,*)Z,P(J),T,( RO(I), I = 1, NGAS )
      P(J)=P(J)*1013.25
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
end

