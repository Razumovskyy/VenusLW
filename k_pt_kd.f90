SUBROUTINE PT_K_COEF(MET,N_SEP,LEV_SEP,NEW_OLD) 
USE INITIAL_PT_KD
!* --------------------------------------------------------------------- *
!*
	REAL*8 VSTART, VFINISH, VS, VFISH,VFISHKA, H 
!*
	CHARACTER COMPONENT*5,new_old*3, KISHKANAME*50
    CHARACTER DIR_NAME*4,MET*2

!*
    COMMON/KSHK/KISHKANAME
	COMMON/R_A/VSTART,SEP(0:20),KISH(10250)
	DIMENSION ICODE(NCOMP),LEV_SEP(0:20)
	save VFISH, H
	DATA ISTART/0/
!*
	IF(ISTART == 0)THEN
			ISTART=1
                  OPEN(97,FILE=KISHKANAME)
				  H=10.D0/20480.D0
				  VFISH=VSTART
    END IF  
		VS=VSTART
	WRITE(*,*)'<<< L-by-L :',VS,' cm**(-1) >>>'
!*....................Cycle on atmosphere MAC calc, levels:
                 DO JN=1,JMAX
!*	WRITE(*,*)'Level No: ',JN,' ',Z(JN),' km'
 DIR_NAME='___.'
IF(JN<10)THEN
WRITE(DIR_NAME(1:1),'(I1)')JN   
ELSE
     IF(JN<100)THEN
WRITE(DIR_NAME(1:2),'(I2)')JN
	 ELSE
WRITE(DIR_NAME(1:3),'(I3)')JN
	 END IF
END IF
OPEN(491,ACCESS='DIRECT', FORM='UNFORMATTED',	&
 RECL = NTH*4, FILE=LINE_PATH//DIR_NAME//MET) !  NTH*4                  < === NTH*4 for OTHER FORTRANs
NW=(VSTART+1.)/10.0
READ(491,REC=NW)RK
CLOSE(491)
RABMA(:,JN)=RK
END DO

!do j=JMAX, 1, -1
!    write(*,*) j, RABMA(1,J) ; pause 777
!end do
! *** Separation - when RABMA is READY. ATTENTION- 1term ONLY!  *** !
WRITE(*,*)LEV_SEP(1),SEP(1),(RABMA(JC,LEV_SEP(1)),JC=1,5)
WRITE(*,*)LEV_SEP(2),SEP(2),(RABMA(JC,LEV_SEP(2)),JC=1,5)
            IK=NTH
             IK_2=IK-2
				  IF(NEW_OLD=='NEW') THEN
              JH=0
             DO J=1,IK_2,2
              JH=JH+1
! ------------------------------------------------- !
KISH(JH)=1
IF(RABMA(J+1,LEV_SEP(2))<=SEP(2).AND.RABMA(J+1,LEV_SEP(1))>=SEP(1))KISH(JH)=2
			END DO
    			   END IF

                         DO JU=1,JH
                          VFISH=VFISH+H+H
           IF(NEW_OLD=='NEW') THEN
                     WRITE(97,*)VFISH,KISH(JU) ! writing
 197      FORMAT(I2)
                         ELSE
                     READ(97,*)VFISH, KISH(JU)  ! reading
           END IF                
                         END DO 
                        IF(NEW_OLD.eq.'NEW') then
           WRITE(97,*)-VFISH,' 0 '   ! end of interval record
                        ELSE
                read(97,*)VFISHKA     ! reading 
                        END IF  

	END


