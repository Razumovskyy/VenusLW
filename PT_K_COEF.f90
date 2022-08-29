SUBROUTINE PT_K_COEF(MET)  !  05 Aug.,2022. 
USE INITIAL_PT
CHARACTER DIR_NAME*4,MET*2
REAL*8 VSTART
COMMON/R_A/VSTART
save DIR_NAME       
!
!* --------------------------------------------------------------------- *
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
OPEN(491, ACCESS='DIRECT', FORM='UNFORMATTED',	&
RECL = NTH*4, FILE=LINE_PATH//DIR_NAME//MET) !                    < === NT*4 for Unix
!### RECL = NTH*4, FILE=LINE_PATH//DIR_NAME//MET) !               < === NT*4 for Unix
NW=(VSTART+1.)/10.0
READ(491,REC=NW)RK
CLOSE(491)
RABMA(:,JN)=RK
END DO
!write(*,*) (RABMA(i,1), i=1,3), (RABMA(i,10), i=1,3); PAUSE 1
END


