SUBROUTINE PT_K_COEF(MET, levelSeparations, newOrOld) ! N_SEP removed (not used)
    USE atmosphere

    IMPLICIT NONE

	REAL(KIND=8) :: vStart, vFinish, VS, vFish, vFishka, H
    INTEGER :: iStart, nw, ik, ik_2, jh, ju, jc, jn

	CHARACTER(LEN=5) :: component
    CHARACTER(LEN=3), INTENT(IN) :: newOrOld 
    CHARACTER(LEN=50) :: kishkaName
    CHARACTER(LEN=4) :: dirName
    CHARACTER(LEN=2) :: MET

	INTEGER, DIMENSION(0:20), INTENT(IN) :: levelSeparations ! LEV_SEP(0:20)
    INTEGER, DIMENSION(0:20) :: SEP(0:20)
    !INTEGER, DIMENSION iCode(maxNumComponents) ! ICODE(NCOMP) unused

    ! if the variable intent(in) it must be an argument to subroutine !

	COMMON /KSHK/ KishkaName
    COMMON /R_A/ vStart, SEP(0:20), KISH(10250)   ! SEP(0:20) - avoided
	
    SAVE vFish, H ! to be explained or removed later
	DATA iStart/0/ ! to be explained or removed later
!*
	IF(ISTART == 0)THEN
		ISTART = 1
        OPEN(97, FILE=kishkaName)
		H = 10.D0 / (numSpecPoints - 1) 
		vFish = vStart
    END IF  
		
    VS=VSTART
	WRITE(*,*)'<<< L-by-L :',VS,' cm**(-1) >>>'
!*....................Cycle on atmospheric levels:
    DO jn=1,numLevels
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
        
        ! reading from direct access file
        OPEN(491,ACCESS='DIRECT', FORM='UNFORMATTED',	&
        RECL = NTH*4, FILE=LINE_PATH//DIR_NAME//MET) !  NTH*4                  < === NTH*4 for OTHER FORTRANs
        
        NW=(VSTART+1.)/10.0
        READ(491, REC=NW)RK
        CLOSE(491)
        
        RABMA(:, JN) = RK
    END DO

!do j=JMAX, 1, -1
!    write(*,*) j, RABMA(1,J) ; pause 777
!end do
! *** Separation - when RABMA is READY. ATTENTION- 1term ONLY!  *** !

    WRITE(*,*) levelSeparations(1),SEP(1),(RABMA(JC,LEV_SEP(1)),JC=1,5)
    WRITE(*,*) levelSeparations(2),SEP(2),(RABMA(JC,LEV_SEP(2)),JC=1,5)
    
    IK = NTH
    IK_2 = IK-2
	IF( newOrOld =='NEW' ) THEN
        JH = 0
    
        DO J=1,IK_2,2
            JH = JH + 1
! ------------------------------------------------- !
            KISH(JH)=1
            IF (RABMA(J+1,LEV_SEP(2))<=SEP(2) .AND. RABMA(J+1,LEV_SEP(1))>=SEP(1)) KISH(JH)=2
		END DO
    END IF

    DO JU=1,JH
        vFish=vFish+H+H
        IF (newOrOld=='NEW') THEN
            WRITE(97,*)vFish,KISH(JU) ! writing
            197      FORMAT(I2)
        ELSE
            READ(97,*)vFish, KISH(JU)  ! reading
        END IF                
    END DO 

    IF(newOrOld .eq. 'NEW') THEN
        WRITE(97,*)-vFish,' 0 '   ! end of interval record
    ELSE
        READ(97,*)vFishka     ! reading 
    END IF  
END


