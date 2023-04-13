SUBROUTINE PT_K_COEF(MET, levelSeparations, newOrOld) ! N_SEP removed (not used)
    USE atmosphere
    USE shared_variables

    IMPLICIT NONE

	REAL(KIND=8) :: vStart, vFinish, VS, vFish, vFishka, H
    INTEGER :: iStart, ik, ik_2, jh, ju, jc, jn 
    INTEGER :: i, j ! loop variables
    INTEGER :: rn ! nw -> rn (record number in direct access files)

	CHARACTER(LEN=5) :: component
    CHARACTER(LEN=3), INTENT(IN) :: newOrOld 
    CHARACTER(LEN=50) :: kishkaName
    CHARACTER(LEN=4) :: dirName
    CHARACTER(LEN=2) :: MET

	INTEGER, DIMENSION(0:20), INTENT(IN) :: levelSeparations ! levelSeparations(0:20)
    INTEGER, DIMENSION(0:20) :: SEP
    !INTEGER, DIMENSION iCode(maxNumComponents) ! ICODE(NCOMP) unused

    ! if the variable intent(in) it must be an argument to subroutine !

	COMMON /KSHK/ kishkaName
	
    SAVE vFish, H ! to be explained or removed later
	DATA iStart/0/ ! to be explained or removed later
!*
	IF ( iStart == 0 ) THEN
		iStart = 1
        OPEN(97, FILE=kishkaName)
		H = 10.D0 / (numSpecPoints - 1) 
		vFish = vStart
    END IF  
		
    VS=vStart
	WRITE(*,*)'<<< L-by-L :',VS,' cm**(-1) >>>'

    DO i = 1, numLevels
        dirName='___.'
        IF (i < 10) THEN
            WRITE(dirName(1:1),'(I1)') i   
        ELSE
            IF ( i < 100) THEN
                WRITE(dirName(1:2),'(I2)') i
	        ELSE
                WRITE(dirName(1:3),'(I3)') i
	        END IF
        END IF
        
        ! reading from direct access file
        OPEN(491, ACCESS='DIRECT', FORM='UNFORMATTED',	&
        RECL = numSpecPoints*4, FILE=LINE_PATH//dirName//MET) !  NTH*4                  < === NTH*4 for OTHER FORTRANs
        
        rn = (vStart + 1.) / 10.0
        READ(491, REC=rn) RK
        CLOSE(491)
        
        RABMA(:, i) = RK
    END DO

!do j=JMAX, 1, -1
!    write(*,*) j, RABMA(1,J) ; pause 777
!end do

    ! *** Separation - when RABMA is READY. ATTENTION- 1term ONLY!  *** !

    WRITE(*,*) levelSeparations(1), SEP(1),(RABMA(JC, levelSeparations(1)), JC=1,5)
    WRITE(*,*) levelSeparations(2), SEP(2),(RABMA(JC, levelSeparations(2)), JC=1,5)
    
    IK = numSpecPoints
    IK_2 = IK - 2
	IF( newOrOld =='NEW' ) THEN
        JH = 0
    
        DO j = 1, IK_2, 2
            JH = JH + 1
! ------------------------------------------------- !
            KISH(JH)=1
            IF (RABMA(j+1, levelSeparations(2))<=SEP(2) .AND. RABMA(j+1, levelSeparations(1))>=SEP(1)) KISH(JH)=2
		END DO
    END IF

    DO JU=1,JH
        vFish = vFish+ 2 * H
        IF (newOrOld .eq. 'NEW') THEN
            WRITE(97,*) vFish, KISH(JU) ! writing
            197      FORMAT(I2)
        ELSE
            READ(97,*) vFish, KISH(JU)  ! reading
        END IF                
    END DO 

    IF(newOrOld .eq. 'NEW') THEN
        WRITE(97,*) -vFish, ' 0 '   ! end of interval record
    ELSE
        READ(97,*) vFishka     ! reading 
    END IF  
END


