! This module would be used instead of COMMON blocks in previous version of the project

MODULE shared_variables
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: numSpecPoints=20481, maxNumComponents=49, maxNumLevels=200
    CHARACTER(LEN=50) :: kishkaName ! name of the kishka file
    
	REAL(KIND=4) :: RK(numSpecPoints), RABMA(numSpecPoints, maxNumLevels) ! absorption coefficient arrays
	!! use selected_real_kind() later

    REAL(KIND=4) :: SEP(0:20) ! overrided indices starting with 1
    REAL(KIND=8) :: KISH(10250) 
!CONTAINS
    ! You can also add subroutines or functions here if needed
END MODULE shared_variables