MODULE atmosphere
	USE shared_variables
	IMPLICIT NONE
	
	CHARACTER(LEN=15) :: LINE_PATH
	CHARACTER(LEN=80) :: atmosphereTitle ! title or description of the atmosphere profile file
	CHARACTER(LEN=5) :: molecName(maxNumComponents) ! array of molecule names presented in atmosphere profile
    INTEGER :: numComponents, numLevels !! use selected_int_kind() later  
    
	REAL(KIND=4) :: height(maxNumLevels),pressure(maxNumLevels),temperature(maxNumLevels) ! gas profile data from file
    REAL(KIND=4) :: concentration(10, maxNumLevels) ! components concentrations

	CONTAINS
	SUBROUTINE read_atm_profile(atmosphereFile) ! Information from atmospheric profile
    	CHARACTER(LEN=20), INTENT(IN) :: atmosphereFile
		INTEGER :: i, j ! loop variables

        OPEN(12, FILE='./Control/atmosphere_check.txt')
		OPEN(11, FILE='./Atmospheres/'//atmosphereFile)

		READ(11,'(A)') atmosphereTitle
		WRITE(*,*) atmosphereTitle
		WRITE(12,*) atmosphereTitle
		READ(11,*) numComponents, numLevels
		
		WRITE(*,*)'The number of gases (numComponents) : ', numComponents
		WRITE(*,*)'The number of levels (numLevels) : ', numLevels
		WRITE(12,*)'The number of gases (numComponents) : ', numComponents
		WRITE(12,*)'The number of levels (numLevels) : ', numLevels
		
		DO i=1, numComponents
			READ(11,'(A)') molecName(i)
		END DO
		
		WRITE(12,*)'Atmospheric components: ',	&
				(molecName(i), i=1, numComponents)
				
		DO j=1,numLevels
			READ(11,*) height(j),pressure(j),temperature(j), (concentration(i,j), i = 1, numComponents)
			WRITE(12,*) height(j),pressure(j),temperature(j),(concentration(i,j), i = 1, numComponents)
		END DO
		
		CLOSE(11)
        CLOSE(12)
  	END SUBROUTINE read_atm_profile
  END MODULE atmosphere
