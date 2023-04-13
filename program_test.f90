PROGRAM test_atmosphere_module
    USE atmosphere
    IMPLICIT NONE
    
    CHARACTER(LEN=20) :: atmosphereFile

    ! Set the name of the atmosphere file you want to test
    atmosphereFile = 'CO2_gas_profile.dat'

    ! Call the read_atm_profile subroutine from the new module
    CALL read_atm_profile(atmosphereFile)
    
    ! Add any other tests or checks here if necessary

END PROGRAM test_atmosphere_module