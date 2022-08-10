# VenusLW
K-terms for radiative transfer in Venus Atmosphere in longwave region

Main programs package for obtaining K-terms "one-by-one"

**Actions**:

1. Fill the file: 'Band_V1-V2.txt'
Example:
--------------------------------
__2                     Current channel number   
_M_                     Band  identification
620.0D0  670.D          Band boundaries (Fortran double precision)
TRP_CO2.330             Atmospheric model which is located in  ./Atmospheres/
--------------------------------

2. Fill the file:  'Chanel_Setting.txt'
Example (thresholds and atmospheric levels- LOWER and UPPER, respectively) :
--------------------------------
  0.5 20                LOWER
  0.1 30                UPPER
--------------------------------
 ============================================================

3.  Run program by typing command ./run.out
After run you get files (for this example with extension '__2'):

*  'Flux-Term__2' -- contains fluxes for KD calculations in  W/m**2:
                  (Fd* and Fu* are not needed)   
 Z (km)      Fd*                Fu*                     Fd               Fu  
------------------------- example ----------------------------------
    0.0  0.1617E+02  0.1642E+02  0.3401E+01  0.3406E+01
    1.0  0.1506E+02  0.1554E+02  0.3180E+01  0.3189E+01
    2.0  0.1400E+02  0.1457E+02  0.2967E+01  0.2975E+01
  ...........
-------------------------------------------------------------------------
* 'PLANCKseparate.__2' - is needed for getting PLanck function in KD calculations.
------------------------- example ----------------------------------
        2001                                     - number of lines
   150.000000000000       0.283174060238745      - temperature in K and value of Planck function
   150.100000000000       0.284363579969140                       ...
   150.200000000000       0.285556513217618     
  ...........
-------------------------------------------------------------------------

Other files 'LENTA.__2',  'speH_separ.__2', 'k_coef.in',  'k_coef.chk' are auxillary for control.

4. 'V-D_tau_SEP_Q.exe' # make it Unix like!!
After running executable you can get needed file with K-terms: 'K(Z).__2'
  Z(km)    K                          P(mb)               T (K)       Ro    
------------------------- example ----------------------------------
   0.0  0.403280E-19   1013.25000000  299.70  0.8085E+21
    1.0  0.402041E-19    904.22430000  293.70  0.7362E+21
    2.0  0.438545E-19    805.22977500  287.70  0.6692E+21
    3.0  0.381895E-19    715.15185000  283.70  0.6029E+21
  ...........
-------------------------------------------------------------------------
Also there are 3 scripts for LBL and KD fluxes and cooling for control
     
