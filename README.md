# **VenusLW**

K-terms for radiative transfer in Venus Atmosphere in longwave region.

## **Branch**: k1-terms-pt

Obtaining first k-terms for each partition interval.

### **Steps:**

1. <ins>Fill the file **"Band_V1-V2.PT"**</ins>. Example:

    - **__2**  ~~~~~~~~~~~~~~~~~~~~~~~~~ current channel number
    - **M**    ~~~~~~~~~~~~~~~~~~~~~~~~~~~ band identification
    - **620.0D0 670.0D0** ~~~~~~~~~~~~~~ band boundaries in $cm^{-1}$
    - **N2** ~~~~~~~~~~~~~~~~~~~~~~~~~~ label for PT-tables
    - **TRP_CO2.330** ~~~~~~~~~~~~~~~~~ atmospheric model in ./Atmospheres/

2. <ins>Fill the file **"Channel_Setting.txt"**</ins>. Set values of thresholds for lower and upper atmospheric levels. Example:

    - **0.5 20** ~~~~~~~~~~~~~~~~~~~~~~~~                  lower
    - **0.1 30** ~~~~~~~~~~~~~~~~~~~~~~~~                 upper

3. <ins>Build project by running in terminal</ins>: `gfortran PT_KD_LW_main_1step.f90 MOD_INI_PT_KD.f90 PT_K_forKD.f90 PT_Dflux_forKD.f90 KD_Planck_separate_H.f90 -o run.out`

    This will create new binary file and substitute the old one automatically.

4. <ins>Execute the binary file by running</ins>: `./a.out`

    After that there will be several output files with extension that will correspond to channel number from "Band_V1-V2.PT" file.

    - **'Flux-Term_2'** is needed for calculation of k-terms in next steps. Fluxes are in $W/m^2$

    | Z (km) | Fd | Fu | Fd | Fu |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.1617E+02 | 0.1642E+02 | 0.3401E+01 | 0.3406E+01 |
    | 1.0 | 0.1506E+02 | 0.1554E+02 | 0.3180E+01 | 0.3189E+01 |
    | 2.0 | 0.1400E+02 | 0.1457E+02 | 0.2967E+01 | 0.2975E+01 |
    | ... | .......... | .......... | .......... | .......... |

    - **'PLANCKseparate.__2'** is needed for calculation the Planck function in KD calculations. First number is the number of lines. Next table -- values of Planck function depending on temperature.

    | temperature | value of Planck function |
    | --- | --- |
    | 150.000000000000 | 0.283174060238745 | 
    | 150.100000000000| 0.284363579969140 |
    | 150.200000000000 | 0.285556513217618|
    | ................ | .................|

    - Other files **'LENTA.__2'**, **'speH_separ.__2'**, **'k_coef.in'**, **'k_coef.chk'** -- are just for control.

5. <ins>Make second build</ins>: `gfortran V-D_tau_SEP_Q V-fromTAUtoFLUX.f90 V-PLANCK_oper.f90 -o run2.out`

    This will create new binary file and substitute the old one automatically.

6. <ins>Execute second binary file by typing</ins> `./run2.out`

    The main output file is **'K(Z).__2'**. Example contents of this file:

    | Z (km) | K | P (mb) | T (K) | Ro () |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.403280E-19  | 1013.25000000 | 299.70 | 0.8085E+21 |
    | 1.0 | 0.402041E-19  | 904.22430000 | 293.70 | 0.7362E+21 |
    | 2.0 | 0.438545E-19 | 805.22977500 | 287.70 | 0.6692E+21 |
    | ... | .......... | .......... | .......... | .......... |
