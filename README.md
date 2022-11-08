# **VenusLW**

K-terms for radiative transfer in Venus Atmosphere in longwave region (10 - 10000 cm$^{-1}$)

## **Branch:** k1-terms-pt

Obtaining first k-terms for a given spectral band.

## **Prerequistes:**

Ubuntu 20.04 (or other linux distro), gfortran, python (for plots)

Precalculated spectral PT-tables located in /srv/PT_TABLES  directory.

### **Build project:**

1. <ins>Build first part of the project, by typing in terminal:</ins>:</br> `gfortran main_1step_pt_kd.f90 mod_pt_kd.f90 k_pt_kd.f90 dflux_pt_kd.f90 planck_separate_kd.f90 -o get1kd.out`
</br>This will create executable file get1kd.out. By executing it one can obtain first k-term in given spectral band.</br></br>

2. <ins>Build another part of the project</ins> for checking the accuracy of obtained k-term:</br> `gfortran main_compare.f90 tau_to_flux_compare.f90 planck_compare.f90 -o compare.out`

### **Steps:**

1. <ins>Fill the file **"band_settings.txt"**</ins>. Example:

    - **__2**  &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;  current channel number
    - **M** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; band identification
    - **620.0D0 670.0D0** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; band boundaries in $cm^{-1}$
    - **N2** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; label for PT-tables
    - **TRP_CO2.330** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; atmospheric model located in *./Atmospheres* directory

2. <ins>Fill the file **"channel_settings.txt"**</ins>. Set values of thresholds for lower and upper atmospheric levels. Example:

    - **0.5 20** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                  lower
    - **0.1 30** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                 upper

3. <ins>Run executable</ins>: `./get1kd.out`

    After that there will be several output files with extension that will correspond to channel number from "band_setting.txt" file.

    - **'Flux-Term_2'** is needed for calculation of k-terms in next steps. Fluxes are in W/m$^{2}$.

    | Z (km) | Fd* | Fu* | Fd | Fu |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.1617E+02 | 0.1642E+02 | 0.3401E+01 | 0.3406E+01 |
    | 1.0 | 0.1506E+02 | 0.1554E+02 | 0.3180E+01 | 0.3189E+01 |
    | 2.0 | 0.1400E+02 | 0.1457E+02 | 0.2967E+01 | 0.2975E+01 |
    | ... | .......... | .......... | .......... | .......... |

    - **'PLANCKseparate.__2'** is needed for calculation the Planck function in KD calculations. First number is the number of lines. Next table -- values of Planck function depending on temperature.

    | Temperature (K) | Value of Planck function (W/m$^{2}$) |
    | --- | --- |
    | 150.000000000000 | 0.283174060238745 | 
    | 150.100000000000| 0.284363579969140 |
    | 150.200000000000 | 0.285556513217618|
    | ................ | .................|

    - Other temporary control files **'LENTA.__2'**, **'speH_separ.__2'**, **'k_coef.in'**, **'k_coef.chk'** are located in *Control* directory.

4. <ins>Execute second binary file by typing</ins> `./compare.out`

    The main output file is **'K(Z).__2'**. Example contents of this file:

    | Z (km) | K | P (mb) | T (K) | Ro () |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.403280E-19  | 1013.25000000 | 299.70 | 0.8085E+21 |
    | 1.0 | 0.402041E-19  | 904.22430000 | 293.70 | 0.7362E+21 |
    | 2.0 | 0.438545E-19 | 805.22977500 | 287.70 | 0.6692E+21 |
    | ... | .......... | .......... | .......... | .......... |
