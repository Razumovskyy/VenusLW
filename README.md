# **VenusLW**

Getting k-distributions for radiative transfer in Venus Atmosphere in longwave region (10 - 10000 cm$^{-1}$)

## **Prerequistes:**

Ubuntu 20.04 (or other linux distro), gfortran, matplotlib (for plots)

Precalculated spectral PT-tables located in /srv/PT_TABLES  directory.

### **Building the project:**

1. Enter the folder `./k1-terms-sbs`. In this part of the project one can find k-terms one-by-one. To build this part of the project, type:</br> `./make.sh`
</br>This will create delete some temporal (previous) output files and make executables: `get1kd.out` and `compare.out`. By executing them one can obtain next k-term in given spectral band.</br></br>

2. The source code that is in folder `./k-terms-all` is for obtaining all k-terms at once by given values of thresholds and atmopsheric layers (see Step 1). To build this part of the project, again, type:</br> `./make.sh`
</br>Executables `./getnkd.out` and `./compare_tot.out` will be ready to run.

### **Steps:**

0. Start by entering the directory `k-terms-sbs`

1. <ins>Fill the file **"band_settings.txt"**</ins>. Example (for the first channel):

    - **__1**  &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;  current channel number
    - **M** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; band identification
    - **620.0D0 670.0D0** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; band boundaries in cm$^{-1}$
    - **N2** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; label for PT-tables
    - **TRP_CO2.330** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; atmospheric model located in *./Atmospheres* directory

2. <ins>Fill the file **"channel_settings.txt"**</ins>. Set values of thresholds for lower and upper atmospheric levels. Example:

    - **0.5 20** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                  lower
    - **0.1 30** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                 upper

3. <ins>Run executable</ins>: `./get1kd.out`

    After that there will be several output files with extension that will correspond to channel number from "band_setting.txt" file.

    - **'Flux-Term_1'** is needed for calculation of k-terms in next steps. Fluxes are in W/m$^{2}$.

    | Z (km) | Fd* | Fu* | Fd | Fu |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.1617E+02 | 0.1642E+02 | 0.3401E+01 | 0.3406E+01 |
    | 1.0 | 0.1506E+02 | 0.1554E+02 | 0.3180E+01 | 0.3189E+01 |
    | 2.0 | 0.1400E+02 | 0.1457E+02 | 0.2967E+01 | 0.2975E+01 |
    | ... | .......... | .......... | .......... | .......... |

    - **'PLANCKseparate.__1'** is needed for calculation the Planck function in KD calculations. First number is the number of lines. Next table -- values of Planck function depending on temperature.

    | Temperature (K) | Value of Planck function (W/m$^{2}$) |
    | --- | --- |
    | 150.000000000000 | 0.283174060238745 | 
    | 150.100000000000| 0.284363579969140 |
    | 150.200000000000 | 0.285556513217618|
    | ................ | .................|

    - Other temporary control files **'LENTA.__1'**, **'speH_separ.__1'**, **'k_coef.in'**, **'k_coef.chk'** are located in `./Control` directory.

4. <ins>Execute second binary file by typing</ins> `./compare.out`

    The main output file is **'K(Z).__1'**. Example contents of this file:

    | Z (km) | K | P (mb) | T (K) | Ro () |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.403280E-19  | 1013.25000000 | 299.70 | 0.8085E+21 |
    | 1.0 | 0.402041E-19  | 904.22430000 | 293.70 | 0.7362E+21 |
    | 2.0 | 0.438545E-19 | 805.22977500 | 287.70 | 0.6692E+21 |
    | ... | .......... | .......... | .......... | .......... |

5. Run python scripts `plot_cooling.py`, `plot_flux_down.py` or `plot_flux_up.py` to plot cooling rates and upward/downward fluxes, to compare the accuracy of the obtained k-term (current channel number in band_settings.txt file).

6. If needed repeat steps 1-5 with other input values in channels_settings.txt file.

7. Run steps 1-5 for next channel numnbers (__2, __3, etc.), by continuously keeping the compromise between number of channels and discrepancy between line-by-line and KD upward fluxes and heating rates. Be care not lose contents of channel_settings file for each confirmed channel. 

8. Suppose that first channel corresponds to lower absorption (recommended). Then the last channel should contain spectral points with greater absorption. Example of the channel_settings file for the last channel:

   - **1e+9 60** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                  lower
   - **1e-9 60** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                 upper

9. Enter the directory `./k-terms-all`. Fill the band_settings.txt file as in step 1, but current channel number should be __N, where N -- is the number of channels from steps 1-8. In channels_settings.txt file add all the different lines from steps 1-8 files. For example:

    - **1e+9 60**
    - **1e-9 60**
    - **1e-12 60**
    - **0.0 30**

10. Run executables: `./getnkd.out` and `compare_tot.out`

11. Run python scripts `plot_cooling.py` and `plot_flux.py` and check consistency of k-distributions obtained on the whole interval.