# **VenusLW**

Getting k-distributions for radiative transfer in Venus Atmosphere in longwave region (10 - 10000 cm$^{-1}$)

## **Prerequistes:**

Ubuntu 20.04 (or other linux distro), gfortran, matplotlib (for plots)

Precalculated spectral PT-tables located in /srv/PT_TABLES/  directory.

### **Building the project:**

1. Enter the folder `./k1-terms-sbs`. In this part of the project one can find k-terms one-by-one. To build this part of the project, type:</br> `./make.sh`
</br>This will delete previous executables, some temporal files and create executables: `get1kd.out` and `compare.out`. By executing them one can obtain next k-term in given spectral band.</br></br>

1. The source code that is in folder `./k-terms-all` is for obtaining all k-terms at once by given values of thresholds and atmopsheric layers (see Step 1). To build this part of the project, again, type:</br> `./make.sh`
</br>Executables `./getnkd.out` and `./compare_tot.out` will be ready to run.

### **Steps:**

0. Start by entering the directory `k-terms-sbs`

1. <ins>Fill the file **"band_settings.txt"**</ins>. Example (for the first channel):

    - **__1**  &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;  current channel number
    - **M** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; band identification (set label you will be comfortable with)
    - **10.0D0 200.0D0** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; band boundaries in cm$^{-1}$
    - **H2** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; extension of particular PT-tables files, which are located in /srv/PT_TABLES/ directory. We set different extensions for different species and concentration profile.
    - **H2O_gas_profile.dat** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; atmospheric model located in *./Atmospheres* directory

2. <ins>Fill the file **"channel_settings.txt"**</ins>. Set values of thresholds for lower and upper atmospheric levels. Example:

    - **0.0 1** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                  first number -- lower threshold (smallest absorption coefficient), second number -- height level number
    - **1e15 50** &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot; &sdot;                  first number -- upper threshold (largest absorption coefficient), second number -- height level number
  
  In this file one can set lower and upper boundaries for spectral points to be selected during one iteration of obtaining k-terms. In the example presented we deliberately set 0.0 and 10$^{15}$ lower and upper thresholds. It means that all spectral points would be selected and one gets only one k-term over the whole interval 10 - 200 cm$^{-1}$. This results in a serious losing of accuracy, so one needs to decrease spectral boundaries in *channel_settings.txt* file and get more k-terms with bigger accuracy.

3. <ins>Run executable</ins>: `./get1kd.out`

    After that there will be several output files with extension that will correspond to channel number from "band_setting.txt" file.

    - **'Flux-Term_1'** is needed for calculation of k-terms in next steps. These file presents two sets of upward and downward fluxes (in W/m$^{2}$). First set (first two columns) is the net flux for not selected spectral points, second set -- for selected points. As it seen from sample figures, because all points were selected first two columns turn out to have zero values.

    | Z (km) | Fd* | Fu* | Fd | Fu |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.00 | 0.00 | 0.3401E+01 | 0.3406E+01 |
    | 1.0 | 0.00 | 0.00 | 0.3180E+01 | 0.3189E+01 |
    | 2.0 | 0.00 | 0.00 | 0.2967E+01 | 0.2975E+01 |
    | ... | .......... | .......... | .......... | .......... |

    - **'PLANCKseparate.__1'** is needed for calculation the Planck function in KD calculations. First number is the number of lines. Next table -- values of Planck function depending on temperature. This file present total Planck function value for selected spectral points and for given temperature.

    | Temperature (K) | Value of Planck function (W/m$^{2}$) |
    | --- | --- |
    | 150.00 | 0.283 | 
    | 150.10| 0.284 |
    | 150.20 | 0.286 |
    | ........... | .........|

    - Other temporary control files **'LENTA.__1'**, **'speH_separ.__1'**, **'k_coef.in'**, **'k_coef.chk'** are located in `./Control` directory. One may refer to these files in case of unexpected results or for simplifying debugging process.

4. <ins>Execute second binary file by typing</ins> `./compare.out`

    The main output file is **'K(Z).__1'**. This file contains the following colomuns: height (km), effective cross-section (cm$^2$), temperature (K), concentraion(molecules/(cm$^2$ * km)). Desired k-term is presented as effective cross-section profile.

    | Z (km) | K | P (mb) | T (K) | Ro () |
    | ------ | -- | -- | -- | -- |  
    | 0.0 | 0.403E-19  | 1013.250 | 299.70 | 0.808E+21 |
    | 1.0 | 0.402E-19  | 904.224 | 293.70 | 0.736E+21 |
    | 2.0 | 0.439E-19 | 805.229 | 287.70 | 0.669E+21 |
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
    - **0.0 1**

10. Run executables: `./getnkd.out` and `compare_tot.out`

11. Run python scripts `plot_cooling.py` and `plot_flux.py` and check consistency of k-distributions obtained on the whole interval.