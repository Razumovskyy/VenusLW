[![Project Status: Active – The project is under active development](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![Version 1.0](https://img.shields.io/badge/version-1.0-blue.svg)](https://github.com/Razumovskyy/VenusLW)

# **VenusLW**

Getting k-distributions for radiative transfer in Venus Atmosphere in longwave region (10 - 6000 cm<sup>-1</sup>). After running the executables one can obtain the next k-term in a given spectral band. 

Note

- obtaining k-terms is an iterative process and requires running this programs several times.
- resulting parameterizations are strongly affected by input atmospheric profile and chosen spectroscopy (look in PTTables repository).

## **Prerequisites:**

- Ubuntu 20.04 (or other Linux distro)
- gfortran compiler
- matplotlib (for plots)
- precalculated spectral PT-tables located in `/srv/PT_TABLES/` directory. PT-tables should be available for CO<sub>2</sub>, SO<sub>2</sub> and H<sub>2</sub>O atmospheric components.

## **Building the project:**

   To build, type: `./make.sh` in the terminal. Make sure you are in the root directory of the project.

## **Steps:**

1. **Fill the input file "band_settings.txt".**
   Example (for the first channel):

    __1 M
    10.0D0 200.0D0
    H2
    H2O_gas_profile.dat
   
2. **Fill the file "channel_settings.txt".** 
Set values of thresholds of absorption coeffcient for lower and upper atmospheric levels. Example:

    0.0 1
    1e15 50

3. **Run the executables**

4. **Run Python scripts** `plot_cooling.py`, `plot_flux_down.py`, or `plot_flux_up.py` to plot cooling rates and upward/downward fluxes, to compare the accuracy of the obtained k-term (current channel number in `band_settings.txt` file).

5. If needed, **repeat steps 1-4** with updated input values in `channels_settings.txt` file to get more accurate k-term.
