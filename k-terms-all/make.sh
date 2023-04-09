#!/bin/bash

find . -name '*.out' -delete
find . -name '*__?' -delete
find . -name '*_*[0-9]_' -delete
find . -name 'k_coef*' -delete
find . -name '*.*[0-9]_' -delete
find . -name '*_[0-9]_*.dat' -delete
find . -name '*.Ven' -delete
find . -name '*VenCO2*' -delete
find . -name 'Fl_N-*' -delete

gfortran PT_LW_main_TOT.f90 MOD_INI_PT_KD.f90 PT_Dflux_forKD.f90 PT_K_forKD.f90 PT_Planck_separate.f90 PT_Slon_AVT_SEPAR.f90 -o separate.out
gfortran V-D_tau_SEP_Q.f90 V-fromTAUtoFLUX.f90 V-PLANCK_oper.f90 -o getnkd.out