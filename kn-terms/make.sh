#!/bin/bash

rm ./Control/*
find . -name '*.out' -type f -exec rm {} \;
find . -name '*__?' -type f -exec rm {} \;
rm CO2_COOLING_* CO2_Fdown_* CO2_Fup_*
 
gfortran PT_LW_main_TOT.f90 MOD_INI_PT_KD.f90 PT_Dflux_forKD.f90 PT_K_forKD.f90 PT_Planck_separate.f90 PT_Slon_AVT_SEPAR.f90 -o getnkd.out
gfortran V-D_tau_SEP_Q.f90 V-fromTAUtoFLUX.f90 V-PLANCK_oper.f90 -o compare_tot.out