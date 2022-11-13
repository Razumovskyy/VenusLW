#!/bin/bash

rm ./Control/*
find . -name '*.out' -type f -exec rm {} \;
find . -name '*__?' -type f -exec rm {} \;
rm _COOLING _Fdown_ _Fup_
 
gfortran main_1step_pt_kd.f90 mod_pt_kd.f90 k_pt_kd.f90 dflux_pt_kd.f90 planck_separate_kd.f90 -o get1kd.out
gfortran main_compare.f90 tau_to_flux_compare.f90 planck_compare.f90 -o compare.out 