#!/bin/bash

DIR="./Control"
if [ -d $DIR ]
then
	if [ "$(ls -A $DIR)" ]; then
        rm ./Control/*
    fi
fi

find . -name '*.out' -type f -exec rm {} \;
find . -name '*__?' -type f -exec rm {} \;

if [ -f "_COOLING" ]; then
    rm _COOLING 
fi

if [ -f "_Fdown_" ]; then
    rm _Fdown_ 
fi

if [ -f "_Fup_" ]; then
    rm _Fup_ 
fi
 
gfortran mod_pt_kd.f90 main_1step_pt_kd.f90 k_pt_kd.f90 dflux_pt_kd.f90 planck_separate_kd.f90 -o get1kd.out
gfortran main_compare.f90 tau_to_flux_compare.f90 planck_compare.f90 -o compare.out 