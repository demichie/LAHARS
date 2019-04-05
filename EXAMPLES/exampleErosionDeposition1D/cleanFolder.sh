#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of IMEX-SfloW2d"

file="exampleErosionDeposition*"

if [ -f $file ] ; then
    rm $file
fi

file="IMEX_SfloW2D.inp"

if [ -f $file ] ; then
    rm $file
fi

file="topography_dem.asc"

if [ -f $file ] ; then
    rm $file
fi


