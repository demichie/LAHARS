#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of IMEX-SfloW2d"

for file in exampleLahar*
do
    if [ -f $file ] ; then
        rm $file
    fi
done

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


