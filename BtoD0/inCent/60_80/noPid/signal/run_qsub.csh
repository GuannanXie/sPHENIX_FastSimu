#!/bin/tcsh
set nevents=${nevents1}
set particle=${particle1}
set runIndex=${runIndex1}

set dir=$PWD
starver SL16d  # this star version should keep same with you complie by ".L toyMcBtoD.C++
echo $STAR
echo "out dir=$dir/out"

cd $dir
root4star -b -q runBtoD.C\(${nevents},\"${dir}\/out\/${particle}_${runIndex}\",\"${particle}\"\)
#root4star -b -q toyMcBtoD.C++\(${nevents},\"\/global\/homes\/x\/xlchen\/work\/Run14\/B_to_D\/FastSimulation\/out\/${particle}_${runIndex}\",\"${particle}\"\)
