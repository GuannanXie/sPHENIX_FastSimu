#!/bin/tcsh
set nevents=${nevents1}
set particle=${particle1}
set runIndex=${runIndex1}

cd /global/homes/x/xlchen/rnc/sPhenix/noPid
starver SL16d  # this star version should keep same with you complie by ".L toyMcBtoD.C++
root4star -b -q runBtoD.C\(${nevents},\"\/global\/homes\/x\/xlchen\/rnc\/sPhenix\/noPid\/out\/${particle}_${runIndex}\",\"${particle}\"\)
#root4star -b -q toyMcBtoD.C++\(${nevents},\"\/global\/homes\/x\/xlchen\/work\/Run14\/B_to_D\/FastSimulation\/out\/${particle}_${runIndex}\",\"${particle}\"\)
