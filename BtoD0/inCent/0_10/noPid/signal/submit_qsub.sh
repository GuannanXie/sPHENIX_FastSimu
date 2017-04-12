#!/bin/bash
nevents=$1
particle=$2 #
startIndex=$3
endIndex=$4
dir=$PWD

[ -d log ] || mkdir log
[ -d csh ] || mkdir csh
[ -d out ] || mkdir out

for((i=$startIndex;i<$endIndex;i++))
do
  csh="csh/${particle}_${i}.csh"
  [ -e $csh ] && echo "-------$csh exists!!!!"
  [ -e $csh ] && exit
  cat run_qsub.csh >$csh
  macro="./$csh"
  err="${dir}/log/${particle}_${i}.err"
  out="${dir}/log/${particle}_${i}.out"
  qsub -hard -l scratchfree=500,h_vmem=4G -o $out -e $err -v nevents1=$nevents,particle1=$particle,runIndex1=$i $macro
  echo "err = $err"
  echo "out = $out"
  echo "macro = $macro"
  echo "-------submitting ${particle}_${i}.csh"
  echo ""
done
