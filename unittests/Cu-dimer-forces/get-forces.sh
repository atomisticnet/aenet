#!/bin/bash

for f in dimer*.xsf.out
do
  d=$(egrep "^ Cu" $f | awk '{d+=(-1.0)^(NR)*$2} END{print d}')
  F=$(egrep "^ Cu" $f | awk '{fx=$5} END{print fx}')
  E=$(awk '/Cohesive energy/{print $4}' $f)
  printf "%4.1f  %10.6f  %10.6f\n" $d $E $F
done
