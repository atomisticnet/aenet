#!/bin/bash

inp=dimer.xsf-in

x0=2.0
x1=10.0

y=0.0
z=0.0

for x in $(seq $x0 0.1 $x1)
do
  fname="dimer-x$(printf "%04.1f" $x).xsf"
  sed "s/%X%/$x/g;s/%Y%/$y/g;s/%Z%/$z/g" $inp > $fname
done

exit 0
