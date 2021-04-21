#!/bin/bash

declare -a K1
K1=(0 2 4 10 14 20 40 60 80 200 400)
declare -a K2
K2=(0 5 10 25 35 50 100 150 200 500 1000)

for G in `seq 1 1`
do
for MOD in "scf"
do
  FILE=$(printf "data/g%02d.dat" $G)
  ./kmst -f $FILE -m $MOD -k ${K1[$G]}
done
done
