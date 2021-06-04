#!/usr/bin/env bash

#array=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#array=(2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 20 21 22 X)
array=( X )

for i in "${array[@]}"
do
   echo ${i}
   array2=(/Sid/tdvorkina/monomers/CenSat/IvanAnnotByHOR/chr${i}*.bed)
   for f in "${array2[@]}"
   do
       echo $f
       python3 build_monoconsensus.py /Sid/tdvorkina/monomers/CenSat/Asm/chr${i}.fasta $f /Sid/tdvorkina/monomers/CenSat/HORdec/cen${i} --hors /Sid/tdvorkina/monomers/CenSat/IvanHORs/HOR${i}.tsv
   done
done
