#!/bin/bash

for MONOMER_ID in 1 2 3 .. 12
do 
    python3 build_matching_heatmap.py $MONOMAR_ID&
    echo "Script on $MONOMER_ID strted"
done

echo "WAIT for python scripts"
