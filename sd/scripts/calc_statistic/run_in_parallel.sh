#!/bin/bash

for MONOMER_ID in 2 3 4 5 6 7 8 9 10 11 12
do 
    python3 build_matching_heatmap.py "$MONOMER_ID" &
    echo "python3 buils_matchng_heatmap.py $MONOMER_ID &"
    echo "Script on $MONOMER_ID strted"
done

echo "WAIT for python scripts"
