#!/bin/bash

echo "filename,triplet_eV,singlet_eV" > excitation_energies.csv

for file in *.log; do
    triplet=$(grep -m 1 "Excited State" "$file" | grep "Triplet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Triplet-A") print $(i+1)}')
    singlet=$(grep "Excited State" "$file" | grep "Singlet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Singlet-A") print $(i+1)}' | head -n 1)
    echo "$file,$triplet,$singlet" >> excitation_energies.csv
done