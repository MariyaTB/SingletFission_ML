#!/bin/bash

csv="excitation_energies.csv"

if [ ! -f "$csv" ]; then
    echo "filename,triplet_eV,singlet_eV" > "$csv"
fi

for file in *.log; do
    if ! grep -q "^$file," "$csv"; then
        triplet=$(grep -m 1 "Excited State" "$file" | grep "Triplet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Triplet-A") print $(i+1)}')
        triplet=$(grep -m 2 "Excited State" "$file" | grep "Triplet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Triplet-A") print $(i+1)}')
        singlet=$(grep "Excited State" "$file" | grep "Singlet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Singlet-A") print $(i+1)}' | head -n 1)
        echo "$file,$triplet1,$triplet1,$singlet" >> "$csv"
    fi
done