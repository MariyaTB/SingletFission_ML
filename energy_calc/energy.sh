
#!/bin/bash

# csv="excitation_energies.csv"

# if [ ! -f "$csv" ]; then
#     echo "filename,triplet1_eV,triplet2_eV,singlet1_eV" > "$csv"
#     echo "" >> "$csv"
# fi

# for file in *.log; do
#     if ! grep -q "^$file," "$csv"; then
#         triplet1=$(grep "Excited State" "$file" | grep "Triplet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Triplet-A") print $(i+1)}' | sed -n 1p)
#         triplet2=$(grep "Excited State" "$file" | grep "Triplet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Triplet-A") print $(i+1)}' | sed -n 2p)
#         singlet1=$(grep "Excited State" "$file" | grep "Singlet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Singlet-A") print $(i+1)}' | sed -n 1p)
#         echo "$file,$triplet1,$triplet2,$singlet1" >> "$csv"
#     fi
# done

#!/bin/bash

csv="excitation_energies.csv"

if [ ! -f "$csv" ]; then
    echo "filename,triplet1_eV,triplet2_eV,singlet1_eV,2T1_minus_S1,2T1_minus_T2" > "$csv"
    echo "" >> "$csv"
fi

for file in *.log; do
    if ! grep -q "^$file," "$csv"; then
        triplet1=$(grep "Excited State" "$file" | grep "Triplet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Triplet-A") print $(i+1)}' | sed -n 1p)
        triplet2=$(grep "Excited State" "$file" | grep "Triplet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Triplet-A") print $(i+1)}' | sed -n 2p)
        singlet1=$(grep "Excited State" "$file" | grep "Singlet-A" | awk '{for(i=1;i<=NF;i++) if($i=="Singlet-A") print $(i+1)}' | sed -n 1p)
        twot1_minus_s1=$(echo "scale=6; 2*${triplet1:-0} - ${singlet1:-0}" | bc)
        twot1_minus_t2=$(echo "scale=6; 2*${triplet1:-0} - ${triplet2:-0}" | bc)
        echo "$file,$triplet1,$triplet2,$singlet1,$twot1_minus_s1,$twot1_minus_t2" >> "$csv"
    fi
done