#!/bin/bash

 
HEADER="# td=(50-50,nstates=3) b3lyp/6-31g"

for log in *.log; do
    base="${log%.log}"
    gjf="${base}.gjf"

    obabel "$log" -o gjf -O "$gjf"

       tmpfile=$(mktemp)
    {
        echo "$HEADER"
        echo
        echo "${base} TD calculation"
        echo "0 1"
       
        cat "$gjf"
    } > "$tmpfile"
    mv "$tmpfile" "$gjf"
done