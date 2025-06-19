#!/bin/bash
for file in *; do
    [[ "$file" == "$(basename "$0")" ]] && continue
    [[ ! -f "$file" ]] && continue
    { 
        echo "# opt b3lyp/6-31g"
        echo ""
        tail -n +4 "$file"
    } > "$file.tmp" && mv "$file.tmp" "$file"
done