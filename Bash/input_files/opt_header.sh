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

# for file in *; doAdd commentMore actions
#     [[ "$file" == "$(basename "$0")" ]] && continue
#     [[ ! -f "$file" ]] && continue
#     sed -i '1i# opt b3lyp/6-31g' "$file"
# done