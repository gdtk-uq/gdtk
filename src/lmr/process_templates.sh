#!/bin/bash

# Simple template processor
# Usage: ./process_templates.sh <file>
# Replaces {{command}} with the output of command

file="$1"

while grep -q '{{.*}}' "$file"; do
    template=$(grep -o '{{[^}]*}}' "$file" | head -n1)
    cmd=$(echo "$template" | sed 's/{{//;s/}}//')
    result=$(eval "$cmd")
    sed -i "s|$template|$result|" "$file"
done
