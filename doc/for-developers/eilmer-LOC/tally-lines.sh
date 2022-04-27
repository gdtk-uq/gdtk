#! /bin/bash
# tally-lines.sh
# copy this script into gdtk/src/ and run from there.

echo "D "
find . -name "*.d" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "C++ "
find . -name "*.cxx" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "hh"
find . -name "*.hh" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "Python"
find . -name "*.py" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "Lua -- total"
find . -name "*.lua" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "Lua -- excluding species-database & sample-data"
find . \( -name "species-database" -prune -o -name "sample-data" -prune \) -o -name "*.lua" -exec wc '{}' \;  | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "LaTex"
find . -name "*.tex" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "ReStructured Text"
find . -name "*.rst" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "ASCIIdoc"
find . -name "*.adoc" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

echo "Markdown"
find . -name "*.md" -exec wc '{}' \; | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'

