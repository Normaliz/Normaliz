#! /bin/sh
set -e
# Limit number of threads for normaliz
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
cd examples
for a in *.py; do
    echo "#################### $a ####################"
    python $a
done
