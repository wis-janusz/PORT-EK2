#!/bin/bash
project_dir='/home/labadmin/repos/PORT-EK 2/projects/test'
for k in 5 7 9 11 13 15 17 19 25
do
    mkdir -p "$project_dir"/input/"$k"mer_indices
    python PORTEKfind.py "$project_dir"/input/MD.fasta "$project_dir"/input --k "$k" --group MD &
    python PORTEKfind.py "$project_dir"/input/N.fasta "$project_dir"/input --k "$k" --group N &
    python PORTEKfind.py "$project_dir"/input/O.fasta "$project_dir"/input --k "$k" --group O &
done