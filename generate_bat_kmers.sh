#!/bin/bash
for k in 5 7 11 15 19 25
do
    mkdir -p output/bat/"$k"mer_indices
    python PORTEKfind.py input/bat/bat.fasta output/bat/"$k"mer_indices/ --k "$k" --group bat &
    python PORTEKfind.py input/bat/EPI_SET_240422qm.fasta output/bat/"$k"mer_indices/ --k "$k" --group hum &
done