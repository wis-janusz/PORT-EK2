#!/bin/bash
for k in 5 7 11 15 19 25
do
    mkdir -p output/deer/"$k"mer_indices
    python PORTEKfind.py input/deer/EPI_SET_240422va.fasta output/deer/"$k"mer_indices/ --k "$k" --group deer &
    python PORTEKfind.py input/deer/EPI_SET_240422rw.fasta output/deer/"$k"mer_indices/ --k "$k" --group humearly &  
    python PORTEKfind.py input/deer/EPI_SET_240422qc.fasta output/deer/"$k"mer_indices/ --k "$k" --group humlate &  
done