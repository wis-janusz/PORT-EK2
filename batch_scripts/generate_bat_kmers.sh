#!/bin/bash
for k in 15
do
    mkdir -p output/bat/"$k"mer_indices
    python PORTEKfind.py input/bat/bat.fasta output/bat/ --k "$k" --group bat --header_format ncbi&
    python PORTEKfind.py input/bat/EPI_SET_240422qm.fasta output/bat/ --k "$k" --group hum --header_format gisaid&  
done