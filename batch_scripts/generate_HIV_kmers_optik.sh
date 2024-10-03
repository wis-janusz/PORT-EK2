#!/bin/bash
for k in 5 7 9 11 13
do
    mkdir -p output/HIV/optik/"$k"mer_indices
    python PORTEKfind.py input/HIV/prepped/M.fasta output/HIV/optik --k "$k" --group M &
    python PORTEKfind.py input/HIV/prepped/nonM.fasta output/HIV/optik --k "$k" --group nonM &  
done