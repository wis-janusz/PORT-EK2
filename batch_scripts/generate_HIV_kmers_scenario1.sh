#!/bin/bash
for k in 5 7 9 11 13 15 17 19 25
do
    mkdir -p output/HIV/scenario1/"$k"mer_indices
    python PORTEKfind.py input/HIV/prepped/M.fasta output/HIV/scenario1 --k "$k" --group M &
    python PORTEKfind.py input/HIV/prepped/nonM.fasta output/HIV/scenario1 --k "$k" --group nonM &  
done