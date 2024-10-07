#!/bin/bash
for k in 7
do
    mkdir -p output/test/"$k"mer_indices
    python PORTEKfind.py input/test/A.fasta output/test/ --k "$k" --group A 
    
done