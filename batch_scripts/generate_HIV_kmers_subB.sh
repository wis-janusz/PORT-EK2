#!/bin/bash
for k in 11
do
    mkdir -p output/HIV/subB/"$k"mer_indices
    python PORTEKfind.py input/HIV/prepped/B_main.fasta output/HIV/subB --k "$k" --group Bmain &
    python PORTEKfind.py input/HIV/prepped/B_out.fasta output/HIV/subB --k "$k" --group Bout &
done