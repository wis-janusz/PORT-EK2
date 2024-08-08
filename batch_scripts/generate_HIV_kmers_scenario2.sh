#!/bin/bash
for k in 5 7 9 11 13
do
    mkdir -p output/HIV/scenario2/"$k"mer_indices
    python PORTEKfind.py input/HIV/prepped/M.fasta output/HIV/scenario2 --k "$k" --group M &
    python PORTEKfind.py input/HIV/prepped/N.fasta output/HIV/scenario2 --k "$k" --group N &
    python PORTEKfind.py input/HIV/prepped/O.fasta output/HIV/scenario2 --k "$k" --group O &
    python PORTEKfind.py input/HIV/prepped/P.fasta output/HIV/scenario2 --k "$k" --group P &

done