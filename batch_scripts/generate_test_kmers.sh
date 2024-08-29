#!/bin/bash
for k in 5 7 9 11 13 15 17 19 25
do
    mkdir -p output/test/"$k"mer_indices
    python PORTEKfind.py input/HIV/prepped/MD.fasta output/test --k "$k" --group MD &
    python PORTEKfind.py input/HIV/prepped/N.fasta output/test --k "$k" --group N &
    python PORTEKfind.py input/HIV/prepped/O.fasta output/test --k "$k" --group O &
done