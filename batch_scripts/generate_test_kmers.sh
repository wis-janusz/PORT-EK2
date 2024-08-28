#!/bin/bash
for k in 5 7 9 11 13 15 17 19 25
do
    mkdir -p output/test/"$k"mer_indices
    python PORTEKfind.py input/HIV/prepped/MA.fasta output/HIV/scenario1 --k "$k" --group MA &
    python PORTEKfind.py input/HIV/prepped/MC.fasta output/HIV/scenario1 --k "$k" --group MC &
    python PORTEKfind.py input/HIV/prepped/MD.fasta output/HIV/scenario1 --k "$k" --group MD &
    python PORTEKfind.py input/HIV/prepped/N.fasta output/HIV/scenario1 --k "$k" --group N &
    python PORTEKfind.py input/HIV/prepped/O.fasta output/HIV/scenario1 --k "$k" --group O &
    python PORTEKfind.py input/HIV/prepped/P.fasta output/HIV/scenario1 --k "$k" --group P & 
done