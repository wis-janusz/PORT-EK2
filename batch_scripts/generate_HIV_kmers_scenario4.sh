#!/bin/bash
for k in 11
do
    mkdir -p output/HIV/scenario4/"$k"mer_indices
    python PORTEKfind.py input/HIV/prepped/MA.fasta output/HIV/scenario4 --k "$k" --group MA &
    python PORTEKfind.py input/HIV/prepped/MB.fasta output/HIV/scenario4 --k "$k" --group MB &
    python PORTEKfind.py input/HIV/prepped/MC.fasta output/HIV/scenario4 --k "$k" --group MC &
    python PORTEKfind.py input/HIV/prepped/MD.fasta output/HIV/scenario4 --k "$k" --group MD &
    python PORTEKfind.py input/HIV/prepped/Mrest.fasta output/HIV/scenario4 --k "$k" --group Mrest &
    python PORTEKfind.py input/HIV/prepped/N.fasta output/HIV/scenario4 --k "$k" --group N &
    python PORTEKfind.py input/HIV/prepped/O.fasta output/HIV/scenario4 --k "$k" --group O &
    python PORTEKfind.py input/HIV/prepped/P.fasta output/HIV/scenario4 --k "$k" --group P &
done