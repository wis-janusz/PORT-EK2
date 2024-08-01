#!/bin/bash
for k in 5 7 9 11 13
do
    mkdir -p output/HIV/scenario3/"$k"mer_indices
    python PORTEKfind2.py input/HIV/prepped/MA.fasta output/HIV/scenario3 --k "$k" --group MA &
    python PORTEKfind2.py input/HIV/prepped/MB.fasta output/HIV/scenario3 --k "$k" --group MB &
    python PORTEKfind2.py input/HIV/prepped/MC.fasta output/HIV/scenario3 --k "$k" --group MC &
    python PORTEKfind2.py input/HIV/prepped/MD.fasta output/HIV/scenario3 --k "$k" --group MD &
    python PORTEKfind2.py input/HIV/prepped/Mrest.fasta output/HIV/scenario3 --k "$k" --group Mrest &
    python PORTEKfind2.py input/HIV/prepped/CRF.fasta output/HIV/scenario3 --k "$k" --group CRF &
done