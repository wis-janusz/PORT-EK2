#!/bin/bash
for k in 5 7 9 11 13 15 17 19 25
do
    mkdir -p output/deer/"$k"mer_indices
    python PORTEKfind.py input/deer/EPI_SET_240422va.fasta output/deer/ --k "$k" --group deer --header_format gisaid&
    python PORTEKfind.py input/deer/EPI_SET_240422rw.fasta output/deer/ --k "$k" --group humearly --header_format gisaid&  
    python PORTEKfind.py input/deer/EPI_SET_240422qc.fasta output/deer/ --k "$k" --group humlate --header_format gisaid&  
done