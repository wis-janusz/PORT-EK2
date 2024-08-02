#!/bin/bash
for k in 15
do
    mkdir -p output/deer/"$k"mer_indices
    python PORTEKfind.py input/deer/EPI_SET_240422va.fasta output/deer/ --k "$k" --group deer --header_format gisaid&
    python PORTEKfind.py input/deer/EPI_SET_240422rw.fasta output/deer/ --k "$k" --group humearly --header_format gisaid&  
    python PORTEKfind.py input/deer/EPI_SET_240422qc.fasta output/deer/ --k "$k" --group humlate --header_format gisaid&  
done