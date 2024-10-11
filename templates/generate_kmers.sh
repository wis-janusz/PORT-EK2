#!/bin/bash
project_dir='' #Fill with path to your project directory.
for k in 15 #Change to a space separated list of desired k-mer lengths.
do
    mkdir -p "$project_dir"/input/"$k"mer_indices
    #Uncomment the following lines and change input file and group names as desired.
    #If you want more that 2 groups just copy and paste more lines below.
    #Type 'python ./portek/portek_get_kmers.py -h' for more info about k-mer finding options

    # python ./portek/portek_get_kmers.py "$project_dir"/input/GROUP1.fasta "$project_dir"/input --k "$k" --group GROUP1 &
    # python ./portek/portek_get_kmers.py "$project_dir"/input/GROUP2.fasta "$project_dir"/input --k "$k" --group GROUP2 &


done