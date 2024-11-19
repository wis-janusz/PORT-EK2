# PORT-EK
Pathogen Origin Recogition Tool (using) Enriched K-mers
Developmental version 2.0 - WORK IN PROGRESS

A tool for identification of genomic variants of virues that arise in connection with host change.
Based on k-mer counting, does not require MSA.
Can highlight changes independent of viral philogeny.
Can be used to find group-enriched variants for other sequences and groups than viruses and hosts.

Installation:
1. Clone this repo,
2. Create a Python virtual enviroment,
3. Install required packages using `pip install -r requirements.txt` ,

Quick usage guide:
1. Execute `python PORTEKrun.py new $project_directory` to create a new project in project_directory.
2. Edit the project configuration file and copy fasta files with input seqences to $project_directory/input.
3. Execute `python PORTEKrun.py find_k $project_directory --max_k $k` where k is the maximum k-mer length of k-mer you want to test, to test different length of k-mers and choose the best one for your data.
4. Execute `python PORTEKrun.py enriched $project_directory --k $k ` to get k-mers enriched in your sample groups, where k is the length of the k-mers.
5. --TBD--


