import argparse
import json

from Bio import SeqIO
from datetime import datetime


parser = argparse.ArgumentParser(
    description="Find all kmers in all sequences from a fasta file. Save count and optionally position matrix."
)
parser.add_argument("in_file", help="path to the .fasta file with input sequences")
parser.add_argument(
    "out_dir",
    help="absolute path to the output directory, which must already exist",
)
parser.add_argument("--k", help="lenght of kmers to find", type=int)
parser.add_argument("--group", help="name of the sample group", type=str)



def _find_kmers_df(seq_list: list, k, out_dir, group):
    seq_done = 1
    for seq in seq_list:
        kmers_dict = {}
        for i in range(len(seq.seq) - k):
            kmer = seq.seq[i : i + k]
            if all(nuc in ["A", "T", "G", "C"] for nuc in kmer):
                if kmer in kmers_dict.keys():
                    kmers_dict[kmer].append(i + 1)
                else:
                    kmers_dict[str(kmer)] = [i + 1]
        
        with open(f"{out_dir}/{group}_{seq.id}.json", mode="w") as out_file:
            json.dump(kmers_dict, out_file)

        print(
            f"Completed {seq.id}, {seq_done} of {len(seq_list)} sequences.",
            sep=" ",
            end="\r",
            flush=True,
        )
        seq_done += 1


def main():
    args = parser.parse_args()
    print(f"Start, {datetime.now()}")

    seq_list = list(SeqIO.parse(args.in_file, format="fasta"))
    for seq in seq_list:
        if "EPI_ISL" in seq.id:
            seq.id = seq.id.split("|")[1]
        else:
            seq.id = seq.id.split("|")[0][:-1]
        seq.seq = seq.seq.upper()

    _find_kmers_df(seq_list, args.k, args.out_dir, args.group)
    print(f"\nDone, {datetime.now()}")


if __name__ == "__main__":
    main()

