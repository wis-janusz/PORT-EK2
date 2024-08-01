import pathlib
from Bio import SeqIO

input_path = "input/HIV"
output_path = "input/HIV/prepped"

in_files = pathlib.Path(input_path).glob("*.fasta")

seq_lists_source = {}

for filename in in_files:
    seq_list = list(SeqIO.parse(filename, format="fasta"))
    for seq in seq_list:
        seq.seq = seq.seq.replace("-", "")
    seq_lists_source[filename.stem.replace("hiv-db_", "")] = seq_list

M_type_A_subtype = ("A", "A1", "A2", "A3", "A4", "A6", "A7", "A8")
M_type_rest_subtype = ("F1", "F2", "G", "H", "J", "K")
M_type = M_type_A_subtype + ("B", "C", "D") + M_type_rest_subtype
recombinants = ("CRF01_AE", "CRF02_AG")
non_M_type = ("N", "O", "P")

subtypes_groups_dict = {
    "M": M_type,
    "nonM": non_M_type,
    "N": ("N"),
    "O": ("O"),
    "P": ("P"),
    "MA": M_type_A_subtype,
    "MB": ("B"),
    "MC": ("C"),
    "MD":("D"),
    "Mrest":M_type_rest_subtype,
    "CRF":recombinants
}

seq_lists_target = {subtype:[] for subtype in subtypes_groups_dict.keys()}


for source_subtype, sequences in seq_lists_source.items():
    for target_subtype in subtypes_groups_dict.keys():
        if source_subtype in subtypes_groups_dict[target_subtype]:
            seq_lists_target[target_subtype].extend(sequences)

for target_subtype, sequences in seq_lists_target.items():
    SeqIO.write(sequences, f"{output_path}/{target_subtype}.fasta", format="fasta")
