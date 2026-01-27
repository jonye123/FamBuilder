#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

def parse_orthogroups(orthogroups_file):
    orthogroups = {}
    with open(orthogroups_file, 'r') as f:
        # skip header line
        header = next(f)
        if header.strip().split()[0] == "HOG":
            offset = 3
        elif header.strip().split()[0] == "Orthogroup":
            offset = 1
        else:
            raise ValueError(f"Unexpected header format: {header.strip()}")
        for line in f:
            parts = line.strip().split()
            og = parts[0]
            members = []
            for part in parts[offset:]:
                members.extend(part.split(','))
            orthogroups[og] = [member.strip() for member in members if member.strip()]
    return orthogroups

def extract_sequences(orthogroups, fasta_dir, output_dir):
    # Create a dictionary to store all sequences
    all_sequences = {}
    
    # Read all FASTA files in the input directory
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fa') or filename.endswith('.fasta') or filename.endswith('.pep') or filename.endswith('.faa'):
            file_path = os.path.join(fasta_dir, filename)
            for record in SeqIO.parse(file_path, 'fasta'):
                all_sequences[record.id] = record

    os.makedirs(output_dir, exist_ok=True)

    # Write sequences for each orthogroup
    for og, members in orthogroups.items():
        output_file = os.path.join(output_dir, f"{og}.fasta")
        with open(output_file, 'w') as out_f:
            for seq_id in members:
                if seq_id in all_sequences:
                    SeqIO.write(all_sequences[seq_id], out_f, 'fasta')
                else:
                    print(f"Warning: Sequence {seq_id} not found in FASTA files")

def main():
    parser = argparse.ArgumentParser(description="Create FASTA files for each orthofamily from OrthoFinder output")
    parser.add_argument("orthogroups_file", help="Path to the Orthogroups.tsv file, either HOG or OG, requires header")
    parser.add_argument("fasta_dir", help="Directory containing input FASTA files")
    parser.add_argument("output_dir", help="Directory to save output FASTA files")
    
    args = parser.parse_args()

    orthogroups = parse_orthogroups(args.orthogroups_file)
    extract_sequences(orthogroups, args.fasta_dir, args.output_dir)
    print(f"FASTA files for each orthofamily have been created in {args.output_dir}")

if __name__ == "__main__":
    main()