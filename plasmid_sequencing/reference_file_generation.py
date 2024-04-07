#!/usr/bin/env python
# coding: utf-8

import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob

BioBrick_sequences_dict = {}
with open("2018_Distribution.fasta", "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Using the sequence identifier (header) as the key
        # and the sequence as the value
        parts = record.id.split('_')

        # Concatenate the first two parts
        BBa_number = '_'.join(parts[:2])
        BioBrick_sequences_dict[BBa_number] = record

BioBrick_backbones_dict = {}
for BioBrick_backbone_path in glob.glob("BioBrick_backbones/*.gbk"):
    with open(BioBrick_backbone_path, "r") as genbank_file:
        for record in SeqIO.parse(genbank_file, "genbank"):
            print(BioBrick_backbone_path + " " + record.id)
            # Using the sequence identifier (header) as the key
            # and the sequence as the value
            BioBrick_backbones_dict[record.id] = record


with open("Burden-o-meter-sequencing.csv", "r", encoding="utf-8-sig") as csv_file:
    reader = csv.DictReader(csv_file)
    for row in reader:

        print(row)
        if row['included'] == 'FALSE':
            continue

        if row['backbone'] == 'NA':
            continue

        print("Creating sequence for " + row['sample.name'])

        backbone_sequence = BioBrick_backbones_dict[row['backbone']]

        #print(str(backbone_sequence) + str(BioBrick_sequences_dict[row['biobrick']]))
        sequence_object = Seq(backbone_sequence.seq + BioBrick_sequences_dict[row['biobrick']].seq)
        sequence_record = SeqRecord(sequence_object, id=row['sample.name'], description="")
        sequence_record.annotations["molecule_type"] = "DNA"

        #print(sequence_record)
        # add annotations back from backbone
        for feature in backbone_sequence.features:
            sequence_record.features.append(feature)

        with open("references/" + row['sample.name'] + ".gbk", "w") as output_file:
            SeqIO.write([sequence_record], output_file, "genbank")
        




