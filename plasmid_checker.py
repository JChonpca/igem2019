#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Thu Mar  28 2019

python script for plasmids

@author: ded
"""

import argparse
import os
import re
from collections import  defaultdict
from Bio import SeqIO
from itertools import combinations
from io import StringIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import datetime


import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("-a", "--adapter_dir", help="directory containing trimmomatic logs", default="Adapter_trim_logs")
parser.add_argument("-s", "--spade_completeion_log", help="completeion log checking for Spades completeing sucessfully", default="Spade_completeion_log_summary.txt")
parser.add_argument("-cd", "--contig_dir", help="directory of contigs in fasta format.", default="assemblies/All_contigs/")
parser.add_argument("-r", "--reference_dir", help="directory containing reference files", default="02_Ref")
parser.add_argument("-ss", "--sign_up_sheet", help="generate command files based on signup-sheet read in", default="SignUpSheet.tsv")
parser.add_argument("-l", "--local", help="flag to run on local machine and not run every aspect. used primarily for testing", action='store_true')
parser.add_argument("--commands_only", help="only generate commands files for adapter, spades, breseq based on input file", action='store_true')
parser.add_argument("--stats_only", help="only generate stats files for adapter, and spades output", action='store_true')

# show help if no arguments passed
# if len(sys.argv) < 2:
#     parser.print_help()
#     sys.exit(1)

args = parser.parse_args()

print("\nScript executed at: " + str(datetime.datetime.now()) + "\nWith the following command line options:")
print("\t" + "\n\t".join(map(str, [str(arg) + "\t" + str(getattr(args, arg)) for arg in vars(args)])))  # log all input options



#Test for script failure conditions
fail_to_print = []
if not os.path.exists("Raw_Reads") or len([x for x in os.listdir("Raw_Reads") if x.endswith(".gz")]) == 0:
    fail_to_print.extend(["Reads should be in `Raw_Reads` directory. Try the following suggested commands from inside the directory containing the raw reads and rerun this script", "\tmkdir Raw_Reads", "\tmv *.gz Raw_Reads/"])

if not os.path.exists(args.sign_up_sheet):
    fail_to_print.append("Sign up sheet not located. Specified by default or at command line to be: `%s`. Please make sure this file exists or specify an alternative at the command line using `-ss` or `--sign_up_sheet` flags") % args.sign_up_sheet

if len(fail_to_print) > 0:
    print("\n".join(map(str, ["\n\n\n!!!!!!\nSCRIPT FAIL!!\n!!!!!!\n\n\n"] + fail_to_print + [os.getcwd()])))
    sys.exit()


#generate lists for warnings and completed tasks and failed tasks  #TODO redo as single dictionary see final print area for more details
warnings_to_print = []
completed_tasks_to_print = []
failed_tasks_to_print = []
warned_tasks_to_print = []





adapter_dict = {}
Contig_dict = {}

def adapter_trimming_stats():
    print("adapter trimming stats")
    for f in os.listdir(args.adapter_dir):

        if f.endswith(".log"):
            name = re.sub(".log$", "", f)
            assert name not in adapter_dict and name not in Contig_dict, "duplicated file name:%s" % f
            Contig_dict[name] = {'R1R2': [], 'R1': [], 'R2': [], 'P1P2U1U2': [], 'P1P2': [], 'P1': [], 'P2': [], 'P1P2U1': [], 'P1P2U2': [],  'P1U1': [], 'P2U2': []}

            with open(args.adapter_dir + "/" + f, "r") as f_in:
                for line in f_in:
                    line=line.rstrip()
                    if re.match("Input Read Pairs", line):
                        line = line.split(" ")
                        # [input, both, forward, reverse, drop]
                        adapter_dict[name] = [line[3], line[6], line[11], line[16], line[19]]

    print("\t".join(map(str, ["Sample", "Raw read pairs", "Both", "Forward only", "Reverse only", "Drop"])))
    for name in sorted(adapter_dict):
        print("\t".join(map(str, [name] + adapter_dict[name])))


def spades_stats():
    # identify failed plasmid assembiles. NOTE this is failed to complete assembly pipeline, NOT failed to assemble any contigs
    print("calculating spade stats")
    with open(args.spade_completeion_log, "r") as sl_in:

        for line in sl_in:
            line=line.rstrip()
            if re.search(".log:1$", line):
                continue
            assert re.search(".log:0$", line), "Spade completion log lists sample that somehow completed multiple times. This does not make sense. \n%s" % line
            assert line.count("/") == 1 and line.count(".log:0") == 1, "Unexpected assembly log format log files expected to be stored in a single directory, and end with '.log'. this file should also have ended in :0 as currently looking at/for failed assemblies.\n%s" % line
            name = re.sub(".log:0","", re.sub(".*\/", "", line))  # strip directory info from line
            name = name.split(".")
            Contig_dict[name[0]]["".join(name[1:])] = "NA"

def contig_read_in():
    #Read contigs into dictionary
    print("reading contigs in")
    for f in os.listdir(args.contig_dir):

        if f.endswith(".contigs.fasta"):
            base = re.sub(".contigs.fasta$", "", f)
            assert len(base.split(".")) == 2, "unable to identify name from contig file: %s\nName expected to be between '.' and file name to end with .contigs.fasta" % f
            name = base.split(".")[1]

            assert base.split(".")[0].split("_")[1] == "reads", "unable to idtnify assembly conditions from contig file: %s\nFormat expected to be 'Raw|Trim_reads_' followed by R|P|U 1|2 each separated by '_' marks" % f
            assembly = "".join(base.split(".")[0].split("_")[2:])

            assert Contig_dict[name][assembly] is not "NA", "Contig file found for assembly previously identified as a failed assembly. This is different than an assembly which produces no contig but successfully completes assembly pipeline. %s" % f

            with open(args.contig_dir + "/" + f, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    Contig_dict[name][assembly].append((record.name, record.seq))

    all_single_contig_dict = {}
    print("\t".join(map(str, ["Sample", 'R1R2', 'R1', 'R2', 'P1P2U1U2', 'P1P2', 'P1', 'P2', 'P1P2U1', 'P1P2U2', 'P1U1', 'P2U2', "tot Fail", "tot None", "tot >1", "tot 1"])))
    for sample in sorted(Contig_dict):
        if args.local:
            # select one of the following to run or not run this block if running/testing locally
            continue
            #pass
        to_print = [sample]
        all_single_contig_dict[sample] = []
        for assembly_type in ['R1R2', 'R1', 'R2', 'P1P2U1U2', 'P1P2', 'P1', 'P2', 'P1P2U1', 'P1P2U2', 'P1U1', 'P2U2']:
            if Contig_dict[sample][assembly_type] == "NA":
                to_print.append("NA")
            else:
                to_print.append(len(Contig_dict[sample][assembly_type]))
                if len(Contig_dict[sample][assembly_type]) == 1:
                    all_single_contig_dict[sample].append(Contig_dict[sample][assembly_type][0][1])
        to_print += [to_print[1:].count("NA"), to_print[1:].count(0), sum(i>1 for i in to_print[1:] if i is not "NA"), to_print[1:].count(1)]
        print("\t".join(map(str, to_print)))

    #check for identical
    testing = 0
    for sample in sorted(Contig_dict):
        if args.local:
            # select one of the following to run or not run this block if running/testing locally
            continue
            #pass

        if sample != "JCM-10" : continue
        print(sample)


        subject_ref_file = "Reference_sequences/JCM-10.fa"





        sample_score_matrix = []
        if len(all_single_contig_dict[sample]) > 1:
            if len(set(all_single_contig_dict[sample])) > 1:
                num_contigs = 0
                fasta_files = []
                for seq in set(all_single_contig_dict[sample]):
                    num_contigs +=1
                    fasta_files.append("contig_check_temp_files/" + sample + "_" + str(num_contigs) + ".fa")
                    if not os.path.exists("contig_check_temp_files/" + sample + "_" + str(num_contigs) + ".fa"):
                        with open("contig_check_temp_files/" + sample + "_" + str(num_contigs) + ".fa", "w") as output:
                            if num_contigs == 2:
                                print(">" + sample + "\n" + seq +seq, file=output)
                            else:
                                print(">" + sample + "\n" + seq, file=output)
                    print(sample + "_" + str(num_contigs) + ".fa" + "\t" + str(len(seq)))


                for c in combinations(fasta_files, 2):
                    print("\t".join(map(str, [sample, c[0], c[1]])))

                    blast_output = NcbiblastnCommandline(query=c[0], subject=c[1], outfmt=5)()[0]
                    blast_result = NCBIXML.read(StringIO(blast_output))

                    perfects = 0
                    for aln in blast_result.alignments:
                        print(len(aln.hsps), "alignments")

                        for hsp in aln.hsps:
                            if (hsp.align_length - hsp.identities == 0) and (hsp.gaps == 0):
                                #print hsp.identities
                                print("following is perfect")
                                perfects += 1
                            print('****Alignment****')
                            print('sequence:', aln.title)
                            print('length:', hsp.align_length)
                            print('e value:', hsp.expect)
                            print(hsp.query)
                            print(hsp.match)
                            print(hsp.sbjct)
                            print(hsp.identities, "matches")
                    print(perfects, "perfects")
                testing +=1
                if testing >= 1:
                    sys.exit()

# fail pairwise2 would have been inside contig_read_in
# skip = True  #temp for testing. this block takes the longest by far
# for sample in sorted(Contig_dict):
#     if skip:
#         continue
#     sample_score_matrix = []
#     if len(all_single_contig_dict[sample]) > 1:
#         if len(set(all_single_contig_dict[sample])) > 1:
#             for c in combinations(set(all_single_contig_dict[sample]), 2):
#                 if len(c[0]) >= len(c[1]):
#                     sample_score_matrix.append(max(([pairwise2.align.globalxx(c[0] + c[0], c[1], score_only=True),
#                                                      pairwise2.align.globalxx(c[0] + c[0], c[1].reverse_complement(), score_only=True)])) / len(c[0]))
#                 else:
#                     sample_score_matrix.append(max(([pairwise2.align.globalxx(c[1] + c[1], c[0], score_only=True),
#                                                      pairwise2.align.globalxx(c[1] + c[1], c[0].reverse_complement(), score_only=True)])) / len(c[1]))
#             print "\t".join(map(str, [sample, sum(sample_score_matrix) / float(len(sample_score_matrix)), sample_score_matrix]))
#         else:
#             print "\t".join(map(str, [sample, 1.0, "ident string"]))
#     else:
#         if len(all_single_contig_dict[sample]) == 1:
#             print "\t".join(map(str, [sample, 1.0, "single positive"]))
#         else:
#             print "\t".join(map(str, [sample, "NA", "No single contig"]))

#print("starting ref checks")
# #check contigs against references
# for ref in os.listdir(args.reference_dir):
#     if args.local:
#         # select one of the following to run or not run this block if running/testing locally
#         continue
#         #pass
#     if ref.endswith('.fa') or ref.endswith('.fasta'):
#         name = re.sub(".fasta$", "", ref)
#         name = re.sub(".fa$", "", name)
#         assert name in Contig_dict, "unmatched reference file name %s" % ref
#         ref_compare_score = []
#         if len(set(all_single_contig_dict[name])) >= 1:
#             best_contig = "NA"
#             with open(args.reference_dir + "/" + ref, "rU") as handle:
#                 fasta_count = 0
#                 for record in SeqIO.parse(handle, "fasta"):
#                     fasta_count += 1
#                     for contig in set(all_single_contig_dict[name]):
#                         if len(contig) >= len(record.seq):
#                             ref_compare_score.append(max(([pairwise2.align.globalxx(contig + contig, record.seq, score_only=True), pairwise2.align.globalxx(contig +contig, record.seq.reverse_complement(), score_only=True)]))/ len(contig))
#                         else:
#                             ref_compare_score.append(max(([pairwise2.align.globalxx(record.seq + record.seq, contig, score_only=True),
#                                                            pairwise2.align.globalxx(record.seq + record.seq, contig.reverse_complement(), score_only=True)])) / len(record.seq))
#                         if ref_compare_score[-1] == max(ref_compare_score):
#                             best_contig = contig
#                 assert best_contig != "NA", "No contig scored best? %s %s" % (name, ref_compare_score)
#                 print max(ref_compare_score)
#                 if len(best_contig) >= len(record.seq):
#                     if pairwise2.align.globalxx(best_contig + best_contig, record.seq, score_only=True) / len(best_contig) > pairwise2.align.globalxx(best_contig + best_contig, record.seq.reverse_complement(), score_only=True) / len(best_contig):
#                         print pairwise2.format_alignment(*pairwise2.align.globalxx(best_contig + best_contig, record.seq, one_alignment_only=True)[0])
#                     else:
#                         print pairwise2.format_alignment(*pairwise2.align.globalxx(best_contig + best_contig, record.seq.reverse_complement(), one_alignment_only=True)[0])
#                 else:
#                     if pairwise2.align.globalxx(record.seq + record.seq, best_contig, score_only=True) / len(record.seq) > pairwise2.align.globalxx(record.seq + record.seq, best_contig.reverse_complement(), score_only=True) / len(record.seq):
#                         print pairwise2.format_alignment(*pairwise2.align.globalxx(record.seq + record.seq, best_contig, one_alignment_only=True)[0])
#                     else:
#                         print pairwise2.format_alignment(*pairwise2.align.globalxx(record.seq + record.seq, best_contig.reverse_complement(), one_alignment_only=True)[0])
#
#                 if fasta_count > 1:
#                     print "\n\n\n!!!!!!\n\n\n!!!!!!!\nWARNING multiple fasta sequences found in file: %s\n\n\n" % ref
#
#                 print "\t".join(map(str, [name, sum(ref_compare_score)/ float(len(ref_compare_score)), ref_compare_score]))
#         else:
#             print "\t".join(map(str, [name, "No single contig to compare"]))
#
#                 #for assembly_type in Contig_dict[]


def spades_commands_gen():
    #check if this will run:
    if len(os.listdir("Trim_Reads/")) == 0:
        print("Spades commands can not be run yet, but commands generated. commands can be run upon completion of adapter trimming commands")

    print("Spades command gen started")

    assembly_dirs_check = ["assemblies", "assemblies/Trim_reads_P1_P2_U1_U2", "assemblies/Trim_reads_P1_P2", "assemblies/Trim_reads_P1", "assemblies/Trim_reads_P2", "assemblies/Trim_reads_P1_P2_U1", "assemblies/Trim_reads_P1_P2_U2", "assemblies/Trim_reads_P1_U1", "assemblies/Trim_reads_P2_U2", "Spades_Logs"]
    for dir_check in assembly_dirs_check:
        if not os.path.exists(dir_check):
            os.makedirs(dir_check)


    with open(args.sign_up_sheet, "r") as SUS, open("spades_commands", "w") as spades_command_out:
        Header_search = True
        for line in SUS:
            line = line.rstrip()
            if Header_search:
                assert re.match("^Researcher\tSample_ID\tReference_file", line), "header line missing in %s" % line
                Header_search = False
                continue
            if re.match("Your name", line):
                continue

            line = line.split("\t")

            Name = line[0].replace(" ", "_") + "." + line[1]
            P1 = Name + "_1P.fastq.gz"
            P2 = Name + "_2P.fastq.gz"
            U1 = Name + "_1U.fastq.gz"
            U2 = Name + "_2U.fastq.gz"

            print("plasmidspades.py -o assemblies/Trim_reads_P1_P2_U1_U2/%s -1 %s -2 %s -s %s -s %s >& Spades_Logs/%s.P1.P2.U1.U2.log" %(Name, P1, P2, U1, U2, Name), file=spades_command_out)
            print("plasmidspades.py -o assemblies/Trim_reads_P1_P2/%s -1 %s -2 %s >& Spades_Logs/%s.P1.P2.log" % (Name, P1, P2, Name), file=spades_command_out)
            print("plasmidspades.py -o assemblies/Trim_reads_P1/%s -s %s >& Spades_Logs/%s.P1.log" % (Name, P1, Name), file=spades_command_out)
            print("plasmidspades.py -o assemblies/Trim_reads_P2/%s -s %s >& Spades_Logs/%s.P2.log" % (Name, P2, Name), file=spades_command_out)
            print("plasmidspades.py -o assemblies/Trim_reads_P1_P2_U1/%s -1 %s -2 %s -s %s >& Spades_Logs/%s.P1.P2.U1.log" % (Name, P1, P2, U1, Name), file=spades_command_out)
            print("plasmidspades.py -o assemblies/Trim_reads_P1_P2_U2/%s -1 %s -2 %s -s %s >& Spades_Logs/%s.P1.P2.U2.log" % (Name, P1, P2, U2, Name), file=spades_command_out)
            print("plasmidspades.py -o assemblies/Trim_reads_P1_U1/%s -s %s -s %s >& Spades_Logs/%s.P1.U1.log" % (Name, P1, U1, Name), file=spades_command_out)
            print("plasmidspades.py -o assemblies/Trim_reads_P2_U2/%s -s %s -s %s >& Spades_Logs/%s.P2.U2.log" % (Name, P2, U2, Name), file=spades_command_out)




def SignUpSheet_read_in():
    with open(args.sign_up_sheet, "r") as SUS:
        all_raw_reads = [x for x in os.listdir("Raw_Reads") if x.endswith(".gz") and not re.match("^" + "Undetermined_S0_L001_", x)]
        used_read_file = []
        # for raw_read_file in os.listdir("Raw_Reads"):
        #     if raw_read_file.endswith(".gz"):
        #         all_raw_reads.append(raw_read_file)
        assert len(all_raw_reads) > 0, "No reads were found in `Raw_Reads` move gzipped raw read files to that directory. Possible suggested command:\n\tmv *.gz Raw_Reads/\nPossible script failure. This should have been caught in failure checks"
        all_ids = []

        Header_search = True
        for line in SUS:
            line = line.rstrip()
            if Header_search:
                assert re.match("^Researcher\tSample_ID\tReference_file",line), "header line missing in %s" %line
                Header_search = False
                continue
            if re.match("Your name", line):
                continue
            line = line.split("\t")
            researcher = line[0].replace(" ", "_")
            sample_ID = line[1]
            reference_file = line[2]
            assert sample_ID not in master_dict[researcher], "duplicated sampleID for single researcher. %s\nmaster_dict:\n%s" % (line, list(master_dict[researcher].keys()))
            assert sample_ID not in all_ids, "duplicated sampleID for different researchers %s\nall ids:\n%s" % (line, all_ids)

            master_dict[researcher][sample_ID] = {"Ref": reference_file, "Raw_reads": []}
            all_ids.append(sample_ID)

            for raw_read in all_raw_reads:
                if re.match("^" + sample_ID +"_", raw_read):
                    master_dict[researcher][sample_ID]["Raw_reads"].append(raw_read)
                    used_read_file.append(raw_read)

        unused_reads = [x for x in all_raw_reads if x not in used_read_file]
        if len(unused_reads) > 0:
            warnings_to_print.append("The following read file exists without a sample ID and will not be trimmed.\n\t" + "\n\t".join(sorted(unused_reads)))

        samples_missing_reads = [" from ".join([s_ID, res]) for res in master_dict for s_ID in master_dict[res] if len(master_dict[res][s_ID]["Raw_reads"]) == 0]
        if len(samples_missing_reads) > 0:
            warnings_to_print.append("No read files were identified for the following samples, and samples will not be analyzed.\n\t" + "\n\t".join(sorted(samples_missing_reads)))

        all_refs = list(set([master_dict[res][s_ID]["Ref"] for res in master_dict for s_ID in master_dict[res]]))
        for ref in all_refs:
            if os.path.exists("02_Ref/" + ref):
                pass
            else:
                if "signupsheet read in" not in warned_tasks_to_print:
                    warned_tasks_to_print.append("signupsheet read in")
                warnings_to_print.append("%s reference file missing. move file to 02_Ref/ directory" % ref)

            if ref.endswith(".gbk") or ref.endswith(".gb") or ref.endswith(".fa") or ref.endswith(".fasta"):
                pass
            else:
                if "signupsheet read in" not in warned_tasks_to_print:
                    warned_tasks_to_print.append("signupsheet read in")
                warnings_to_print.append("%s reference file is of unknown type. rename or verify file is correctly formated genbank or fasta file" % ref)

        if "signupsheet read in" not in warned_tasks_to_print:
            completed_tasks_to_print.append("signupsheet read in")


def adapter_trim_commands_gen():
    if not os.path.exists("Trim_Reads"): os.makedirs("Trim_Reads")
    if not os.path.exists("Adapter_trim_logs"): os.makedirs("Adapter_trim_logs")
    with open("adapter_trim_commands", "w") as adapter_command_out:
        for researcher in master_dict:
            for sample_ID in master_dict[researcher]:
                sample_reads = ["Raw_Reads/" + x for x in master_dict[researcher][sample_ID]["Raw_reads"]]
                sample_reads = " ".join(sample_reads)
                print("trimmomatic PE %s -baseout Trim_Reads/%s.fastq.gz ILLUMINACLIP:illumina_truseq_6bp_dual_barcode_PE.fasta:4:30:10 MINLEN:30 >& Adapter_trim_logs/%s.log" % (
                sample_reads, sample_ID, sample_ID), file=adapter_command_out)

    if not os.path.exists("illumina_truseq_6bp_dual_barcode_PE.fasta"):
        warnings_to_print.append(
            "adapter fasta file not found. run the following command to allow adapter trimming:\n\tcp /corral-repl/utexas/breseq/genomes/adapters/illumina_truseq_6bp_dual_barcode_PE.fasta .")
        warned_tasks_to_print.append("adapter command generation")
    else:
        completed_tasks_to_print.append("adapter command generation")


def breseq_commands_gen():
    if not os.path.exists("03_Output"): os.makedirs("03_Output")
    if not os.path.exists("04_Logs"): os.makedirs("04_Logs")

    with open("breseq_commands", "w") as breseq_command_out, open("gdtools_commands", "w") as gd_command_out:
        gdtools_compare_dict = defaultdict(lambda: defaultdict(list))
        for researcher in master_dict:

            if not os.path.exists("03_Output/" + researcher): os.makedirs("03_Output/" + researcher)

            for sample_ID in master_dict[researcher]:

                gdtools_compare_dict[researcher][master_dict[researcher][sample_ID]["Ref"]].append(sample_ID + ".gd")

                trim_read_files = " ".join(["Trim_Reads/" + sample_ID + x for x in ["_1P.fastq.gz", "_2P.fastq.gz", "_1U.fastq.gz"]])  # , "_2U.fastq.gz"  # eliminated as several are 0 and dont have quick fix yet

                print("breseq -o 03_Output/%s/%s -r 02_Ref/%s %s >& 04_Logs/%s.%s.log.txt"  % (researcher, sample_ID, master_dict[researcher][sample_ID]["Ref"], trim_read_files, researcher, sample_ID), file=breseq_command_out)

        #TODO work this into run on tacc
        print("The following commands can not be run as a script:", file=gd_command_out)
        for researcher in gdtools_compare_dict:
            print(researcher, file=gd_command_out)
            for ref in gdtools_compare_dict[researcher]:
                print("gdtools compare -r ../../../02_Ref/%s -o %s.%s.html %s" % (ref, researcher, ref, " ".join(sorted(gdtools_compare_dict[researcher][ref]))), file=gd_command_out)

    #TODO need stuff here about sucess/failure

    #Check if commands can run
    # if len(os.listdir("Trim_Reads/")) == 0:  # TODO change this to warn commands can't be run until after adapters complete rather than looking for reads
    #     print "Breseq commands can not be generated at this time as Trimmed reads folder is currently empty. Likely cause of this is adapter commands not having ben run yet. rerun script after adapters are trimmed"
    #     return
    #
    # print "breseq command gen started"
    # gdtools_compare_dict = {}
    # if not os.path.exists("03_Output_Trim"): os.makedirs("03_Output_Trim")
    # if not os.path.exists("04_Logs"): os.makedirs("04_Logs")
    #
    # with open(args.sign_up_sheet, "r") as f_in, open("breseq_commands", "w") as command_out, open("gdtools_commands", "w") as gd_command_out:
    #
    #     Header_search = True
    #     Samples_lacking_reads = []
    #     for line in f_in:
    #         line = line.rstrip()
    #
    #         if Header_search:
    #             assert re.match("^Researcher\tSample_ID\tReference_file",line), "header line missing in %s" % line
    #             Header_search = False
    #             continue
    #         if re.match("Your name", line):
    #             continue

            # line = line.split("\t")
            # Name = line[0].replace(" ", "_")
            # Base_sample_ID = line[1]
            # ref_file = line[2]
            #
            # if Name not in gdtools_compare_dict:
            #     gdtools_compare_dict[Name] = {}
            #
            # if ref_file not in gdtools_compare_dict[Name]:
            #     gdtools_compare_dict[Name][ref_file] = []
            #
            # gdtools_compare_dict[Name][ref_file].append(Base_sample_ID+".gd")
            #
            # trim_reads = ["Trim_Reads/" + x for x in os.listdir("Trim_Reads") if re.match("^" + Base_sample_ID + "_", x)]
            # if len(trim_reads) == 0:
            #     Samples_lacking_reads.append(Base_sample_ID)
            #     continue
            # trim_reads = " ".join(trim_reads)
            #
            # Output_dir = Name + Base_sample_ID
            # if not os.path.exists("03_Output_Trim/" + Output_dir):os.makedirs("03_Output_Trim/" + Output_dir)
            # print>>command_out, "breseq -o 03_Output_Trim/%s -r 02_Ref/%s %s >& 04_Logs/%s.%s.Trim.log.txt" % (Output_dir, ref_file, trim_reads, Name, Base_sample_ID)

        # print>>gd_command_out, "The following commands can not be run as a script:"
        # for researcher in gdtools_compare_dict:
        #     print>>gd_command_out, researcher
        #     for ref_entry in gdtools_compare_dict[researcher]:
        #         print>>gd_command_out, "gdtools compare -r ../../../02_Ref/%s -o %s.%s.html %s" % (ref_entry, researcher, ref_entry, " ".join(gdtools_compare_dict[researcher][ref_entry]))

        # if len(Samples_lacking_reads) > 0:
        #     print "The following samples did not have breseq commands generated as reads were not identified corresponding to the sample:"
        #     print "\t" + "\n\t".join(map(str, Samples_lacking_reads))


master_dict = defaultdict(lambda: defaultdict())



SignUpSheet_read_in()  # must be called to allow anything else to run
if args.commands_only:
    adapter_trim_commands_gen()
    breseq_commands_gen() # won't work until trimming done
    spades_commands_gen() # commands won't work until trimming done

if args.stats_only:
    pass



if args.local:
    #comment out parts you dont want to run
    #adapter_trim_commands_gen()
    #breseq_commands_gen() # won't work until trimming done
    #spades_commands_gen() # commands won't execute until trimming done
    #adapter_trimming_stats()
    #spades_stats()
    #contig_read_in() #requires adapter trimming stats to run
    pass

# else: # run all elements
#     adapter_trimming_stats()
#     breseq_commands_gen()
#     spades_stats()
#     contig_read_in()




#TODO rework as dictionary of tasks: d[task] == {} -- sucess ... d[fail/warn] = [failures/warnings] use asserts to make sure task isn't in d, and warnings last thing generated. fail should exit function always

if len(completed_tasks_to_print) > 0:
    print("The following tasks were successfully completed:")
    print("\t"+"\n\t".join(map(str, completed_tasks_to_print)))

if len(failed_tasks_to_print) > 0:
    print("The following tasks failed:")
    print("\t"+"\n\t".join(map(str, failed_tasks_to_print)))

if len(warned_tasks_to_print) > 0:
    assert len(warnings_to_print) > 0, "Task warned without warning being generated. Script problem. \n%s" % warned_tasks_to_print
    print("The following tasks generated warnings:")
    print("\t"+"\n\t".join(map(str, warned_tasks_to_print)))



if len(warnings_to_print) != 0:
    print("\n".join(map(str, ["\n\n\n!!!!!!\nSCRIPT WARNINGS!!\n!!!!!!\n\n\n"] + warnings_to_print)))
else:
    print("Script completed successfully without generating any warnings.")
