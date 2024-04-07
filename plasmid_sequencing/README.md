## Overview

This directory contains scripts and reference sequences used to analyze sequencing
data generated to validate and verify BioBrick plasmids.

`Burden-o-meter-sequencing.csv` contains metadata about the sequenced samples.

The file `data.csv` defines the reference files and read files used by
to predict mutations with _breseq_ in each sample. Reads are available from the SRA.
Reference sequences are provided here under `referecnes`. How they are generated is
described below.

## Reference File Generation

The script `download_distribution.py` downloads the BioBrick part sequences from the iGEM 
Registry website. The descriptions of these parts are in the files in `plates`. It creates 
the output file `2018_Distribution.fasta`.

The script `reference_file_generation.py` uses these part sequences and the vector sequences
in `BioBrick_backbones` to generate the full sequences of each plasmid that was analyzed.
Be aware that minor differences due to different assembly RFCs and CDS versus noncoding
parts are not accounted for by this script! The output is in `references`.
