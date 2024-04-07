## Samples

The file `data.csv` defines the reference files and read files used by
to predict mutations with _breseq_ in each sample. Reads are available from the SRA.

`JEB1203` is the burden monitor strain we used for our assays.

`JEB1216` is a colony pick of the AddGene strain #66073 submitted by the Ellis lab.

## Reference files

`NC_010473.1_Eco_K-12_DH10B.gbk` is the _E. coli_ chromosome reference sequence.

`AddGene-66073-pAH63-with-gfp-insert.gbk` is the sequence assembled from the unmapped reads for
this plasmid (whole sequence not available from Addgene).

## Conclusions

* Both strains have the expected insertion of the plasmid at the Lamba attP site.
    * This gives two junctions between the plasmid and the genome
    * Additionally, the circular plasmid junction is identified.
    * The coverage of the inserted plasmid matches the genome, so
      there is no evidence of tandem duplications that can sometimes occur.
* In addition, they have:
    * Two "mutations" which are assigning ambiguous bases to the correct bases (20,896, 142,348)
    * An intergenic mutation at 2,875,910
    * An intergenic +T insertion at 4,272,972
* JEB1216 has two ADDITIONAL mutations
    * A synonymous mutation in ECDH10B_RS19585 (DUF1198 domain‑containing protein) at 3,937,758
    * What appears to be an IS element insertion in the mdtL gene which may be associated
      with a change in antibiotic resistance (it is an MDR transporter). The duplicated target site
      is 3988356-3988364 (9 bases). The transposone is "IS4-like element ISVsa5 family transposase" This
      mutation may be polymorphic and present in ~70-90% of the sample.
