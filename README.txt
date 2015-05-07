Identify V and J genes using conserved sites and apply calculated v ties to V identification.

Main file:
identification.m

Sub functions:
calculateList.m
identificationJ.m
identificationV.m

Databases:
germlines.mat
HHV.fasta

How to use:
identification(file_name), file_name = name of fasta file that contains sequences you want to identify

Outputs:
noFoundJ.txt: order number of sequences without J
noFoundV.txt: order number of sequences with J but without V
intactCDR3.txt: order number of sequences with both J and V
Jgene.txt: order number and identified J gene groups
Vgene.txt: order number and identified V genes, number of mismatches, number of aligned nucleotides and indels (true/false)
CDR3.txt: order number and CDR3 nucleotides and amino acids
numbers.txt: order number and positions of J and V anchor
Vties_{database}_{aligned length}_{mutation frequency}_{v ties threshold}.txt: v tie pairs calculated based on {database}, {aligned length}, {mutation frequency} and {v ties threshold}