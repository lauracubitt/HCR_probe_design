# HCR_probe_design
Generates probes for multiplex experiment from multiple genes for third generation HCR

Requires:
1. Needs RNA fasta files to be downloaded from NCBI as format: genename_datasets/ncbi_dataset/data/rna.fna
2. Requires genelist.csv as input with names of genes, format: needs to have a Header in cell A1 and then list below 
3. Change species on line 26 as required - this is for Blast function only 

Note: this will run a lot quicker if you remove the Blast function (which will check that the probes don't match any other genes in the genome). To do this comment out lines 148-198.
