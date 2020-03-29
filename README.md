### BLAST_cytB_16S
#### This repo contains scripts I use to make Genus/species calls from raw ampliconSeq data.

##### blastn can give you very odd results if you use the full database. For example, you might get Homo sapiens and Homo neanderthalensis hits for a human sample. The custom_mito_DB.py script makes a custom database with species know to be in CA plus specific invasive species that you can add in.

##### The ampliconseq_genotyper.py is a script that identifies unique amplicons that contain specific primer sequences that you input. Then it blasts the paired-end reads to a custom mito database. It outputs an excel file with 2 tabs. One with calls for each target and one that includes all unique reads for each target with at least 2 identical copies. This tab can be used to examine individual calls in more detail.


