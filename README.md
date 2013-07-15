# Alternative mRNA transcripts in Zebrafish

## Aim
To discover novel gene models that haven’t been annotated yet

## Solution
- combine RNA-Seq, alternative splicing and proteomics 
- discover alternative transcripts
- create a non-annotated protein database
- combine annotated with non-annotated protein db
- search for novel peptide evidence in novel transcripts

## Basic Steps
- RNA-Seq experiments
- map RNA-Seq reads against genome and transcriptome to identify junctions (RUM)
- mass spectometry experiments
- identify and characterize proteins from primary sequence databases (mascot)
- identify novel peptides (custom pipeline)

## Data
- 6 datasets (3 x 1day and 3 x 3days experiments)
- annotated mrnas from UCSC (lineformat)
- annotated  proteins sequences from UCSC (fasta)
- genome sequences (fasta)
- high quality junctions from RUM (bed)
- MS2 data that contain peak lists and experimental parameters (mgf)

## Algorithm
- use the clean database for mascot searches
- discover non-annotated proteins
- find the peptides of the non-annotated proteins
- produce new database joining annotated and non-annotated proteins
- search against them to identify novel peptides
- visualize the candidates in UCSC genome browser

## Pipeline

- clean the raw db files (only before creating the db on mascot server)
ruby scripts/clean_files.rb data/mrnas_ucsc_lineformat.txt data/proteins_ucsc_db.fa data/danRer7.fa results/mrnas_ucsc_lineformat_clean.txt results/proteins_ucsc_db_clean.fa

- find unannotated proteins
ruby scripts/find_unannotated_proteins.rb results/mrnas_ucsc_lineformat_clean.txt data/1d-30oct-zfish-ucsc-db.dat results/proteins_ucsc_db_clean.fa data/junctions_high-quality.bed data/danRer7.fa results/1d-30oct/peptides-1d-30oct.bed results/1d-30oct/mrnas_unannotated.bed results/1d-30oct/mrnas_unannotated_seqs.fa results/1d-30oct/mrnas_unannotated_lineformat.txt results/1d-30oct/proteins_unannotated.fa

- run mascot for unannotated db 

- create bed file with the peptides of the unannotated proteins
ruby scripts/peptides_to_bed.rb results/1d-30oct/mrnas_unannotated_lineformat.txt data/1d-30oct-zfish-unannotated-db.dat results/1d-30oct/proteins_unannotated.fa results/1d-30oct/peptides_on_unannotated_proteins-1d-30oct.bed

- create the db on mascot server

- join the db files (mrnas_lineformat and proteins)
cat results/mrnas_ucsc_lineformat_clean.txt > results/mrnas_full_db_lineformat.txt
grep -v "#" results/mrnas_unannotated_lineformat.txt >> results/mrnas_full_db_lineformat.txt
cat results/proteins_ucsc_db_clean.fa results/proteins_unannotated.fa >> results/proteins_full_db.fa

- run mascot for joined db

- find unexplained peptides with evidence
ruby scripts/find_unexplained_peptides.rb data/1d-30oct-zfish-full-db.dat results/mrnas_full_db_lineformat.txt data/mrnas_ucsc_lineformat.txt results/1d-30oct/unexplained_peptides_of_unannotated_proteins-1d-30oct.bed results/1d-30oct/unique_peptides_of_unannotated_proteins-1d-30oct.bed results/1d-30oct/nonunique_peptides_of_unannotated_proteins-1d-30oct.bed results/1d-30oct/unexplained_peptides_of_unannotated_proteins_per_spectra-1d-30oct.bed


## TODO list
- review the peptide candidates
- add more unit tests
- run all datasets
- TPP for data validation
- visualization of genomic annotation information with DAS
