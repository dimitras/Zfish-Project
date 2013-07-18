# Alternative mRNA transcripts in Zebrafish

## Aim
To discover novel gene models and peptides from alternative splicing that haven’t been annotated yet

## The idea
- combine genomics and proteomics 
- explore for novel gene models
- identify novel proteins
- search for novel peptide evidence in novel transcripts

## Basic Steps
- RNA-Seq experiments
- map RNA-Seq reads against genome and transcriptome to identify junctions (RUM)
- translate mRNAs to proteins to create a novel protein database from RNA-Seq
- mass spectometry experiments
- identify novel proteins (mascot)
- identify and justify novel peptides (custom pipeline)

## Data
- post fertilized developing zebrafish embryos 
- for both RNA-Seq and MS2
	- 1 day experiment with 4 replicates
	- 3 days experiment with 3 replicates
	- 5 days experiment with 5 replicates

- mrnas from UCSC (lineformat)
- proteins sequences from UCSC (fasta)
- genome sequences (fasta)
- RNA-Seq short reads
- MS2 data that contain peptide spectrum matches (mgf)
- high quality junctions from RUM (bed)
- protein/peptide matches from mascot

## Steps
- create protein database for mascot searches
- create novel proteins
- identify peptides of the novel proteins 
- double check by searching for novel peptides in a database that contains both annotated and novel proteins
- visualize the candidates in UCSC genome browser

## Pipeline

- clean the raw db files (only before creating the db on mascot server)
ruby scripts/clean_files.rb data/mrnas_ucsc_lineformat.txt data/proteins_ucsc_db.fa data/danRer7.fa results/mrnas_ucsc_lineformat_clean.txt results/proteins_ucsc_db_clean.fa

- find novel proteins
ruby scripts/find_unannotated_proteins.rb results/mrnas_ucsc_lineformat_clean.txt data/1d-30oct-zfish-ucsc-db.dat results/proteins_ucsc_db_clean.fa data/junctions_high-quality.bed data/danRer7.fa results/1d-30oct/peptides-1d-30oct.bed results/1d-30oct/mrnas_unannotated.bed results/1d-30oct/mrnas_unannotated_seqs.fa results/1d-30oct/mrnas_unannotated_lineformat.txt results/1d-30oct/proteins_unannotated.fa

- run mascot for novel proteins database 

- create bed file with the peptides of the novel proteins
ruby scripts/peptides_to_bed.rb results/1d-30oct/mrnas_unannotated_lineformat.txt data/1d-30oct-zfish-unannotated-db.dat results/1d-30oct/proteins_unannotated.fa results/1d-30oct/peptides_on_unannotated_proteins-1d-30oct.bed

- create the db on mascot server

- join the db files (mrnas_lineformat and proteins)
cat results/mrnas_ucsc_lineformat_clean.txt > results/mrnas_full_db_lineformat.txt
grep -v "#" results/mrnas_unannotated_lineformat.txt >> results/mrnas_full_db_lineformat.txt
cat results/proteins_ucsc_db_clean.fa results/proteins_unannotated.fa >> results/proteins_full_db.fa

- run mascot for joined db

- find novel peptides with evidence
ruby scripts/find_unexplained_peptides.rb data/1d-30oct-zfish-full-db.dat results/mrnas_full_db_lineformat.txt data/mrnas_ucsc_lineformat.txt results/1d-30oct/unexplained_peptides_of_unannotated_proteins-1d-30oct.bed results/1d-30oct/unique_peptides_of_unannotated_proteins-1d-30oct.bed results/1d-30oct/nonunique_peptides_of_unannotated_proteins-1d-30oct.bed results/1d-30oct/unexplained_peptides_of_unannotated_proteins_per_spectra-1d-30oct.bed


## TODO list
- review the peptide candidates
- add more unit tests
- run all datasets
- TPP for data validation
- visualization of genomic annotation information with DAS
