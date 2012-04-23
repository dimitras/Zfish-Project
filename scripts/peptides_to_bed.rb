# USAGE for unannotated: ruby scripts/peptides_to_bed.rb results/mrnas_unannotated_lineformat.txt results/mascot_peptides.txt results/proteins_unannotated.fa results/peptides_on_unannotated_proteins.bed

# USAGE for full db: ruby scripts/peptides_to_bed.rb results/mrnas_full_db_lineformat.txt results/mascot_peptides_full_db.txt results/proteins_full_db.fa results/peptides_full_db.bed

#!/usr/bin/env ruby

$LOAD_PATH << './lib'

require 'mrna'
require 'mrna_parser'
require 'peptide'
require 'peptide_parser'
require 'fap'


puts "Initializing (Getting command arguments)"

ucsc_mrnas_lineformat_file = ARGV[0]
mascot_peptides_file = ARGV[1]
proteins_seqs_file = ARGV[2]
peptides_bed_file = ARGV[3]

ucsc_mrnas_lineformat_mrnap = MrnaParser.open(ucsc_mrnas_lineformat_file)
mascot_peptides_pep = PeptideParser.open(mascot_peptides_file)
proteins_seqs_fap = FastaParser.open(proteins_seqs_file) 
peptides_bedp = File.open(peptides_bed_file,"w")


puts "Create bed file from peptides on unannotated proteins"
bed_entries = {}
peptides_bedp.puts "track name=\"Peptides on Unannotated Proteins\" description=\"Peptides on Unannotated Proteins\" visibility=2 itemRgb=\"On\" useScore=1"
mascot_peptides_pep.each do |peptide|
 	peptide.mrna = ucsc_mrnas_lineformat_mrnap.mrna_by_prot_accno(peptide.prot_acc)
	protein = proteins_seqs_fap.entry_by_id(peptide.prot_acc)
	# if (((protein.length + 1) * 3) - peptide.mrna.cds_length) > 3 # +1 because mrna stop codon is not translated to protein
	# 	next
	# end
	entry = peptide.to_bed()
	# checking for duplicates bed entries
	if !bed_entries.has_key?(entry)
		peptides_bedp.puts entry
		bed_entries[entry] = nil
	end
end
peptides_bedp.close
