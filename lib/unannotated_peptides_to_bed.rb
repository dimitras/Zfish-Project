# USAGE: ruby unannotated_peptides_to_bed.rb ../out/unannotated_mrnas_refseqs.txt ../out/mascot_unannotated_peptides.txt ../out/unannotated_proteins.fa ../out/unannotated_peptides.bed

#!/usr/bin/env ruby
require 'mrna'
require 'mrna_parser'
require 'peptide'
require 'peptide_parser'
require 'fap'


puts "Initializing (Getting command arguments)"

ucsc_mrnas_refseq = ARGV[0] # ../out/unannotated_mrnas_refseqs.txt
mascot_peptides = ARGV[1] # ../out/mascot_unannotated_peptides.txt
ucsc_protein_refseq = ARGV[2] # ../out/unannotated_proteins.fa
peptides_bed = ARGV[3] # unannotated peptides in bed format for ucsc # /out/unannotated_peptides.bed

ucsc_mrnas_refseq_mrnap = MrnaParser.open(ucsc_mrnas_refseq)
mascot_peptides_pep = PeptideParser.open(mascot_peptides)
ucsc_protein_refseq_fap = FastaParser.open(ucsc_protein_refseq) 
peptides_bedp = File.open(peptides_bed,"w")


puts "Create bed file from peptides"
bed_entries = {}
peptides_bedp.puts "track name=\"Unannotated Peptides\" description=\"Unannotated Peptides\" visibility=2 itemRgb=\"On\""
mascot_peptides_pep.each do |peptide|
 	peptide.mrna = ucsc_mrnas_refseq_mrnap.mrna(peptide.prot_acc)
	protein = ucsc_protein_refseq_fap.entry_by_id(peptide.prot_acc)
	if (((protein.length + 1) * 3) - peptide.mrna.cds_length) > 3 # +1 because mrna stop codon is not translated to protein
		next
	end
	entry = peptide.to_bed
	# checking for duplicates bed entries
	if !bed_entries.has_key?(entry)
		peptides_bedp.puts entry
		bed_entries[entry] = nil
	end
end
peptides_bedp.close

