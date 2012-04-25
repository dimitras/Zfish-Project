#!/usr/bin/env ruby

# USAGE: ruby scripts/find_unexplained_peptides.rb data/1d-16dec-zfish-full-db.dat results/mrnas_full_db_lineformat.txt results/unexplained_peptides_of_unannotated_proteins.bed results/unique_peptides_of_unannotated_proteins.bed results/nonunique_peptides_of_unannotated_proteins.bed

$LOAD_PATH << './lib'

require 'peptide'
require 'peptide_parser'
require 'mrna_parser'


full_db_dat_file = ARGV[0]
mrnas_lineformat_file = ARGV[1]
unexplained_peptides_of_unannotated_proteins_bed_file = ARGV[2]
unique_peptides_of_unannotated_proteins_bed_file = ARGV[3]
multiple_peptides_of_unannotated_proteins_bed_file = ARGV[4]

full_db_dat = PeptideParser.open(full_db_dat_file)
mrnas_lineformat_mrnap = MrnaParser.open(mrnas_lineformat_file)

stronger_pep_for_unann = {}
max_score_for_ann = {}
peptide_times_found = Hash.new { |h,k| h[k] = 0 }
# peptide_protein_combination = Hash.new { |h,k| h[k] = 0 }

full_db_dat.each do |peptide|
	protein_is_annotated = true
	
	if peptide.prot_acc.split("-")[3] == 'unann'
		protein_is_annotated = false
	end
	
	if protein_is_annotated
		if !max_score_for_ann.has_key?(peptide.peptide) || (peptide.score > max_score_for_ann[peptide.peptide])
			max_score_for_ann[peptide.peptide] = peptide.score
		end
	else
		peptide_times_found[peptide.peptide] += 1 # pos na kano loop gia kathe protein? thelo na vro oles tis proteines pou paratireitai to peptidio
		if !stronger_pep_for_unann.has_key?(peptide.prot_acc) || (peptide.score > stronger_pep_for_unann[peptide.prot_acc].score)
			peptide.mrna = mrnas_lineformat_mrnap.mrna_by_prot_accno(peptide.prot_acc)
			stronger_pep_for_unann[peptide.prot_acc] = peptide
		end
	end
end
full_db_dat.rewind

unique_peps_for_unann_bed_output = File.open(unique_peptides_of_unannotated_proteins_bed_file,"w")
nonunique_peps_for_unann_bed_output = File.open(multiple_peptides_of_unannotated_proteins_bed_file,"w")
unique_peps_for_unann_bed_output.puts "track name=\"Unique Peptides on Unannotated Proteins\" description=\"Unique Peptides on Unannotated Proteins\" visibility=2 itemRgb=\"On\" useScore=1"
nonunique_peps_for_unann_bed_output.puts "track name=\"Non-Unique Peptides on Unannotated Proteins\" description=\"Non-Unique Peptides on Unannotated Proteins\" visibility=2 itemRgb=\"On\" useScore=1"

full_db_dat.each do |peptide|
	protein_is_annotated = true
	if peptide.prot_acc.split("-")[3] == 'unann'
		protein_is_annotated = false
	end
	
	if protein_is_annotated == false
		peptide.mrna = mrnas_lineformat_mrnap.mrna_by_prot_accno(peptide.prot_acc)
		if peptide_times_found[peptide.peptide] == 1
			unique_peps_for_unann_bed_output.puts peptide.to_bed
		else
			nonunique_peps_for_unann_bed_output.puts peptide.to_bed
		end
	end
end
unique_peps_for_unann_bed_output.close
nonunique_peps_for_unann_bed_output.close


unexplained_peptides_bed_output = File.open(unexplained_peptides_of_unannotated_proteins_bed_file,"w")
unexplained_peptides_bed_output.puts "track name=\"Unexplained highest-scored Peptides on Unannotated Proteins\" description=\"Unexplained highest-scored Peptides on Unannotated Proteins\" visibility=2 itemRgb=\"On\" useScore=1"
stronger_pep_for_unann.each_value do |peptide|
	if !max_score_for_ann.has_key?(peptide.peptide) || (max_score_for_ann[peptide.peptide] < peptide.score)
		if peptide.genomic_starts.length >= 2
			unexplained_peptides_bed_output.puts peptide.to_bed
		end
	end	
end
unexplained_peptides_bed_output.close

