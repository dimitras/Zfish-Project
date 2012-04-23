#!/usr/bin/env ruby

# USAGE: ruby scripts/find_unexplained_peptides.rb data/1d-16dec-zfish-full-db.dat results/mrnas_full_db_lineformat.txt results/unexplained_peptides_of_unannotated_proteins.bed results/unique_peptides_of_unannotated_proteins.bed results/nonunique_peptides_of_unannotated_proteins.bed

$LOAD_PATH << './lib'

require 'rubygems'
require 'mascot/dat'
require 'mascot/dat/peptides'
require 'mascot/dat/psm'
require 'peptide'
require 'mrna_parser'


full_db_dat_file = ARGV[0]
mrnas_lineformat_file = ARGV[1]
unexplained_peptides_of_unannotated_proteins_bed_file = ARGV[2]
unique_peptides_of_unannotated_proteins_bed_file = ARGV[3]
multiple_peptides_of_unannotated_proteins_bed_file = ARGV[4]

@full_db_dat = Mascot::DAT.open(full_db_dat_file, true)
mrnas_lineformat_mrnap = MrnaParser.open(mrnas_lineformat_file)

@peptides = @full_db_dat.peptides
stronger_pep_for_unann = {}
max_score_for_ann = {}
peptide_times_found = {}


@peptides.each do |psm|
	peptide_times_found[psm.pep] = 0
	psm.proteins.each do |protein_info|
		protein_accno = protein_info[0]
		start = protein_info[2]
		stop = protein_info[3]
		multiplicity = protein_info[4]
		protein_is_annotated = true
		
		if protein_accno.split("-")[3] == 'unann'
			protein_is_annotated = false
		end
		
		if protein_is_annotated
			if !max_score_for_ann.has_key?(psm.pep) || (psm.score > max_score_for_ann[psm.pep])
				max_score_for_ann[psm.pep] = psm.score
			end
		else
			peptide_times_found[psm.pep] += 1
			if !stronger_pep_for_unann.has_key?(protein_accno) || (psm.score > stronger_pep_for_unann[protein_accno].score)
				peptide = Peptide.new(protein_accno, psm.pep, start-1, stop-1, multiplicity, psm.score)
				peptide.mrna = mrnas_lineformat_mrnap.mrna_by_prot_accno(peptide.prot_acc)
				stronger_pep_for_unann[protein_accno] = peptide
			end
		end
	end
end


unique_peps_for_unann_bed_output = File.open(unique_peptides_of_unannotated_proteins_bed_file,"w")
nonunique_peps_for_unann_bed_output = File.open(multiple_peptides_of_unannotated_proteins_bed_file,"w")
@peptides.rewind
@peptides.each do |psm|
	psm.proteins.each do |protein_info|
		protein_accno = protein_info[0]
		start = protein_info[2]
		stop = protein_info[3]
		multiplicity = protein_info[4]
		protein_is_annotated = true
		
		if protein_accno.split("-")[3] == 'unann'
			protein_is_annotated = false
		end
		
		if protein_is_annotated == false
			peptide = Peptide.new(protein_accno, psm.pep, start-1, stop-1, multiplicity, psm.score)
			peptide.mrna = mrnas_lineformat_mrnap.mrna_by_prot_accno(peptide.prot_acc)
			if peptide_times_found[psm.pep] == 1
				unique_peps_for_unann_bed_output.puts peptide.to_bed
			else
				nonunique_peps_for_unann_bed_output.puts peptide.to_bed
			end
		end
	end
end
unique_peps_for_unann_bed_output.close
nonunique_peps_for_unann_bed_output.close


unexplained_peptides_bed_output = File.open(unexplained_peptides_of_unannotated_proteins_bed_file,"w")
stronger_pep_for_unann.each_value do |peptide|
	if !max_score_for_ann.has_key?(peptide) || (max_score_for_ann[peptide] < peptide.score)
		unexplained_peptides_bed_output.puts peptide.to_bed
	end	
end
unexplained_peptides_bed_output.close





