# USAGE: ruby scripts/find_unannotated_proteins.rb data/mrnas_ucsc_lineformat.txt results/mascot_peptides.txt data/proteins_ucsc_db.fa data/junctions_high-quality.bed data/danRer7.fa results/peptides.bed results/mrnas_unannotated.bed results/mrnas_unannotated_seqs.fa results/mrnas_unannotated_lineformat.txt results/proteins_unannotated.fa results/ucsc_mrnas_lineformat_cleanedup.txt

#!/usr/bin/env ruby

$LOAD_PATH << './lib'

require 'mrna'
require 'mrna_parser'
require 'peptide'
require 'peptide_parser'
require 'fap'
require 'junction'
require 'bed_parser'
require 'track'
require 'rubygems'
require 'bio'
require 'genome_fap'


# Input files
ucsc_mrnas_lineformat = ARGV[0] # mrnas_ucsc_lineformat.txt
mascot_peptides = ARGV[1] # mascot_peptides.txt
ucsc_protein_fasta_file = ARGV[2] # proteins_ucsc_db.fa
rum_junctions_hq_bed = ARGV[3] # junctions_high-quality.bed
danRer7_fasta_file = ARGV[4] # danRer7.fa

# Output files
peptides_bed = ARGV[5] # peptides.bed
unannotated_mrnas_bed = ARGV[6] # mrnas_unannotated.bed
unannotated_mrnas_seqs = ARGV[7] # mrnas_unannotated_seqs.fa
unannotated_mrnas_lineformat_file = ARGV[8] # mrnas_unannotated_lineformat.txt
unannotated_proteins = ARGV[9] # proteins_unannotated.fa


puts "Initializing (Reading input files)"
ucsc_mrnas_lineformat_mrnap = MrnaParser.open(ucsc_mrnas_lineformat)
mascot_peptides_pep = PeptideParser.open(mascot_peptides)
ucsc_proteins_fap = FastaParser.open(ucsc_protein_fasta_file) 
danRer7_fap = GenomeFastaParser.new(danRer7_fasta_file)
rum_junctions_hq_bedp = BedParser.open(rum_junctions_hq_bed)


puts "Create bed file from peptides"
peptides_bed_output = File.open(peptides_bed,"w")
peptides_bed_output.puts "track name=\"Peptides\" description=\"Peptides\" visibility=2 itemRgb=\"On\""
mascot_peptides_pep.each do |peptide|
 	peptide.mrna = ucsc_mrnas_lineformat_mrnap.mrna_by_prot_accno(peptide.prot_acc)
	protein = ucsc_proteins_fap.entry_by_id(peptide.prot_acc)
	# if peptide.mrna.cds_length != (protein.length + 1) * 3 then # +1 because mrna stop codon is not translated to protein
	# 	next
	# end
 	peptides_bed_output.puts peptide.to_bed
end
peptides_bed_output.close


puts "Read mrnas and create a track"
mrna_track = Track.new
ucsc_mrnas_lineformat_mrnap.each do |mrna|
	mrna_track.insert_genomic_region(mrna)
end


puts "Scanning junctions to find unannotated ones and create alternative gene models that would explain them"
unannotated_mrnas_seqs_output = File.open(unannotated_mrnas_seqs, "w")
unannotated_mrnas_seqs_fap = GenomeFastaParser.open(unannotated_mrnas_seqs)
unannotated_proteins_output = File.open(unannotated_proteins,"w")
unannotated_mrnas_bed_output = File.open(unannotated_mrnas_bed,"w")
unannotated_mrnas_lineformat_output = File.open(unannotated_mrnas_lineformat_file,"w")
unannotated_mrnas_bed_output.puts "track name=\"Alternative gene models\" description=\"Alternative gene models\" visibility=2 itemRgb=\"On\""
unannotated_mrnas_lineformat_output.puts "# danRer7.unannotated.name\tdanRer7.unannotated.chrom\tdanRer7.unannotated.strand\tdanRer7.unannotated.txStart\tdanRer7.unannotated.txEnd\tdanRer7.unannotated.cdsStart\tdanRer7.unannotated.cdsEnd\tdanRer7.unannotated.exonCount\tdanRer7.unannotated.exonStarts\tdanRer7.unannotated.exonEnds\tdanRer7.unannotated.name2\tdanRer7.unannotated.product\tdanRer7.unannotated.mrnaAcc\tdanRer7.unannotated.protAcc"
alt_mrna_checklist = {}
count = 0
cds_seq = ''
rum_junctions_hq_bedp.each do |junction|
# 	puts "Get mRNAs that overlap the junction"
	overlapping_mrnas = mrna_track.regions_that_contain(junction)
	
	
# 	puts "Check if the junction is annotated or not"
	annotated = false
	overlapping_mrnas.each do |mrna|
		junction_start_exon_idx = mrna.exon_index(junction.start)
		junction_stop_exon_idx = mrna.exon_index(junction.stop)
		if junction_start_exon_idx != nil && junction_stop_exon_idx != nil && (junction_stop_exon_idx - junction_start_exon_idx) == 1
# 			puts "Junction is annotated"
			annotated = true
			break
		end
	end
	
	
# 	puts "If the junction is not annotated then create alternative mRNAs"
	if annotated == false
		overlapping_mrnas.each do |mrna|
			junction_start_exon_idx = mrna.exon_index(junction.start)
			junction_stop_exon_idx = mrna.exon_index(junction.stop)
			
			if junction_start_exon_idx != nil && junction_stop_exon_idx != nil && (junction_stop_exon_idx - junction_start_exon_idx) > 1 && junction.start >= mrna.cds_start && junction.stop <= mrna.cds_stop
				# puts "Creating alternative mRNA that could explain the junction"
				alternative_exon_starts = []
				alternative_exon_stops = []
				for i in 0..mrna.exon_starts.length-1
					if (i <= junction_start_exon_idx) || (i >= junction_stop_exon_idx)
						alternative_exon_starts << mrna.exon_starts[i]
						alternative_exon_stops << mrna.exon_stops[i]
					end
				end
				
				# puts "Check for mrna duplicates, create the alternative mrnas and create the bed file and the fasta file with alternative (unannotated) mrna sequences"
				if alt_mrna_checklist.has_key?(alternative_exon_starts.join(",") + "\t" + alternative_exon_stops.join(","))
					next
				else
					alternative_mrna = Mrna.new(
						mrna.name + "-unann-" + count.to_s,
						mrna.chrom,
						mrna.strand,
						mrna.start,
						mrna.stop,
						mrna.cds_start,
						mrna.cds_stop,
						alternative_exon_starts.length,
						alternative_exon_starts,
						alternative_exon_stops,
						mrna.genename,
						mrna.product,
						mrna.accno + "-unann-" + count.to_s,
						mrna.prot_accno + "-unann-" + count.to_s
					)
					
					alt_mrna_checklist[alternative_exon_starts.join(",") + "\t" + alternative_exon_stops.join(",")] = nil
					# puts ">> Creating bed file for alternative mrnas"
					unannotated_mrnas_bed_output.puts alternative_mrna.cds_to_bed
					# puts ">> Creating fasta file for alternative mrna sequences"
					fasta_entry = alternative_mrna.cds_to_fasta(danRer7_fap)
					unannotated_mrnas_seqs_output.puts fasta_entry
					# puts ">> Creating lineformat for alternative mrna sequences"
					lineformat_entry = alternative_mrna.cds_to_lineformat
					unannotated_mrnas_lineformat_output.puts lineformat_entry
					# puts ">> Creating fasta file with the translated non-annotated mrnas (unannotated proteins)"
					protein = alternative_mrna.cds_to_protein(danRer7_fap)
					unannotated_proteins_output.puts protein
					count = count + 1
				end
			end
		end
	end
end
unannotated_mrnas_seqs_output.close
unannotated_mrnas_lineformat_output.close
unannotated_proteins_output.close
unannotated_mrnas_bed_output.close

