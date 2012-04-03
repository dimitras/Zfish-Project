# USAGE: ruby runner.rb ../data/zfish_data/ucsc_mrnas.txt ../data/zfish_data/mascpep.txt ../data/zfish_data/UCSC_protein_refseqs.fa ../data/zfish_output/peptides.bed ../data/zfish_data/junctions_high-quality.bed ../data/zfish_output/unannotated_mrnas.bed ../data/zfish_data/danRer7.fa ../data/zfish_output/unannotated_mrnas_seqs.fa ../data/zfish_output/unannotated_proteins.fa > ../data/zfish_output/puts.log

# USAGE IF: have_unannotated = true : ruby runner.rb ../data/zfish_data/ucsc_mrnas.txt ../data/zfish_output/mascot_unannotated_peptides.txt ../data/zfish_data/UCSC_protein_refseqs.fa ../data/zfish_output/unannotated_peptides.bed > ../data/zfish_output/puts.log

#!/usr/bin/env ruby
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

have_unannotated = true

puts "Initializing (Getting command arguments)"

ucsc_mrnas_refseq = ARGV[0] # ../data/zfish_data/ucsc_mrnas.txt
mascot_peptides = ARGV[1] # ../data/zfish_data/mascpep.txt
ucsc_protein_refseq = ARGV[2] # ../data/zfish_data/UCSC_protein_refseqs.fa
peptides_bed = ARGV[3] # peptides in bed format for ucsc # ../data/zfish_output/peptides.bed

ucsc_mrnas_refseq_mrnap = MrnaParser.open(ucsc_mrnas_refseq)
mascot_peptides_pep = PeptideParser.open(mascot_peptides)
ucsc_protein_refseq_fap = FastaParser.open(ucsc_protein_refseq) 
peptides_bedp = File.open(peptides_bed,"w")

if have_unannotated == false
	rum_junctions_hq_bed = ARGV[4] # ../data/zfish_data/junctions_high-quality.bed
	unannotated_mrnas_bed = ARGV[5] # ../data/zfish_output/unannotated_mrnas.bed
	danRer7_refseq = ARGV[6] # ../data/zfish_data/danRer7.fa
	unannotated_mrnas_seqs = ARGV[7] # ../data/zfish_output/unannotated_mrnas_seqs.fa
	unannotated_proteins = ARGV[8] # ../data/zfish_output/unannotated_proteins.fa

	rum_junctions_hq_peptides_bedp = BedParser.open(rum_junctions_hq_bed)
	unannotated_mrnas_bed_outp = File.open(unannotated_mrnas_bed,"w")
	danRer7_refseq_fap = GenomeFastaParser.new(danRer7_refseq)
	unannotated_mrnas_seqs_outp = File.open(unannotated_mrnas_seqs, "w")
	unannotated_mrnas_seqs_fap = GenomeFastaParser.open(unannotated_mrnas_seqs)
	unannotated_proteins_outp = File.open(unannotated_proteins,"w")
end


puts "Create bed file from peptides"
peptides_bedp.puts "track name=\"Peptides\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\""
mascot_peptides_pep.each do |peptide|
 	peptide.mrna = ucsc_mrnas_refseq_mrnap.mrna(peptide.prot_acc)
	protein = ucsc_protein_refseq_fap.entry_by_id(peptide.prot_acc)
	if peptide.mrna.cds_length != (protein.length + 1) * 3 then # +1 because mrna stop codon is not translated to protein
		next
	end
 	peptides_bedp.puts peptide.to_bed
end
peptides_bedp.close


if have_unannotated == false

	puts "Create bed file for non-annotated junctions"
	mrna_track = Track.new
	ucsc_mrnas_refseq_mrnap.each do |mrna|
		mrna_track.insert_genomic_region(mrna)
	end


	puts "Scanning junctions to find unannotated ones and create alternative gene models that would explain them" 
	alt_mrna_checklist = {}
	count = 0
	cds_seq = ''
	unannotated_mrnas_bed_outp.puts "track name=\"Alternative gene models\" description=\"Alternative gene models\" visibility=2 itemRgb=\"On\""
	rum_junctions_hq_peptides_bedp.each do |junction|
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
					
					puts "Check for mrna duplicates, create the alternative mrnas and create the bed file and the fasta file with alternative (unannotated) mrna sequences"
					if alt_mrna_checklist.has_key?(alternative_exon_starts.join(",") + "\t" + alternative_exon_stops.join(","))
						next
					else
						alternative_mrna = Mrna.new(
							mrna.name + "." + count.to_s,
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
							mrna.accno,
							mrna.prot_accno
						)
						
						alt_mrna_checklist[alternative_exon_starts.join(",") + "\t" + alternative_exon_stops.join(",")] = nil
						# puts ">> Creating bed file for alternative mrnas"
						unannotated_mrnas_bed_outp.puts alternative_mrna.cds_to_bed
						# puts ">> Creating bed file for alternative mrna sequences"
						cds_seq = alternative_mrna.cds_to_fasta(danRer7_refseq_fap)
						unannotated_mrnas_seqs_outp.puts cds_seq
						count = count + 1
					end
				end
			end
		end
	end
	unannotated_mrnas_seqs_outp.close
	unannotated_mrnas_bed_outp.close


	puts "Create fasta file with the translated non-annotated mrnas"
	unannotated_mrnas_seqs_fap.each do |entry|
		mrna = Bio::Sequence::NA.new(entry.seq)
		protein = mrna.translate
		unannotated_proteins_outp.puts protein.to_fasta(entry.header.to_s, 60)
	end
	unannotated_proteins_outp.close

end