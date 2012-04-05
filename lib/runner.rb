# USAGE: ruby runner.rb ../data/ucsc_mrnas.txt ../data/mascpep.txt ../data/UCSC_protein_refseqs.fa ../out/peptides.bed ../data/junctions_high-quality.bed ../out/unannotated_mrnas.bed ../data/danRer7.fa ../out/unannotated_mrnas_seqs.fa ../out/unannotated_mrnas_refseqs.txt ../out/unannotated_proteins.fa > ../out/puts.log

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


puts "Initializing (Getting command arguments)"

ucsc_mrnas_refseq = ARGV[0] # ../data/ucsc_mrnas.txt
mascot_peptides = ARGV[1] # ../data/mascpep.txt
ucsc_protein_refseq = ARGV[2] # ../data/UCSC_protein_refseqs.fa
peptides_bed = ARGV[3] # peptides in bed format for ucsc # ../out/peptides.bed
rum_junctions_hq_bed = ARGV[4] # ../data/junctions_high-quality.bed
unannotated_mrnas_bed = ARGV[5] # ../out/unannotated_mrnas.bed
danRer7_refseq = ARGV[6] # ../data/danRer7.fa
unannotated_mrnas_seqs = ARGV[7] # ../out/unannotated_mrnas_seqs.fa
unannotated_mrnas_refseqs = ARGV[8] # ../out/unannotated_mrnas_refseqs.txt
unannotated_proteins = ARGV[9] # ../out/unannotated_proteins.fa

ucsc_mrnas_refseq_mrnap = MrnaParser.open(ucsc_mrnas_refseq)
mascot_peptides_pep = PeptideParser.open(mascot_peptides)
ucsc_protein_refseq_fap = FastaParser.open(ucsc_protein_refseq) 
peptides_bedp = File.open(peptides_bed,"w")
rum_junctions_hq_peptides_bedp = BedParser.open(rum_junctions_hq_bed)
unannotated_mrnas_bed_outp = File.open(unannotated_mrnas_bed,"w")
danRer7_refseq_fap = GenomeFastaParser.new(danRer7_refseq)
unannotated_mrnas_seqs_outp = File.open(unannotated_mrnas_seqs, "w")
unannotated_mrnas_seqs_fap = GenomeFastaParser.open(unannotated_mrnas_seqs)
unannotated_mrnas_refseqs_outp = File.open(unannotated_mrnas_refseqs,"w")
unannotated_proteins_outp = File.open(unannotated_proteins,"w")


puts "Create bed file from peptides"
peptides_bedp.puts "track name=\"Peptides\" description=\"Peptides\" visibility=2 itemRgb=\"On\""
mascot_peptides_pep.each do |peptide|
 	peptide.mrna = ucsc_mrnas_refseq_mrnap.mrna(peptide.prot_acc)
	protein = ucsc_protein_refseq_fap.entry_by_id(peptide.prot_acc)
	if peptide.mrna.cds_length != (protein.length + 1) * 3 then # +1 because mrna stop codon is not translated to protein
		next
	end
 	peptides_bedp.puts peptide.to_bed
end
peptides_bedp.close


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
unannotated_mrnas_refseqs_outp.puts "# danRer7.unannotated.name\tdanRer7.unannotated.chrom\tdanRer7.unannotated.strand\tdanRer7.unannotated.txStart\tdanRer7.unannotated.txEnd\tdanRer7.unannotated.cdsStart\tdanRer7.unannotated.cdsEnd\tdanRer7.unannotated.exonCount\tdanRer7.unannotated.exonStarts\tdanRer7.unannotated.exonEnds\tdanRer7.unannotated.name2\tdanRer7.unannotated.product\tdanRer7.unannotated.mrnaAcc\tdanRer7.unannotated.protAcc"

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
						mrna.name + "-" + count.to_s,
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
					# puts ">> Creating fasta file for alternative mrna sequences"
					cds_seq = alternative_mrna.cds_to_fasta(danRer7_refseq_fap)
					unannotated_mrnas_seqs_outp.puts cds_seq
					# puts ">> Creating refseq for alternative mrna sequences"
					mrnarefseq_entry = alternative_mrna.cds_to_mrnarefseq
					unannotated_mrnas_refseqs_outp.puts mrnarefseq_entry
					# puts ">> Creating fasta file with the translated non-annotated mrnas (unannotated proteins)"
					protein = alternative_mrna.cds_to_protein(danRer7_refseq_fap)
					unannotated_proteins_outp.puts protein
					count = count + 1
				end
			end
		end
	end
end
unannotated_mrnas_seqs_outp.close
unannotated_mrnas_refseqs_outp.close
unannotated_proteins_outp.close
unannotated_mrnas_bed_outp.close


# puts "Create fasta file with the translated non-annotated mrnas"
# unannotated_mrnas_seqs_fap.each do |entry|
# 	mrna = Bio::Sequence::NA.new(entry.seq)
# 	protein = mrna.translate
# 	unannotated_proteins_outp.puts protein.to_fasta(entry.header.to_s, 60)
# end
# unannotated_proteins_outp.close
