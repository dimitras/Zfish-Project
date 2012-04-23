#!/usr/bin/env ruby

# USAGE: ruby scripts/clean_files.rb data/mrnas_ucsc_lineformat.txt data/proteins_ucsc_db.fa data/danRer7.fa results/mrnas_ucsc_lineformat_clean.txt results/proteins_ucsc_db_clean.fa

# clean mrnas ucsc lineformat file of mrnas that their length is not equal to the relative protein length

$LOAD_PATH << './lib'

require 'mrna'
require 'mrna_parser'
require 'fap'
require 'genome_fap'


puts "Initializing (Getting command arguments)"
mrnas_lineformat_file = ARGV[0]
proteins_fasta_file = ARGV[1]
danRer7_fasta_file = ARGV[2]
mrnas_clean_lineformat_file = ARGV[3]
proteins_clean_fasta_file = ARGV[4]

# Input
mrnas_lineformat_mrnap = MrnaParser.open(mrnas_lineformat_file)
proteins_fap = FastaParser.open(proteins_fasta_file) 
danRer7_gfap = GenomeFastaParser.new(danRer7_fasta_file)

# Output
mrnas_clean_lineformat_output = File.open(mrnas_clean_lineformat_file,"w")
proteins_clean_fasta_output = File.open(proteins_clean_fasta_file,"w")


count = 0
mrnas_clean_lineformat_output.puts "# danRer7.unannotated.name\tdanRer7.unannotated.chrom\tdanRer7.unannotated.strand\tdanRer7.unannotated.txStart\tdanRer7.unannotated.txEnd\tdanRer7.unannotated.cdsStart\tdanRer7.unannotated.cdsEnd\tdanRer7.unannotated.exonCount\tdanRer7.unannotated.exonStarts\tdanRer7.unannotated.exonEnds\tdanRer7.unannotated.name2\tdanRer7.unannotated.product\tdanRer7.unannotated.mrnaAcc\tdanRer7.unannotated.protAcc"
mrnas_lineformat_mrnap.each do |mrna|
	if mrna.prot_accno != nil && mrna.prot_accno != ''
		protein = proteins_fap.entry_by_id(mrna.prot_accno)
		if (((protein.length + 1) * 3) == mrna.cds_length)
			new_mrna = Mrna.new(
				mrna.name + "-ann-" + count.to_s,
				mrna.chrom,
				mrna.strand,
				mrna.start,
				mrna.stop,
				mrna.cds_start,
				mrna.cds_stop,
				mrna.exon_count,
				mrna.exon_starts,
				mrna.exon_stops,
				mrna.genename,
				mrna.product,
				mrna.accno + "-ann-" + count.to_s,
				mrna.name + "-ann-" + count.to_s
			)
			new_protein = new_mrna.cds_to_protein(danRer7_gfap)
			
			proteins_clean_fasta_output.puts new_protein
			mrnas_clean_lineformat_output.puts new_mrna.to_lineformat
			
			count = count + 1
		end
	end
end
mrnas_clean_lineformat_output.close
proteins_clean_fasta_output.close
