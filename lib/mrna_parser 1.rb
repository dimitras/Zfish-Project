require 'mrna'

class MrnaParser
	attr_accessor :filename, :count

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename) #if filename
		@index = {}
		create_index
		@count = 0
	end

	def increment()
		@count += 1
	end
		
	def self.open(filename)
		mrnap = MrnaParser.new(filename)
		if block_given? # call it with a block and not with each()
			mrnap.each do |mrna|
				yield mrna
			end
		else
			return mrnap
		end
	end

	def each()
		@filehandle.each do |line|
			self.increment
			yield line_parse(line.chomp)
		end
	end
	
	def create_index()
		@filehandle.readline # skip first line which is the header
		temp_pos = @filehandle.pos
		@filehandle.each do |line|
			prot_accno = line.split("\t")[13].chomp
			@index[prot_accno.to_s] = @filehandle.pos - line.length
		end
		@filehandle.pos = temp_pos
	end
	
	def mrna(accno)
		@filehandle.pos = @index[accno]		
		line = @filehandle.readline.chomp
		return line_parse(line)
	end
	
	def line_parse(line)
		(name,chrom,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts_str,exon_stops_str,genename,product,accno,prot_accno) = line.split("\t")
		# deduct 1 because ucsc database is 1-based, meaning that ucsc gives the stop positions increased by 1
		stop = stop.to_i - 1
		cds_stop = cds_stop.to_i - 1
		exon_starts = exon_starts_str.split(",").map{|pos| pos.to_i}
		exon_stops = exon_stops_str.split(",").map{|pos| pos.to_i - 1}
		return Mrna.new(name,chrom,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts,exon_stops,genename,product,accno,prot_accno)
	end
	
end

# 
# mrnap = MrnaParser.open(ARGV[0])
# mrnap.each do |mrna|
# 	mrna.create_cds_exons
# 	mrna.cds_length
# 	puts "name:" + "\t" + mrna.name + "\t" + "protein:" + "\t" + mrna.prot_accno.to_s + "\t" + mrna.cds_length.to_s
# end
# puts "===> NP_991214"
# puts "cds length: " + mrnap.mrna("NP_991214").cds_length.to_s
# puts "exons\t" + mrnap.mrna("NP_991214").exon_starts.join(" ").to_s + "\t" + mrnap.mrna("NP_991214").exon_stops.join(" ").to_s
# puts "cds exons\t" + mrnap.mrna("NP_991214").cds_exon_starts.join(" ").to_s + "\tAND\t" + mrnap.mrna("NP_991214").cds_exon_stops.join(" ").to_s
# puts "<=== NP_991214"
