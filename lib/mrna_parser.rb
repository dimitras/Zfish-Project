require 'mrna'

class MrnaParser
	attr_accessor :filename, :count

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename) #if filename
		@index = {}
		@index_by_prot_accno = {}
		create_indexes
		@count = 0
	end

	def increment()
		@count += 1
	end
		
	def self.open(filename)
		mrnap = MrnaParser.new(filename)
		if block_given?
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
	
	def create_indexes()
		@filehandle.readline # skip first line which is the header
		temp_pos = @filehandle.pos
		@filehandle.each do |line|
			line_length = line.length
			(name,chrom,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts_str,exon_stops_str,genename,product,accno,prot_accno) = line.chomp!.split("\t")
# 			if @index_by_prot_accno.has_key?(prot_accno.to_s)
# 				puts prot_accno.to_s + "\t" + cds_start.to_s + "\t" + cds_stop.to_s 
# 			end
			@index[accno.to_s] = @filehandle.pos - line_length
			@index_by_prot_accno[prot_accno.to_s] = @filehandle.pos - line_length
			
		end
		@filehandle.pos = temp_pos
	end

	def mrna(accno)
		@filehandle.pos = @index[accno]		
		line = @filehandle.readline.chomp
		return line_parse(line)
	end
	
	def mrna_by_prot_accno(accno)
		pos = @index_by_prot_accno[accno.to_s]
		if pos == nil
			return nil
		else 
			@filehandle.pos = pos
			line = @filehandle.readline.chomp
			return line_parse(line)
		end
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

