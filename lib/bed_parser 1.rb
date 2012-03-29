require 'junction'

class BedParser
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
		bedp = BedParser.new(filename)
		if block_given? # call it with a block and not with each()
			bedp.each do |junction|
				yield junction
			end
		else
			return bedp
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
			chrom_start = line.split("\t")[1].chomp
			@index[chrom_start.to_s] = @filehandle.pos - line.length	
		end
		@filehandle.pos = temp_pos
	end
	
	def junction(start)
		@filehandle.pos = @index[start]
		line = @filehandle.readline.chomp
		return line_parse(line)
	end
	
	def line_parse(line)
		(chrom, start, stop, name, score, strand, graph_start, graph_stop, rgb, exons_num, exons_lengths_str, exons_relative_starts_str) = line.split("\t")
		exons_lengths = exons_lengths_str.split(",").map{|pos| pos.to_i}
		exons_relative_starts = exons_relative_starts_str.split(",").map{|pos| pos.to_i}	
		return Junction.new(chrom, start, stop, name, score, strand, graph_start, graph_stop, rgb, exons_num, exons_lengths, exons_relative_starts)
	end
	
end