class GenomeFastaParser
	attr_accessor :filename

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename)
		@entrypos_by_id = {}
		create_seq_index
	end
	
	def create_seq_index()
		@filehandle.each do |line|
			if line[0..0] == ">"
				header = line[1..line.length-1].chomp
				@entrypos_by_id[header] = @filehandle.pos
			end
		end
		@filehandle.rewind
	end

	def byte_offset(header, genomic_location) 
		pos = @entrypos_by_id[header] + genomic_location + (genomic_location / 50).to_i
		return pos
	end

	def region_seq(header, start, stop)
		seq = ''
		start_pos = byte_offset(header, start)
		@filehandle.pos = start_pos
		region_length = stop - start + 1
		@filehandle.each do |line|
			line.chomp!
			if line.length < region_length
				seq = seq + line
				region_length = region_length - line.length
			else
				seq = seq + line[0..region_length-1]
				return seq
			end
		end
	end

	def self.open(filename)
		fp = GenomeFastaParser.new(filename)
		if block_given?
			fp.each do |entry|
				yield entry
			end
		else
			return fp
		end
	end
	
	def each()
		header,seq = '',''
		@filehandle.each do |line|
			if ((line[0..0] == ">") and (header == ''))
				header = line.chomp
				header.slice!(0)
				seq = ''
			elsif (line[0..0] == ">")
				yield Genome_Entry.new(header, seq)
				header = line.chomp
				header.slice!(0)
				seq = ''
			else
				seq += line.chomp
			end
		end
		yield Genome_Entry.new(header, seq)
	end

end

class Genome_Entry
	attr_accessor :header, :seq

	def initialize(header, seq)
	       @header = header
	       @seq = seq
	end
end
