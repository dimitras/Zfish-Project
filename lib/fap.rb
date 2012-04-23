class FastaParser
	attr_accessor :filename, :count, :entrypos

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename)
		@count = 0
		@entrypos = []
		@entrypos_by_id = {}
		create_index
	end
	
	def increment()
		@count += 1
	end
	
	def first()
		return self.entry(1)
	end

	def last()
		return self.entry(@entrypos.length)
	end

	def entry(entryno)
		@filehandle.pos = @entrypos[entryno-1] # using 1-based numbering scheme
		line = @filehandle.readline.chomp
		line =~ />(.+)/
		header, seq = $1, ''
		@filehandle.each do |line|
			line = line.chomp
			if (line[0..0] == ">")
				break
			end
			seq = seq + line
		end
		return Entry.new(header, seq)
	end
	
	def create_index()
		@filehandle.each do |line|
			if line[0..0] == ">" then
				header = line.chomp.split('.')[0]
				header.slice!(0).to_s
				@entrypos_by_id[header] = @filehandle.pos - line.length
				@entrypos << @filehandle.pos - line.length
				self.increment
			end
		end
		@filehandle.rewind
	end
	
	def entry_by_id(header)
		if @entrypos_by_id.has_key?(header) then
			@filehandle.pos = @entrypos_by_id[header] # @entrypos_by_id is a HASH
			@filehandle.readline
			seq = ''
			@filehandle.each do |line|
				if line[0..0] == ">" then
					break
				end
				seq = seq + line.chomp
			end
			return Entry.new(header, seq)
		else
			return nil
		end
	end
	
	def self.open(filename)
		fp = FastaParser.new(filename)
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
				yield Entry.new(header, seq)
				header = line.chomp
				header.slice!(0)
				seq = ''
			else
				seq += line.chomp
			end
		end
		yield Entry.new(header, seq)
	end
end

class Entry
	attr_accessor :header, :seq

	def initialize(header, seq)
	       @header = header
	       @seq = seq
	end

	def length()
		return seq.length
	end
end
