require 'peptide'
# require 'mrna'

class PeptideParser
	attr_accessor :filename

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename)
		@index = {}
		create_index
	end
	
	def self.open(filename)
		peptidep = PeptideParser.new(filename)
		if block_given?
			peptidep.each do |peptide|
				yield peptide
			end
		else
			return peptidep
		end
	end

	def each()
		@filehandle.each do |line|
			yield line_parse(line.chomp)
		end
	end

	def line_parse(line)
		(prot_acc, peptide, start, stop, multiplicity) = line.split("\t")
		return Peptide.new(prot_acc.to_s, peptide, start, stop, multiplicity)
	end

	def create_index()
		temp_pos = @filehandle.pos
		@filehandle.each do |line|
			prot_acc = line.split("\t")[0].chomp
			@index[prot_acc.to_s] = @filehandle.pos - line.length
		end
		@filehandle.pos = temp_pos
	end
	
	def peptide(accno)
		@filehandle.pos = @index[accno]
		line = @filehandle.readline.chomp
		return line_parse(line)
	end

end


# mrnap = MrnaParser.open(ARGV[0])
# peptidep = PeptideParser.open(ARGV[1])
# bep = File.open(ARGV[2],"w")

# peptidep.each do |peptide|
#  	peptide.mrna = mrnap.mrna(peptide.prot_acc)
#  	bep.puts peptide.to_bed()
# end
