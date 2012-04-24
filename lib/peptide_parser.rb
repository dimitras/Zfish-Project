require 'rubygems'
require 'mascot/dat'
require 'mascot/dat/peptides'
require 'mascot/dat/psm'
require 'peptide'

class PeptideParser
	def initialize(filename)
		@dat = Mascot::DAT.open(filename, true)
		@peptides = @dat.peptides
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
		@peptides.each do |psm|
			psm.proteins.each do |protein_info|
				protein_accno = protein_info[0].split(".")[0]
				start = protein_info[2]
				stop = protein_info[3]
				multiplicity = protein_info[4]
				yield Peptide.new(protein_accno, psm.pep, start-1, stop-1, multiplicity, psm.score)
			end
		end
	end

	def create_index()
		@peptides.each do |psm|
			psm.proteins.each do |protein_info|
				protein_accno = protein_info[0].split(".")[0]
				start = protein_info[2]
				stop = protein_info[3]
				multiplicity = protein_info[4]
				peptide = Peptide.new(protein_accno, psm.pep, start-1, stop-1, multiplicity, psm.score)
				@index[protein_accno] = peptide
			end
		end
		rewind
	end

	def peptide(prot_acc)
		peptide = @index[prot_acc]
		return peptide
	end

	def rewind()
		@peptides.rewind
	end
end
