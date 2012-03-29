require 'rubygems'
require "mascot/dat"
require 'mascot/dat/peptides'
require 'mascot/dat/psm'

class MascotPeptides

	@dat = Mascot::DAT.open(ARGV[0], true)
	@peptides = @dat.peptides
	@peptides.each do |psm|
		psm.proteins.each do |protein_info|
			protein_name = protein_info[0].split(".")[0]
			start_pos = protein_info[2]
			stop_pos = protein_info[3]
			multiplicity = protein_info[4]
			puts protein_name + "\t" + psm.pep + "\t" + start_pos.to_s + "\t" + stop_pos.to_s + "\t" + multiplicity.to_s
		end
	end

end
