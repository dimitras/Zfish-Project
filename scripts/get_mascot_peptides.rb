# USAGE : ruby scripts/mascot_peptides.rb data/1d-zf-16de.dat > results/mascot_peptides.txt
# USAGE : ruby scripts/mascot_peptides.rb data/1d-zf-16de-full-db.dat > results/mascot_peptides_full_db.txt

#!/usr/bin/env ruby

$LOAD_PATH << './lib'

require 'rubygems'
require 'mascot/dat'
require 'mascot/dat/peptides'
require 'mascot/dat/psm'


@dat = Mascot::DAT.open(ARGV[0], true)
@peptides = @dat.peptides
@peptides.each do |psm|
	psm.proteins.each do |protein_info|
		protein_name = protein_info[0]
		start_pos = protein_info[2] - 1
		stop_pos = protein_info[3] - 1
		multiplicity = protein_info[4]
		puts protein_name + "\t" + psm.pep + "\t" + start_pos.to_s + "\t" + stop_pos.to_s + "\t" + multiplicity.to_s + "\t" + psm.score.to_s
	end
end
