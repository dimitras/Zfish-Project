require 'bed_parser'
require "test/unit"
 
class TestBedParser < Test::Unit::TestCase
 	
	def setup
		@junction_bedp = BedParser.open("../data/zfish_data/junctions_high-quality.bed")
	end
	
	def test_bedfile
		assert_kind_of(BedParser, @junction_bedp)
	end

	def test_junction
		assert_equal("chr1", @junction_bedp.junction("4759").chrom)
		assert_equal("+", @junction_bedp.junction("4759").strand)
	end

	def test_to_bed
		assert_equal("chr19	24975573	24975690	1	1	+	24975573	24975690	0,255,127	2	50,50	0,67", @junction_bedp.junction("24975573").to_bed)
	end

end
# OK