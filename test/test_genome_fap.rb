require "genome_fap"
require "mrna"
require "test/unit"
 
class TestGenomeFastaParser < Test::Unit::TestCase
 	
	def setup
	   @fp = GenomeFastaParser.open("../data/zfish_data/danRer7.fa")
	   @mrna = Mrna.new("NM_001109719", "chr1", "+", 2469741, 2673265, 2469990, 2673265, 4, [2469741,2546653,2642202,2673100], [2470146,2546811,2642593,2673265], "zgc:171629", "glypican-6 precursor", "NM_001109719", "NP_001103189")
	end

	def test_genomefile()
		assert_kind_of(GenomeFastaParser, @fp)
	end

	def test_mrnas()
		assert_kind_of(Mrna, @mrna)
	end

	def test_byte_offset() 
		assert_equal(59054243, @fp.byte_offset(@mrna.chrom, @mrna.cds_start))
	end	

	def test_region_seq()
		assert_equal('ATGGTGAAGACACCTGTCGTGGTTTTTACTT', @fp.region_seq(@mrna.chrom, @mrna.cds_start, @mrna.cds_start+30))
	end
 
end