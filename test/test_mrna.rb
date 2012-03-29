require 'mrna'
require 'fap'
require "test/unit"
 
class TestMrna < Test::Unit::TestCase
 	
	def setup
		@mrna = Mrna.new("NM_001109719", "chr1", "+", 2469741, 2673265, 2469990, 2673265, 4, [2469741,2546653,2642202,2673100], [2470146,2546811,2642593,2673265], "zgc:171629", "glypican-6 precursor", "NM_001109719", "NP_001103189")
		@danRer7_refseq_fap = FastaParser.new("../data/zfish_data/danRer7.fa")
	end
	
	def test_mrnas
		assert_kind_of(Mrna, @mrna)
	end

	def test_fap
		assert_kind_of(FastaParser, @danRer7_refseq_fap)
	end

	def test_mrna
		assert_equal("NP_001103189", @mrna.prot_accno)
	end
	
	def test_cds_length
		assert_equal(874, @mrna.cds_length)
	end

	def test_cds_to_bed
		assert_equal("chr1\t2469990\t2673266\tNM_001109719\t1\t+\t2469990\t2673266\t0,255,0\t4\t157,159,392,166\t0,76663,172212,203110", @mrna.cds_to_bed)
	end

	def test_cds_to_fasta
		assert_equal("acacacacacaca", @mrna.cds_to_fasta(@danRer7_refseq_fap))
	end

end
#OK