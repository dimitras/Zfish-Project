require 'peptide'
require "test/unit"
 
class TestPeptide < Test::Unit::TestCase
 	
	def setup
		@peptide = Peptide.new("NP_001032499", "DTEVR", 42, 46, 1)
		@mrna = Mrna.new("NM_001037422", "chr7", "-", 30447188, 30452289, 30447714, 30452101, 2, [30447188,30451964], [30448681,30452289], "fbxl8", "F-box/LRR-repeat protein 8", "NM_001037422", "NP_001032499")
		@peptide.mrna = @mrna
	end
	
	def test_peptides
		assert_kind_of(Peptide, @peptide)
	end

	def test_mrnas
		assert_kind_of(Mrna, @mrna)
	end
	
	def test_peptide
		assert_equal("NP_001032499", @peptide.prot_acc)
		@peptide.genomic_locations
		assert_equal(30448681, @peptide.genomic_starts.first)
	end
	
	def test_to_bed
		assert_equal("chr7\t30448681\t30451979\tDTEVR\t1\t-\t30448681\t30451979\t255,0,0\t2\t1,15\t0,3283", @peptide.to_bed)
	end

end
#OK