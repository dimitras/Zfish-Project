$LOAD_PATH << './lib'

require 'mrna_parser'
require 'test/unit'
 
class TestMrnaParser < Test::Unit::TestCase
 	
	def setup
		@mrnap = MrnaParser.open("data/mrnas_ucsc_lineformat.txt")
		@cleanedmrnas = MrnaParser.open("results/mrnas_ucsc_lineformat_clean.txt")
		@mrnap.each do |entry|  end
	end

	def test_mrnafile
		assert_kind_of(MrnaParser, @mrnap)
		assert_kind_of(MrnaParser, @cleanedmrnas)
	end
	
	def test_mrna_by_mrna_acc
		assert_equal("chr1", @mrnap.mrna("NM_131426").chrom)
		assert_equal("+", @mrnap.mrna("NM_131426").strand)
		assert_equal("50322024 50323684 50327722 50376641 50384688 50384995 50387281 50388021 50392530 50393547", @mrnap.mrna("NM_131426").cds_exon_starts.join(" "))
		assert_equal("50322230 50323750 50327849 50376773 50384781 50385108 50387443 50388128 50392578 50393581", @mrnap.mrna("NM_131426").cds_exon_stops.join(" "))
	end

	def test_mrna_by_protein_acc
		assert_equal("chr1", @mrnap.mrna_by_prot_accno("NP_571501").chrom)
		assert_equal("+", @mrnap.mrna_by_prot_accno("NP_571501").strand)
	end

	def test_count
		assert_equal(15561, @mrnap.count())
	end

	def test_line_parse
		line = "NM_131426	chr1	+	50321633	50410568	50322024	50393582	11	50321633,50323684,50327722,50376641,50384688,50384995,50387281,50388021,50392530,50393547,50409289,	50322231,50323751,50327850,50376774,50384782,50385109,50387444,50388129,50392579,50393588,50410568,	lef1	lymphoid enhancer-binding factor 1	NM_131426	NP_571501"
		assert_equal("NM_131426", @mrnap.line_parse(line).accno)
	end	
end
#OK