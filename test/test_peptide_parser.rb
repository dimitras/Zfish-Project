require 'peptide_parser'
require "test/unit"
 
class TestPeptideParser < Test::Unit::TestCase
 	
	def setup
		@peptidep = PeptideParser.open("../data/zfish_data/mascpep.txt")
	end
	
	def test_peptidefile
		assert_kind_of(PeptideParser, @peptidep)
	end

	def test_peptide
		assert_equal("NP_571501", @peptidep.peptide("NP_571501").prot_acc)
	end
	
	def test_line_parse
		line = "NP_001032499	DTEVR	42	46	1"
		assert_equal("NP_001032499", @peptidep.line_parse(line).prot_acc)
	end

end
#OK