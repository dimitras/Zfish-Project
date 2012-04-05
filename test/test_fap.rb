require 'fap'
require "test/unit"
 
class TestFastaParser < Test::Unit::TestCase
 	
	def setup
		@fap = FastaParser.open("../data/UCSC_mrnaseqs.fa")
		counter = 0
	   	@fap.each do |entry|
	   		if counter == 1
		   		@entry = entry
		   	end
	   		counter +=1
	   	end
	end

	def test_fapfile
		assert_kind_of(FastaParser, @fap)
	end

	def test_entries_number
		assert_equal(14674, @fap.count())
	end

	def test_entry
		assert_equal("NM_131426", @fap.first.header)
		assert_equal("NM_001128401", @fap.last.header)
		assert_equal("NM_001110522", @fap.entry(2).header)
		assert_equal("NM_001110522", @fap.entry_by_id("NM_001110522").header)
		assert_equal(2819, @fap.first.seq.length)
	end	

	def test_each
		assert_equal("NM_001110522", @entry.header)
	end
	
end
#OK