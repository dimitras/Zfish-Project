require 'junction'
require "test/unit"
 
class TestJunction < Test::Unit::TestCase
 	
	def setup
		@junction = Junction.new("chr1", 4759, 5227, 2, 2, "+", 4759, 5227, "16,78,139", 2, [50,50], [0,418])
	end
	
	def test_junctions
		assert_kind_of(Junction, @junction)
	end

	def test_junction
		assert_equal(4759, @junction.start)
	end
	
	def test_to_bed
		assert_equal("chr1\t4759\t5227\t2\t2\t+\t4759\t5227\t16,78,139\t2\t50,50\t0,418", @junction.to_bed)
	end

end
#OK