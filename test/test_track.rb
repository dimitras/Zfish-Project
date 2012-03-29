require 'track'
require "test/unit"
 
class TestTrack < Test::Unit::TestCase
 	
	def setup
		@mrna_track = Track.new
		@junction_track = Track.new
		@junction = Junction.new("chr1", 4759, 5227, 2, 2, "+", 4759, 5227, "16,78,139", 2, [50,50], [0,418])
	end

	def test_tracks
		assert_kind_of(Track, @mrna_track)
	end
	
	def test_regions_that_contain_genomic_region
		regions = ["chr1", 4759, 5227, 2, 2, "+", 4759, 5227, "16,78,139", 2, [50,50], [0,418]]
		assert_equal(regions, @mrna_track.regions_that_contain(@junction))		
	end

end
#OK