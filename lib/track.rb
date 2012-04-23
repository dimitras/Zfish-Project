require 'mrna'
require 'junction'

class Track # a track is a collection of genomic regions

	def initialize()
		@track_list = {}
		@sorted = false
	end

	def insert_genomic_region(genomic_region)
		key = genomic_region.chrom + "|" + genomic_region.strand
		if !@track_list.has_key?(key) then
			@track_list[key] = []
		end
		@track_list[key] << genomic_region
		@sorted = false
	end
	
	def sort()
		@track_list.each_key do |key|
			@track_list[key].sort! { |a,b| a.start <=> b.start }
		end
		@sorted = true
	end

	def regions_that_contain(genomic_region)
		regions = []
		key = genomic_region.chrom + "|" + genomic_region.strand
		if @track_list.has_key?(key)		
			@track_list[genomic_region.chrom + "|" + genomic_region.strand].each do |region|
				if region.start <= genomic_region.start && genomic_region.stop <= region.stop then
					regions << region
				end
			end
		end
		return regions
	end
end
