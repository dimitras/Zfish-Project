class Locus
	attr_accessor :chrom, :strand, :start, :stop

	def initialize(chrom, strand, start, stop)
		@chrom = chrom
		@strand = strand
		@start = start
		@stop = stop
	end
end
