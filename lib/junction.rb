class Junction
	attr_accessor :chrom, :start, :stop, :name, :score, :strand, :graph_start, :graph_stop, :rgb, :exons_num, :exons_lengths, :exons_relative_starts

	def initialize(chrom, start, stop, name, score, strand, graph_start, graph_stop, rgb, exons_num, exons_lengths, exons_relative_starts)
		@chrom = chrom
		@start = start.to_i
		@stop = stop.to_i
		@name = name
		@score = score
		@strand = strand
		@graph_start = graph_start.to_i
		@graph_stop = graph_stop.to_i
		@rgb = rgb
		@exons_num = exons_num.to_i
		@exons_lengths = exons_lengths
		@exons_relative_starts = exons_relative_starts
	end

	# extend the junction in both directions => not used method
	def extend(upstream,downstream)
		if @strand == "+"
			@start = @start - upstream
			@graph_start = @graph_start - upstream
			@stop = @stop + downstream
			@graph_stop = @graph_stop + downstream
			@exons_relative_starts[1] += upstream
			@exons_lengths[0] += upstream
			@exons_lengths[1] += downstream 
		else
			@start = @start - downstream
			@graph_start = @graph_start - downstream
			@stop = @stop + upstream
			@graph_stop = @graph_stop + upstream
			@exons_relative_starts[1] += downstream
			@exons_lengths[0] += downstream
			@exons_lengths[1] += upstream 
		end
	end

	def to_bed()
		bed_entry = [@chrom, @start, @stop, @name, @score, @strand, @graph_start, @graph_stop, @rgb, @exons_num, @exons_lengths.join(","), @exons_relative_starts.join(",")].join("\t")
		return bed_entry
	end
end
