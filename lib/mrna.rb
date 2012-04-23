require 'rubygems'
require 'bio'

class Mrna
	attr_accessor :name, :chrom, :strand, :start, :stop, :cds_start, :cds_stop, :exon_count, :exon_starts, :exon_stops, :genename, :product, :accno, :prot_accno, :coding

	def initialize(name,chrom,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts,exon_stops,genename,product,accno,prot_accno)
		@name = name
		@chrom = chrom
		@strand = strand
		@start = start.to_i
		@stop = stop.to_i
		@cds_start = cds_start.to_i
		@cds_stop = cds_stop.to_i
		@exon_count = exon_count.to_i
		@exon_starts = exon_starts
		@exon_stops = exon_stops
		@genename = genename
		@product = product
		@accno = accno
		@prot_accno = prot_accno
		@cds_exon_starts = []
		@cds_exon_stops = []
		if (@cds_stop > @cds_start)
			@coding = true
		else
			@coding = false
		end
	end
	
	def cds_seq_from_genome(genome_fap)
		cds_seq = Bio::Sequence::NA.new("")
		for i in 0..self.cds_exon_starts.length-1
			cds_seq = cds_seq + genome_fap.region_seq(@chrom, self.cds_exon_starts[i], self.cds_exon_stops[i])
			if @strand == "-"
				cds_seq = cds_seq.complement
			end
		end
		return cds_seq
	end
	
	def cds_to_fasta(genome_fap)
		cds_seq = cds_seq_from_genome(genome_fap)
		return ">" + @name.to_s + "\n" + cds_seq
	end

	def cds_to_protein(genome_fap)
		mrna_seq = Bio::Sequence::NA.new(cds_seq_from_genome(genome_fap))
		protein = mrna_seq.translate
		return protein.to_fasta(@name, 50)
	end

	def cds_to_lineformat()
		exon_stops_incremented = []
		for i in 0..@cds_exon_stops.length-1
			exon_stops_incremented[i] = @cds_exon_stops[i] + 1
		end
		entry = [@name, @chrom, @strand, @cds_start, @cds_stop + 1, @cds_start, @cds_stop + 1, @cds_exon_starts.length, @cds_exon_starts.join(","), exon_stops_incremented.join(","), @genename, @product, @accno, @name].join("\t")
		return entry
	end

	def to_lineformat()
		exon_stops_incremented = []
		for i in 0..@exon_stops.length-1
			exon_stops_incremented[i] = @exon_stops[i] + 1
		end
		entry = [@name, @chrom, @strand, @start, @stop + 1, @cds_start, @cds_stop + 1, @exon_starts.length, @exon_starts.join(","), exon_stops_incremented.join(","), @genename, @product, @accno, @prot_accno].join("\t")
		return entry
	end
	
	def cds_to_bed()
		rgb = "0,255,0"
		if @strand == "-" then
			rgb = "255,0,0"
		end
		
		cds_exons_lengths = []
		relative_starts = []
		for i in 0..self.cds_exon_starts.length-1
			cds_exons_lengths << (self.cds_exon_stops[i] - self.cds_exon_starts[i] + 1)
			relative_starts << self.cds_exon_starts[i] - self.cds_exon_starts.first
		end

		bed_entry = [@chrom, @cds_start, @cds_stop+1, @name, 1, @strand, @cds_start, @cds_stop+1, rgb, @cds_exon_starts.length, cds_exons_lengths.join(","), relative_starts.join(",")].join("\t")
		return bed_entry
	end

	def exon_index(genomic_position) # from junction genomic region
		for i in 0..@exon_starts.length-1
			if @exon_starts[i] <= genomic_position && genomic_position <= @exon_stops[i]
				return i # return the index of the exons array
			end
		end
		return nil
	end

	def cds_length()
		if @coding
			cdseq_length = 0
			for i in 0..self.cds_exon_starts.length-1
				cdseq_length += self.cds_exon_stops[i] - self.cds_exon_starts[i] + 1
			end
			return cdseq_length
		else
			return nil
		end
	end
	
	def create_cds_exons()
		if @cds_start != nil && @cds_stop != nil && @exon_starts.length > 0 && @exon_stops.length > 0
			@cds_exon_starts, @cds_exon_stops = cut_regions_within(@exon_starts,@exon_stops,@cds_start,@cds_stop)
		end
	end
	
	def cds_exon_starts()
		if @cds_exon_starts.length == 0
			create_cds_exons
		end
		return @cds_exon_starts
	end
	
	def cds_exon_stops()
		if @cds_exon_stops.length == 0
			create_cds_exons
		end
		return @cds_exon_stops
	end
	
	def cut_regions_within(starts, stops, start, stop) # cuts the two arrays (starts,stops) to create the parts between the two positions (start,stop)
		in_region = false
		region_starts = []
		region_stops = []
		for i in 0..starts.length-1
			if starts[i].to_i <= start.to_i && start.to_i <= stops[i].to_i
				in_region = true
				region_starts << start
				if starts[i].to_i <= stop.to_i && stop.to_i <= stops[i].to_i
					region_stops << stop
				else
					region_stops << stops[i]
				end
			elsif in_region == true
				if starts[i].to_i <= stop.to_i && stop.to_i <= stops[i].to_i 
					in_region = false
					region_starts << starts[i]
					region_stops << stop
				else
					region_starts << starts[i]
					region_stops << stops[i]
				end
			end
		end
		return region_starts, region_stops
	end
end
