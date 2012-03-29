require 'mrna'

class Peptide
	attr_accessor :prot_acc, :peptide, :start, :stop, :multiplicity, :mrna, :genomic_starts, :genomic_stops

	def initialize(prot_acc, peptide, start, stop, multiplicity)
		@prot_acc = prot_acc.to_s
		@peptide = peptide.to_s
		@start = start.to_i
		@stop = stop.to_i
		@multiplicity = multiplicity.to_i
		@mrna = nil
		@genomic_start = nil
		@genomic_stop = nil
		@genomic_starts = []
		@genomic_stops = []
	end

	def length()
		return @stop - @start + 1
	end

	def genomic_locations()
		if self.genomic_starts.length == 0 && self.genomic_stops.length == 0 then
			find_genomic_locations()
		end
		# puts "genomic_locations\t" + self.genomic_starts.to_s + "\t" + self.genomic_stops.to_s + "\n"
		return self.genomic_starts, self.genomic_stops
	end

	def find_genomic_locations()
		if @mrna != nil then
			(@genomic_starts, @genomic_stops) = cut_regions_within(self.mrna.cds_exon_starts,self.mrna.cds_exon_stops,self.genomic_start,self.genomic_stop)
		end
	end
	
	def cut_regions_within(starts, stops, start, stop) # cuts the two arrays (starts,stops) to create the parts between the two positions (start,stop)
		# puts
		# puts "======> cut regions within start\t" + start.to_s + "\tand stop\t" + stop.to_s + "\twithin starts\t" + starts.join(" ").to_s + "\tand stops\t" + stops.join(" ").to_s
		in_region = false
		region_starts = []
		region_stops = []
		for i in 0..starts.length-1
			if starts[i].to_i <= start.to_i && start.to_i <= stops[i].to_i then
				in_region = true
				region_starts << start
				if starts[i].to_i <= stop.to_i && stop.to_i <= stops[i].to_i then
					region_stops << stop
					# puts "region:\t" + region_starts.to_s + "\t" + region_stops.to_s + "\n"
					return region_starts, region_stops
				else
					region_stops << stops[i]
				end
			elsif in_region == true then
				if starts[i].to_i <= stop.to_i && stop.to_i <= stops[i].to_i then
					in_region = false
					region_starts << starts[i]
					region_stops << stop
					# puts "region:\t" + region_starts.to_s + "\t" + region_stops.to_s + "\n"
					return region_starts, region_stops
				else
					region_starts << starts[i]
					region_stops << stops[i]
				end
			end
		end
	end
	
	def genomic_start()
		if @genomic_start == nil then
			# puts "strand: " + self.mrna.strand
			if self.mrna.strand == "+" then
				@genomic_start = to_genomic_location(@start * 3 - 3, self.mrna.strand)
			else
				@genomic_start = to_genomic_location(@stop * 3, self.mrna.strand)
			end
		end
		# puts ":: " + @prot_acc + ">>> genomic_start\t" + @genomic_start.to_s + "\n"
		return @genomic_start
	end
	
	def genomic_stop()
		if @genomic_stop == nil then
			# puts "strand: " + self.mrna.strand
			if self.mrna.strand == "+" then
				@genomic_stop = to_genomic_location(@stop * 3, self.mrna.strand)
			else
				@genomic_stop = to_genomic_location(@start * 3 - 3, self.mrna.strand)
			end
		end
		# puts ":: " + @prot_acc + ">>> genomic_stop\t" + @genomic_stop.to_s + "\n"
		return @genomic_stop
	end
	
	def to_genomic_location (position, strand) # argument: peptide_pos*3-3 or peptide_pos*3
		if strand == "+" then
			for j in 0..self.mrna.cds_exon_starts.length - 1	
				cds_exon_length = self.mrna.cds_exon_stops[j].to_i - self.mrna.cds_exon_starts[j].to_i + 1
				# puts
				# puts "FOR CDS EXON " + j.to_s + " :\tstart: " + self.mrna.cds_exon_starts[j].to_s + "\tstop: " + self.mrna.cds_exon_stops[j].to_s
				# puts "position\t" + position.to_s + "\tto_genomic_location => cds exon length\t" + cds_exon_length.to_s + "\tof cds-exon with start pos\t" + self.mrna.cds_exon_starts[j].to_s + "\tand stop pos\t" + self.mrna.cds_exon_stops[j].to_s + "\n"
				if cds_exon_length - position > 0 then # check if the position fits in this exon
					genomic_pos = self.mrna.cds_exon_starts[j].to_i + position # genomic_pos is the sum of the start of the exon with the part that didn't fit in the previous exons
					# puts "genomic_pos IS in the exon ...\t" + genomic_pos.to_s + "\n"
					return genomic_pos
				else
					position -= cds_exon_length # deduct cds_exon_length from the position if cds_exon_length doesn't fit in the exon
				end
			end
		else
			for m in (self.mrna.cds_exon_stops.length - 1).downto(0)
				cds_exon_length = self.mrna.cds_exon_stops[m].to_i - self.mrna.cds_exon_starts[m].to_i + 1
				# puts
				# puts "FOR CDS EXON " + m.to_s + " :\tstart: " + self.mrna.cds_exon_starts[m].to_s + "\tstop: " + self.mrna.cds_exon_stops[m].to_s
				# puts "position\t" + position.to_s + "\tto_genomic_location => cds exon length\t" + cds_exon_length.to_s + "\tof cds-exon with start pos\t" + self.mrna.cds_exon_starts[m].to_s + "\tand stop pos\t" + self.mrna.cds_exon_stops[m].to_s + "\n"
				if cds_exon_length - position > 0 then # check if the position fits in this exon
					genomic_pos = self.mrna.cds_exon_stops[m].to_i - position # genomic_pos is the deduction of the part that didn't fit in the previous exons from the start of the exon
					# puts "genomic_pos IS in the exon ...\t" + genomic_pos.to_s + "\n"
					return genomic_pos
				else
					position -= cds_exon_length # deduct cds_exon_length from the position if cds_exon_length doesn't fit in the exon
				end
			end
		end
	end

	def to_bed()
		(starts,stops) = genomic_locations()
		rgb = "0,255,0"
		exon_lengths = []
		relative_starts = []
		if self.mrna.strand == "-" then
			rgb = "255,0,0"
		end
		for i in 0..starts.length-1
			exon_lengths << (stops[i] - starts[i] + 1)
			relative_starts << starts[i] - starts.first
		end	
		bed_entry = [self.mrna.chrom, starts.first, stops.last + 1, self.peptide, 1, self.mrna.strand, starts.first, stops.last + 1, rgb, starts.count, exon_lengths.join(","), relative_starts.join(",")].join("\t")
		
		# puts "bed_entry\t" + bed_entry.to_s + "\n"
		# puts "************************************************"
		return bed_entry
	end
end