input = File.open(ARGV.shift,"r")

header = "track type=WIG name=\"FST Wiggle file\" description=\"FST data\""
puts header

current_scaffold = "none"

last_pos = 1

while (line =input.gets)

	next unless line.include?("scaffold")

	scaffold,start,stop,score = line.strip.split("\t")

	next if scaffold == current_scaffold && start.to_i == last_pos

	last_pos = stop.to_i

	if scaffold != current_scaffold
		puts "variableStep chrom=#{scaffold}"
		current_scaffold = scaffold
	end

	
	puts "#{start}\t#{score}"
	#puts "#{stop}\t#{score}"


end

input.close
