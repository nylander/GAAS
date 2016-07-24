seen = []
this_seq = nil

IO.readlines(ARGV.shift).each do |line|
	line.strip!

	seq = line.split("\t")[0]

	this_seq = seq if this_seq.nil?
	if this_seq != seq
		seen = []
		this_seq = seq
	end

	next if seen.include?(line)
	seen << line
	puts line

end
