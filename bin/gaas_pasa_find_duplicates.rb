#!/usr/bin/ruby

seen = []
this_seq = nil
counter = 0

file = File.open(ARGV.shift,"r")

while (line = file.gets)

	line.strip!
	next if line.length == 0

	seq = line.split("\t")[0]

	this_seq = seq if this_seq.nil?
	if this_seq != seq
		seen = []
		this_seq = seq
	end

	if seen.include?(line)
		counter += 1
	else
		seen << line
		puts line
	end

end

warn "#{counter} lines removed!"
