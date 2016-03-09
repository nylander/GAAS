IO.readlines(ARGV.shift).each do |line|

	if line.match(/^\#.*$/)
		puts line
	else

		e = line.strip.split("\t")

		next if e.length < 7

		phase_data = e[7]

		new_phase = case phase_data
			when "." then "."
			when "0" then "0"
			when "1" then "2"
			when "2" then "1"
		end

		e[7] = new_phase

		puts e.collect{|entry| entry + "\t" }.join.strip
	end
end
