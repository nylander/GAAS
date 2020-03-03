
tcounter = 10000
ccounter = 10000

this_mrna = nil

lines = IO.readlines(ARGV.shift)

lines.each do |line|

	elements = line.strip.split("\t")

	if line.include?("EuGene") && line.include?("mRNA") || line.include?("CDS")
		
		elements = line.strip.split("\t")
		features = {}
		elements[-1].split(";").each{|e| features[e.split("=")[0]] = e.split("=")[1]}


		if elements[2] == "mRNA"
			tcounter += 1
			puts "#{elements[0]}\teugene\tmatch\t#{elements[3]}\t#{elements[4]}\t#{elements[5]}\t#{elements[6]}\t.\tID=#{features['ID']};Name=Transcript#{tcounter}"
			this_mrna = features["ID"]
		elsif elements[2] = "CDS"
			ccounter += 1
			puts "#{elements[0]}\teugene\tmatch_part\t#{elements[3]}\t#{elements[4]}\t#{elements[5]}\t#{elements[6]}\t.\tID=#{features['ID']};Name=CDS#{ccounter},Parent=#{this_mrna}"
		end
			
	end

end
