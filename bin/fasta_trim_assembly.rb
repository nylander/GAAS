require 'bio'

fasta = Bio::FastaFormat.open(ARGV.shift,"r")

fasta.each_entry do |entry|

	leading = entry.naseq.slice(/^n*/).length
	trailing = entry.naseq.reverse.slice(/^n*/).length
	length = entry.nalen
	
	

	puts entry.naseq.subseq(leading+1,length-(trailing)).to_fasta(entry.definition,80)

end
