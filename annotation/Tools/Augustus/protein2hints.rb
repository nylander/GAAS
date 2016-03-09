#!/usr/bin/ruby
# == NAME
# gff2hints.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# Converts Maker-generated evidence alignments from GFF3 format to Augustus-compatbile hint format
#
# == OPTIONS
#  -h,--help::                  Show help
#  -i,--infile=INFILE::         input file
#  -o,--outfile=OUTFILE::       output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, mphoeppner@gmail.com

require 'rdoc/usage'
require 'optparse'
require 'ostruct'
require 'logger'


### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-h","--help","Display the usage information") {RDoc::usage}
opts.on("-i","--infile", "=INFILE","Input") {|argument| options.infile = argument }
opts.on("-o","--outfile", "=OUTFILE","Output") {|argument| options.outfile = argument }

opts.parse! rescue RDoc::usage('usage')

options.outfile ? output_stream = File.new(options.outfile,'w') : output_stream = $stdout

lines = IO.readlines(options.infile)

evidence = lines[0].split("\t")[2]
priority = 1
type_of_evidence = nil

if evidence == "expressed_sequence_match"
	type_of_evidence = "E"
	priority = 4
elsif evidence == "protein_match"
	type_of_evidence = "P"
	priority = 2
end

abort "No evidence type identified" unless type_of_evidence
 
lines.each do |line|

  next unless line.include?("match_part")
  
  seq_region,source,feature,start,stop,score,strand,phase,tag_values = line.strip.split("\t")
  group = tag_values.split(";").find{|e| e.include?("Parent")}.split("=")[-1].strip.gsub(/\:/, '_').gsub(/\./, '-')


   output_stream.puts "#{seq_region}\tb2h\tCDSpart\t#{start}\t#{stop}\t0\t.\t.\tgroup=#{group};src=#{type_of_evidence};pri=#{priority}"

end



