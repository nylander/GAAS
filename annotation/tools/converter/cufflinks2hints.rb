#!/usr/bin/ruby
# == NAME
# cufflinks2hints.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# Converts a Cufflinks-formatted GTF file into Augustus compatible exon hints
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


IO.readlines(options.infile).each do |line|
  
  next unless line.include?("Cufflinks")
  
  elements = line.strip.split("\t")
  
  seq_region,from,to,strand = elements[0],elements[3],elements[4],elements[6]
  transcript_id = elements[-1].split(";").find{|e| e.include?("transcript_id")}.gsub(/\"/, '').split(" ")[1]
  
  output_stream.puts "#{seq_region}\tb2h\tep\t#{from}\t#{to}\t.\t#{strand}\t.\tgrp=#{transcript_id};pri=4,src=W"
  
  
end

output_stream.close



