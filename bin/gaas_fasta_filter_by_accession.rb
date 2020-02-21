#!/usr/bin/ruby
# == NAME
# this_script.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# 
#
# == OPTIONS
#  -h,--help::                  Show help
#  -f,--fasta=FASTA::           FASTA file to filter
#  -k,--killlist=KILLLIST::     List of accession numbers to kill from file
#  -o,--outfile=OUTFILE::       output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, mphoeppner@gmail.com

require 'rubygems'
require 'bio'
require 'rdoc/usage'
require 'optparse'
require 'ostruct'
require 'logger'

### Define modules and classes here

### Get the script arguments and open relevant files

options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-h","--help","Display the usage information") {RDoc::usage}
opts.on("-k","--killlist", "=KILLLIST","Killlist") {|argument| options.kill = argument }
opts.on("-f","--fasta", "=FASTA","FASTA file") {|argument| options.fasta = argument }
opts.on("-o","--outfile", "=OUTFILE","Output") {|argument| options.outfile = argument }

opts.parse! rescue RDoc::usage('usage')

options.outfile ? output_stream = File.new(options.outfile,'w') : output_stream = $stdout

list = IO.readlines(options.kill).collect{|e| e.strip }

ff = Bio::FastaFormat.open(options.fasta)

ff.each_entry do |entry|

	next if list.include?(entry.definition.strip.split(" ")[0])

	output_stream.puts entry.to_s

end

output_stream.close
