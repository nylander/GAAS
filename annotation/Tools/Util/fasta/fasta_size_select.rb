#!/usr/bin/ruby
# == NAME
# fasta_size_select.rb
#
# == USAGE
#  ./fasta_size_select.rb [ -h | --help ]
#                    [ -i | --infile ] | [ -s | --size ]
# == DESCRIPTION
# A script that uses BioRuby to find fasta entries above <size>
#
# == OPTIONS
#  -h,--help::                  Show help
#  -i,--infile=INFILE::         a Fasta file 
#  -s,--size=SIZE::      		Minimum size to keep

#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, marc.hoeppner@web.de

require 'rdoc/usage'
require 'optparse'
require 'ostruct'
require 'logger'
require 'bio'

### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-h","--help","Display the usage information") {RDoc::usage}
opts.on("-i","--infile", "=INFILE","Input") {|argument| options.infile = argument }
opts.on("-s","--size", "=SIZE","Min Size") {|argument| options.size = argument.to_i }

opts.parse! rescue RDoc::usage('usage')

log = Logger.new(File.new('this_script.log', File::WRONLY | File::TRUNC | File::CREAT))
log_level = Logger::INFO # or: DEBUG, WARN, FATAL, UNKNOWN

log.info('Script this_script.rb started')
log.info('Options:')
log.info(options.to_yaml)

ff = Bio::FastaFormat.open(options.infile)

ff.each_entry do | entry|

	next if entry.naseq.length < options.size
	
	puts entry.to_s
	
end
