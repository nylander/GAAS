#!/usr/bin/ruby
# == NAME
# rfam2apollo.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# Converts output from Rfam/Infernal to full GFF3 format for WebApollo
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

id_counter = 100000

lines.each do |line|

  id_counter += 1
  
  e = line.strip.split("\t")
  
  rfam_features = {}
  e[-1].split(";").collect{|f| rfam_features[f.split("=")[0]] = f.split("=")[1] }
  
  output_stream.puts "#{e[0]}\t#{e[1]}\tgene\t#{e[3]}\t#{e[4]}\t#{e[5]}\t#{e[6]}\t#{e[7]}\t#{e[8]}"
  output_stream.puts"#{e[0]}\t#{e[1]}\tncRNA\t#{e[3]}\t#{e[4]}\t#{e[5]}\t#{e[6]}\t#{e[7]}\tID=RFAMT#{id_counter};Parent=#{rfam_features['ID']};Name=#{rfam_features['ID']};description=#{rfam_features['rfam-acc']} (#{rfam_features['rfam-id']})"
  output_stream.puts"#{e[0]}\t#{e[1]}\texon\t#{e[3]}\t#{e[4]}\t#{e[5]}\t#{e[6]}\t#{e[7]}\tID=RFAME#{id_counter};Parent=RFAMT#{id_counter}"
  
end

output_stream.close



