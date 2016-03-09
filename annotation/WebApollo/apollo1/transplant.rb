#!/usr/bin/ruby
# == NAME
# transplant.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# Transplants a WebApollo installation
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

require 'optparse'
require 'ostruct'

### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-s","--species", "=SPECIES","Name of species to transplant") {|argument| options.species = argument }
opts.on("-h","--help","Display the usage information") {
  puts opts
  exit
}

opts.parse! 

home = ENV['HOME']

SERVER = "bils-web.imbim.uu.se"
USER = "root"
LOCAL_PROJECT = "/big/webapollo/projects/#{options.species}"
LOCAL_DATA = "/big/data/#{options.species}"
TRANSFER = "/big/transfer/#{options.species}"
REMOTE_PROJECT = "/databases/webapollo/#{options.species}"
REMOTE_DATA = "/databases/data/#{options.species}"

# Create local folders
system("mkdir -p #{TRANSFER}")
system("mkdir -p #{LOCAL_PROJECT}")

# Copy remote folders/files
# Genome sequence
system("scp #{USER}@#{SERVER}:#{REMOTE_PROJECT}/genome.fa #{LOCAL_PROJECT}/")
abort "No genome sequence has been transferred from remote host (expecting: genome.fa)" unless File.exists?("#{LOCAL_PROJECT}/genome.fa")

# User db, if dumped
system("scp #{USER}@#{SERVER}:#{REMOTE_PROJECT}/userdb.sql #{LOCAL_PROJECT}/")

# Data
system("scp -r #{USER}@#{SERVER}:#{REMOTE_DATA}/* #{TRANSFER}")

# Deploy local installation
system("ruby #{home}/git/code/WebApollo/new_species.rb -s #{options.species} -f #{LOCAL_PROJECT}/genome.fa")

# Copy the old data to new location
system("cp -R #{TRANSFER}/Annotations* #{LOCAL_DATA}/annotations/")
system("cp -R #{TRANSFER}/track* #{LOCAL_DATA}/data/")
system("cp -R #{TRANSFER}/blat* #{LOCAL_DATA}/annotations/")
system("cp -R #{TRANSFER}/bam #{LOCAL_DATA}/data/")

system("chown -R #{ENV['USER']}:tomcat #{LOCAL_DATA}")
