#!/usr/bin/ruby
# == NAME
# sync_user_db.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# Parses a WebApollo SQL file and updates the local user database
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
opts.on("-i","--infile", "=INFILE","Input") {|argument| options.infile = argument }
opts.on("-s","--species", "=SPECIES","Species to use") {|argument| options.species = argument }
opts.on("-h","--help","Display the usage information") {
  puts opts
  exit
}

opts.parse! 

abort "Must provide species name!" unless options.species
abort "Must provide an input file!" unless options.infile

home = ENV['HOME']

WEB_APOLLO_SCRIPT="/big/webapollo/build/#{options.species}/tools/user/add_user.pl"
# "#{home}/webapollo/build/template/tools/user/add_user.pl"

options.infile ? input_stream = File.open(options.infile) : input_stream = $stdin

while (line = input_stream.gets)

	next unless line.match(/^[0-9]/)
	
	user_id,name,pass = line.strip.split("\t")
	
	next if user_id == "1"

	warn "Adding user #{name} with password #{pass} to database web_apollo_users_#{options.species}"

	system("perl #{WEB_APOLLO_SCRIPT} -D web_apollo_users_#{options.species} -U web_apollo_users_admin -P web_apollo_users_admin -u #{name} -p #{pass}")
	
end



