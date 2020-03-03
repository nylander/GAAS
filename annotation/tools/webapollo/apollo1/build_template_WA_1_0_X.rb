#!/usr/bin/ruby
# == NAME
# build.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A script to build a new WebApollo template installation (version 1.0 and up)
# You only need to do this once - the new_species.rb script will then copy this template
# to create new build projects. 
#
# == OPTIONS
#  -h,--help::                  Show help
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
options.clean = false
opts.on("-c","--[no]clean","Clean up project") { options.clean = true }
opts.on("-h","--help","Display the usage information") { 
	puts opts 
	exit 
}

opts.parse! 

### Usernames, passwords and locations

user = ENV['USER']
home = ENV['HOME']
build_dir = ENV['APOLLO_BUILD_DIR'] or abort "Environment variable APOLLO_BUILD_DIR not set"


config = {
	:web_apollo_source => "/big/webapollo/release/live",# The source code of webapollo
	:web_apollo_build => "#{build_dir}/template",	# The location where this WA project is to be build
}

### File targets

config_files = [ "sample_config.properties" , "sample_config.xml" , "sample_blat_config.xml", "sample_hibernate.xml" ]

### The workflow

if options.clean == true
	system("rm -Rf #{config[:web_apollo_build]}")

else
  
  	File.directory?(config[:web_apollo_source]) or abort "Could not find the reference git clone (expected: #{config[:web_apollo_source]})"
  	abort "Template already exist. Use flag '--clean' to remove!" if File.directory?(config[:web_apollo_build]) 
  
	# Create the folder where the data is to be stored
	system("mkdir -p #{config[:web_apollo_build]}")
	# Create a copy of the webapollo code for this installation
  	system("mkdir -p #{config[:web_apollo_build]}")
	system("cp -R #{config[:web_apollo_source]}/* #{config[:web_apollo_build]}")
	# Copy the template config files into place
	config_files.each do |file|
		system("cp #{config[:web_apollo_build]}/#{file} #{config[:web_apollo_build]}/#{file.gsub(/sample_/, '')}")
	end
	# Build the webapollo installation
        Dir.chdir(config[:web_apollo_build]) do
		system("./apollo deploy")
	end
	
end	

