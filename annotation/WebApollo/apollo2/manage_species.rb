#!/usr/bin/ruby
# == NAME
# build.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A script to build a new WebApollo installation based on pre-built template
#
# == OPTIONS
#  -h,--help::                  Show help
#  -s,--species=SPECIES::       Name of the species
#  -f,--fasta=FASTA::		Genome sequence 
#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'

### Define modules and classes here

def tomcat_group

        user_groups = `id`
        user_groups.include?("tomcat")

end

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
options.clean = false
opts.on("-s","--species", "=SPECIES","Species name") {|argument| options.species = argument }
opts.on("-f","--fasta", "=FASTA","Fasta file") {|argument| options.fasta = argument }
opts.on("-c","--[no]clean","Clean up project") { options.clean = true }
opts.on("-w","--wa_installation", "=wa_installation_name","path to the WA installation") {|argument| options.wa_installation_name = argument }
opts.on("-h","--help","Display the usage information") { 
	puts opts 
	exit 
}

opts.parse! 

raise "You are not member of the tomcat group and thus cannot deploy WebApollo!" unless tomcat_group

@species = options.species or abort 'No species name provided'
@organism = options.species.split('_')[0].capitalize + ' ' + options.species.split('_')[-1]
@wa_installation_name = options.wa_installation_name or abort 'Name of the WebApollo installation not provided'
if options.clean != true
	@fasta = options.fasta or abort 'No genome sequence provided'
end

# Custom CSS styles needed for WA
CSS_STRING = ".plus-cigarM {\nbackground-color: green; /* color for plus matches */\n}\n\n.minus-cigarM {\nbackground-color: blue; /* color for minus matches */\n}\n"

### Usernames, passwords and locations

user = ENV['USER']
home = ENV['HOME']

build_dir = ENV['APOLLO_BUILD_DIR'] or abort "Environment variable APOLLO_BUILD_DIR not set"
data_dir = ENV['APOLLO_DATA_DIR'] or abort "Environment vairable APOLLO_DATA_DIR not set"

web_apollo_storage = "#{data_dir}/#{@species}"		# Folder tree where data is stored

config = {
	:web_apollo_path => "#{build_dir}/#{@wa_installation_name}",	# The location where this WA installation has been built
	:tool_dir => "/opt/ucsc",			# Location of blat and FaToNib
}

### The workflow

if options.clean == true
        puts "Are you really sure to remove the species folder ? [y|n]:"
	selection = gets.chomp
        if(selection.downcase == "y")
		puts "Cleaning database"
		system("#{config[:web_apollo_path]}/docs/web_services/examples/groovy/delete_annotations_from_organism.groovy  -destinationurl http://localhost:8080/#{@wa_installation_name} -organismname #{@species}")
	
		puts "Cleaning webapollo folder"
  		system("rm -Rf #{web_apollo_storage}")

		puts "Cleaning finished"
	else
		puts "Fine we let everythong as it was."
	end
else
	#check if instalation already existing
	if File.directory?("#{web_apollo_storage}") 
		puts "Instalation for species #{@species} already exits !"
		exit
	else
		# Create the folder where the data is to be stored
		puts "Create folders"
		system("mkdir -p #{web_apollo_storage}")
		system("chgrp -R tomcat #{web_apollo_storage}") # Must be owned by tomcat group
		
		# Load genome assembly
		puts "Loading genome assembly"
		system("#{config[:web_apollo_path]}/bin/prepare-refseqs.pl --fasta #{@fasta} --out #{web_apollo_storage}")


	  	# Create custom CSS style sheet
	  	f = File.new("#{web_apollo_storage}/custom.css","w")
	  	f.puts CSS_STRING
	  	f.close
	  
		# Build Blat database
		puts "Build Blat database"
		system("#{config[:tool_dir]}/faToTwoBit #{@fasta} #{web_apollo_storage}/blat.2bit")
	
	end
end	

