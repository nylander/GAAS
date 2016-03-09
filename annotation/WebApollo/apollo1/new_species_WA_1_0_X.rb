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
#  -g,--gff=GFF::       	Annotation file to process (optional)
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

def parse_config_xml(file,config)

	parsed_content = []
	f = File.open(file,"r")
	while (line = f.gets)

		if line.include?("ENTER_DATASTORE_DIRECTORY_HERE")
			line.gsub!(/ENTER_DATASTORE_DIRECTORY_HERE/, config[:data_dir])
		elsif line.include?("jdbc:postgresql:web_apollo_users")
			line.gsub!(/jdbc:postgresql:web_apollo_users/, "jdbc:postgresql://localhost/#{config[:web_apollo_user_db]}")
		elsif line.include?("ENTER_USER_DATABASE_USERNAME")
			line.gsub!(/ENTER_USER_DATABASE_USERNAME/, config[:web_apollo_admin])
		elsif line.include?("ENTER_USER_DATABASE_PASSWORD")
			line.gsub!(/ENTER_USER_DATABASE_PASSWORD/, config[:web_apollo_admin_pw])
		elsif line.include?("ENTER_PATH_TO_REFSEQS_JSON_FILE")
			line.gsub!(/ENTER_PATH_TO_REFSEQS_JSON_FILE/, "#{config[:data_dir]}/refSeqs.json")
		elsif line.include?("ENTER_ORGANISM")
			line.gsub!(/ENTER_ORGANISM/ , config[:organism] )
		elsif line.include?("/config/translation_tables")
			line.gsub!(/ncbi_1_/, "ncbi_#{config[:translation_table]}_")
		end

		parsed_content << line
	end
	f.close

	parsed_content.compact!

	o = File.new(file, "w+")
	o.puts parsed_content.join("\n")
	o.close

end

def parse_hibernate_xml(file,config)

	parsed_content = []
	f = File.open(file,"r")
	while (line = f.gets)
		if line.include?("ENTER_DATABASE_CONNECTION_URL")
			line.gsub!(/ENTER_DATABASE_CONNECTION_URL/, "jdbc:postgresql://localhost/#{config[:chado_db]}")
		elsif line.include?("ENTER_USERNAME")
			line.gsub!(/ENTER_USERNAME/, "#{config[:chado_db_user]}")
		elsif line.include?("ENTER_PASSWORD")
			line.gsub!(/ENTER_PASSWORD/, "#{config[:chado_db_user_pw]}")
		end

		parsed_content << line			

	end
	f.close
	parsed_content.compact!

	o = File.new(file, "w+")
	o.puts parsed_content.join("\n")
	o.close

end

def parse_blat_xml(file,config)


	lines = IO.readlines(file)
	parsed_lines = []

	lines.each do |line|
		if line.include?("ENTER_PATH_TO_BLAT_BINARY")
			line.gsub!(/ENTER_PATH_TO_BLAT_BINARY/, "#{config[:tool_dir]}/blat")
		elsif line.include?("ENTER_PATH_FOR_TEMPORARY_DATA")
			line.gsub!(/ENTER_PATH_FOR_TEMPORARY_DATA/, '/tmp')
		elsif line.include?("ENTER_PATH_TO_BLAT_DATABASE")
			line.gsub!(/ENTER_PATH_TO_BLAT_DATABASE/, "#{config[:annotation_dir]}/blat.2bit")
		elsif line.include?("ENTER_ANY_BLAT_OPTIONS")
			line.gsub!(/ENTER_ANY_BLAT_OPTIONS/, '-minScore=100 -minIdentity=60')
		end

		parsed_lines << line.strip

	end
	
	f = File.new(file,"w+")
	f.puts parsed_lines.join("\n")
	f.close

end

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
options.clean = false
opts.on("-s","--species", "=SPECIES","Species name") {|argument| options.species = argument }
opts.on("-g","--gff", "=GFF","Annotation") {|argument| options.gff = argument }
opts.on("-f","--fasta", "=FASTA","Fasta file") {|argument| options.fasta = argument }
opts.on("-c","--[no]clean","Clean up project") { options.clean = true }
opts.on("-t","--translation_table","=TRANSLATION_TABLE","NCBI translation table to use") { |argument| options.translation_table = argument }
opts.on("-h","--help","Display the usage information") { 
	puts opts 
	exit 
}

opts.parse! 

raise "You are not member of the tomcat group and thus cannot deploy WebApollo!" unless tomcat_group

@species = options.species or abort 'No species name provided'
@organism = options.species.split('_')[0].capitalize + ' ' + options.species.split('_')[-1]
@fasta = options.fasta or abort 'No genome sequence provided'
options.translation_table ? @translation_table = options.translation_table : @translation_table = "1"

# Custom CSS styles needed for WA
CSS_STRING = ".plus-cigarM {\nbackground-color: green; /* color for plus matches */\n}\n\n.minus-cigarM {\nbackground-color: blue; /* color for minus matches */\n}\n"

### Usernames, passwords and locations

user = ENV['USER']
home = ENV['HOME']
user_db_admin = ENV['WA_USER_DB_ADMIN'] or abort "Environment variable WA_PGUSER not set"
user_db_admin_pw = ENV['WA_WA_USER_DB_ADMIN_PW'] or abort "Environment variable WA_PGPASS not set"
wa_website_username = ENV['WA_WEBSITE_USERNAME'] or abort "Environment variable WA_WEBSITE_USERNAME not set"
wa_website_password = ENV['WA_WEBSITE_PASSWORD'] or abort "Environment variable WA_WEBSITE_PASSWORD not set"

build_dir = ENV['APOLLO_BUILD_DIR'] or abort "Environment variable APOLLO_BUILD_DIR not set"
data_dir = ENV['APOLLO_DATA_DIR'] or abort "Environment vairable APOLLO_DATA_DIR not set"

web_apollo_storage = "#{data_dir}/#{@species}"		# Folder tree where data is stored
tomcat_apps = "/var/lib/tomcat/webapps"			        # Deployment location

config = {

	:pguser => "#{user_db_admin}",			# Username of the SQL WA admin
	:pgpass => "#{user_db_admin_pw}",			# PW of the SQL WA admin
	:web_apollo_admin => "#{wa_website_username}",		# Username of admin for WA website
	:web_apollo_admin_pw => "#{wa_website_password}",			# Password of WA website admin
	:web_apollo_user_db => "web_apollo_users_#{@species}",	# Name of the user database for this WA project
	:web_apollo_source => "#{build_dir}/template_march2015",    	# The pre-built template install
	:web_apollo_build => "#{build_dir}/#{@species}",	# The location where this WA project is to be build
	:web_apollo_storage => "#{web_apollo_storage}",      	# Location where data is stored
	:annotation_dir => "#{web_apollo_storage}/annotations",	# Where annotation data is stored
	:translation_table => @translation_table,		# NCBI translation table to use
	:data_dir => "#{web_apollo_storage}/data",		# Where the Jbrowse data is stored
	:tool_dir => "/opt/ucsc",			# Location of blat and FaToNib
	:organism => @organism,					# Species name
	:chado_db => "bils_chado",				# Name of the chado db to use for storing annotations (optional)
	:chado_db_user => "bils_chado_user",			# User with write permissions in chado db
	:chado_db_user_pw => "bils_chado_user"			# PW for chado user	
}

### File targets

config_files = [ "sample_config.properties" , "sample_config.xml" , "sample_blat_config.xml", "sample_hibernate.xml" ]

### The workflow

if options.clean == true
  	system("rm -Rf #{config[:web_apollo_storage]}")
	system("rm -Rf #{config[:web_apollo_build]}")
	system("rm -R temp")
	system("psql -d template1 -U #{config[:pguser]} -c \"DROP DATABASE IF EXISTS #{config[:web_apollo_user_db]}\"")

else
	File.directory?(config[:web_apollo_source]) or raise "No template installation found!"

	# Create the folder where the data is to be stored
	system("mkdir -p #{config[:web_apollo_storage]}")
	system("mkdir -p #{config[:annotation_dir]}")
	system("mkdir -p #{config[:data_dir]}")
	system("chgrp -R tomcat #{config[:web_apollo_storage]}") # Must be owned by tomcat group
	system("mkdir -p #{config[:web_apollo_build]}")
	# Create a temporary folder
	system("mkdir -p temp")
	# Create a copy of the webapollo template code for this installation
  	system("mkdir -p #{config[:web_apollo_build]}")
	system("cp -R #{config[:web_apollo_source]}/* #{config[:web_apollo_build]}")
	
	# Create a new database for this genome
	system("psql -d template1 -U #{config[:pguser]} -c \"DROP DATABASE IF EXISTS #{config[:web_apollo_user_db]}\"")
  	system("psql -d template1 -U #{config[:pguser]} -c \"CREATE DATABASE #{config[:web_apollo_user_db]}\"")
	# Load the database schema
	system("psql -d #{config[:web_apollo_user_db]} -U #{config[:pguser]} -f #{config[:web_apollo_source]}/tools/user/user_database_postgresql.sql")
	# Add website admin user
	system("#{config[:web_apollo_build]}/tools/user/add_user.pl -D #{config[:web_apollo_user_db]} -U #{config[:pguser]} -P #{config[:pgpass]} -u #{config[:web_apollo_admin]} -p #{config[:web_apollo_admin_pw]}") 
	# Load genome assembly
	system("#{config[:web_apollo_build]}/tools/user/extract_seqids_from_fasta.pl -p Annotations- -i #{@fasta} -o temp/seqids.txt")
	system("#{config[:web_apollo_build]}/tools/user/add_tracks.pl -D #{config[:web_apollo_user_db]} -U #{config[:pguser]} -P #{config[:pgpass]} -t temp/seqids.txt")
	system("#{config[:web_apollo_build]}/tools/user/set_track_permissions.pl -D #{config[:web_apollo_user_db]} -U #{config[:pguser]} -P #{config[:pgpass]} -u #{config[:web_apollo_admin]} -t temp/seqids.txt -a")
	system("#{config[:web_apollo_build]}/bin/prepare-refseqs.pl --fasta #{@fasta} --out #{config[:data_dir]}")
	# LOADING ANNOTATIONS ?

	# Add webapollo plugin to the genome browser
	system("#{config[:web_apollo_build]}/client/apollo/bin/add-webapollo-plugin.pl -i #{config[:data_dir]}/trackList.json")
	
	# Update the config files
	system("echo jbrowse.data=#{config[:data_dir]} > #{config[:web_apollo_build]}/config.properties")
	system("echo datastore.directory=#{config[:annotation_dir]} >> #{config[:web_apollo_build]}/config.properties")
	system("echo database.url=jdbc:postgresql://localhost/#{config[:web_apollo_user_db]} >> #{config[:web_apollo_build]}/config.properties")
	system("echo database.username=#{config[:pguser]} >> #{config[:web_apollo_build]}/config.properties")
	system("echo database.password=#{config[:pgpass]} >> #{config[:web_apollo_build]}/config.properties")
	system("echo organism=#{@organism} >> #{config[:web_apollo_build]}/config.properties")
	system("echo tracks.refseqs=#{config[:data_dir]}/seq/refSeqs.json >> #{config[:web_apollo_build]}/config.properties ")
	system("echo tracks.data=#{config[:data_dir]} >> #{config[:web_apollo_build]}/config.properties")
	
	parse_config_xml("#{config[:web_apollo_build]}/config.xml",config)
	parse_hibernate_xml("#{config[:web_apollo_build]}/hibernate.xml",config)
	parse_blat_xml("#{config[:web_apollo_build]}/blat_config.xml",config)

  	# Create custom CSS style sheet
  	f = File.new("#{config[:data_dir]}/custom.css","w")
  	f.puts CSS_STRING
  	f.close
  
	# Build Blat database
	system("#{config[:tool_dir]}/faToTwoBit #{@fasta} #{config[:annotation_dir]}/blat.2bit")

	# Build the webapollo war file
  	Dir.chdir(config[:web_apollo_build]) do
    		system("./apollo deploy")
  	end
	
	# Copy the packaged WebApollo installation to the Tomcat folder
	system("cp #{config[:web_apollo_build]}/target/apollo*.war #{tomcat_apps}/#{@species}.war")
	
end	

