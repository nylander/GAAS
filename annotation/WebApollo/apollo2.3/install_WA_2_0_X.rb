#!/usr/bin/ruby
# == NAME
# build.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A script to build a new WebApollo installation (version 2.0 and up)
# You only need to do this once - Managing species will then be done directly from web interface
#
# == OPTIONS
#  -h,--help::                  Show help
#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Jacques Dainat, jacques.dainat@nbis.se

require 'optparse'
require 'ostruct'

####################################
### Define modules and classes here

def tomcat_group

        user_groups = `id`
        user_groups.include?("tomcat")

end

def parse_apollo_config(file,config)
	parsed_content = []
	f = File.open(file,"r")
	while (line = f.gets)
		if line.include?("url =")
			line.gsub!(/"jdbc:postgresql:\/\/localhost\/apollo.*"/, "\"jdbc:postgresql://localhost/#{config[:web_apollo_db]}\"")
		elsif line.include?("username =")
			line.gsub!(/<CHANGEME>/, "#{config[:pguser]}")
		elsif line.include?("password =")
			line.gsub!(/<CHANGEME>/, "#{config[:pgpass]}")
		end

		parsed_content << line			

	end
	f.close
	parsed_content.compact!

	o = File.new(file, "w+")
	o.puts parsed_content
	o.close

end

####################################################
### Get the script arguments and open relevant files

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
options.clean = false
opts.on("-n","--name", "=NAME","name of the installation") {|argument| options.name = argument }
opts.on("-c","--[no]clean","Clean up the installation") { options.clean = true }
opts.on("-t","--table","=INTEGER","Codon table used. By default it is the table 1") { |argument| options.codonTable = Integer(argument) }
opts.on("-p","--path","=PATH","Define the path to the WA folder to be installed. By default we are looking at: /big/webapollo/release/live_WA_2") { |argument| options.path = argument }
opts.on("-h","--help","Display the usage information") { 
	puts opts 
	exit 
}

opts.parse! 

raise "You are not member of the tomcat group and thus cannot deploy WebApollo!" unless tomcat_group

@name = options.name or abort 'No Webapollo installation name provided'

###################
### Usernames, passwords and locations

user = ENV['USER']
home = ENV['HOME']
build_dir = ENV['APOLLO_BUILD_DIR'] or abort "Environment variable APOLLO_BUILD_DIR not set"
user_db_admin = ENV['WA_DB_ADMIN'] or abort "Environment variable WA_PGUSER not set"
user_db_admin_pw = ENV['WA_DB_ADMIN_PW'] or abort "Environment variable WA_PGPASS not set"
tomcat_apps = "/var/lib/tomcat/webapps"			        # Deployment location



#Define the codon table to use
unless defined?(options.codonTable).nil? # If define we give the default path
	if (options.codonTable<1 || options.codonTable>25)
		puts "The --table option expects an integer between 1 and 25"
		exit
	end
	if (options.codonTable>1) #We have to add the information to the otpion file
		@addCodonString=options.codonTable
	end
end

#Define path to the WA folder
if (defined?(options.path)).nil? # If not define we give the default path
        path="/big/webapollo/release/live_WA_2"
else
        path=options.path
end

config = {
	:web_apollo_source => "#{path}",# The source code of webapollo
	:web_apollo_build => "#{build_dir}/#{@name}",	# The location where this WA project is to be build
	
	:web_apollo_db => "web_apollo_#{@name}",	# Name of the user database for this WA project	
	:pgpass => "#{user_db_admin_pw}",			# PW of the SQL WA admin
	:pguser => "#{user_db_admin}"				# Username of the SQL WA admin
}


### The workflow
if options.clean == true
	
	puts "Are you really sure to erase the whole webapollo installation ( #{tomcat_apps}/#{@name}, #{config[:web_apollo_build]},DB:#{config[:web_apollo_db]} ) ? [y|n]:"
        selection = gets.chomp
	if(selection.downcase == "y")
		puts "Lets remove everything"
		
		puts "Cleaning webapollo folder"
		system("sudo rm -Rf #{config[:web_apollo_build]}")
		puts "Cleaning database"
#		system("psql -d template1 -U #{config[:pguser]} -c \"DROP DATABASE IF EXISTS #{config[:web_apollo_db]}\"")
		system("sudo su - postgres -c \"dropdb #{config[:web_apollo_db]}\"")
		puts "Cleaning tomcat folder"
		system("rm -f #{tomcat_apps}/#{@name}.war")
		system("sudo rm -Rf #{tomcat_apps}/#{@name}")
		puts "Cleaning finished"
	else
		puts "Fine we let everything as it was."
	end
else
  
  	File.directory?(config[:web_apollo_source]) or abort "Could not find the reference folder (expected: #{config[:web_apollo_source]})"
  	abort "Installation already exist. Use flag '--clean' to remove!" if File.directory?(config[:web_apollo_build]) 
  
	# Create the folder where the data is to be stored
	puts "Create folders"
	system("mkdir -p #{config[:web_apollo_build]}")

	# Create a copy of the webapollo code for this installation
	puts "Create a copy of the webapollo source"
	system("cp -R #{config[:web_apollo_source]}/* #{config[:web_apollo_build]}")

	# Copy the template config files into place and configurate it
	puts "rename the config file and configure it properly" 
	system("cp #{config[:web_apollo_build]}/sample-postgres-apollo-config.groovy #{config[:web_apollo_build]}/apollo-config.groovy")
	parse_apollo_config("#{config[:web_apollo_build]}/apollo-config.groovy",config)
	if(defined?(addCodonString).nil?)
		open("#{config[:web_apollo_build]}/apollo-config.groovy", 'a') { |f|
		  f << "// default apollo settings\n"
		  f << "apollo {\n"
		  f << "get_translation_code = #{@addCodonString}\n"
		  f << "}\n"
		}
	end	
	
	#Create a new database for this installation
	puts "Create the database #{config[:web_apollo_db]}"
	system("sudo su - postgres -c \"createdb -E UTF-8 -O #{config[:pguser]} #{config[:web_apollo_db]}\"")


	# Build the webapollo installation
	Dir.chdir(config[:web_apollo_build]) do
		system("./apollo deploy")
	end
	
	# Copy the packaged WebApollo installation to the Tomcat folder
	system("cp #{config[:web_apollo_build]}/target/apollo*.war #{tomcat_apps}/#{@name}.war")

	#Final message:
	puts "Congratulation ! Installation done. You can retrieve your installation on the webpage http://annotation-prod.scilifelab.se:8080/#{@name}"
	
end	

