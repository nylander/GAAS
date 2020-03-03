#!/usr/bin/ruby

require 'fileutils'

# A rake file to set up a new WebApollo instance
GENOME = ENV['genome'].downcase
abort "Must provide a genome name or horrible things will happen" unless GENOME
SEQUENCE = ENV['sequence']
abort "Must provide a genome sequence!" unless SEQUENCE
ANNOTATION = ENV['gff']
abort "Must provide an annotation!" unless ANNOTATION

WEB_APOLLO_INSTALL_DIR = ENV['install_dir'] || "/databases/webapollo/WebApollo-2014-04-03" # the directory with the installation files
WEB_APOLLO_DIR = "/usr/share/tomcat7/webapps/#{GENOME}"


# Determine where all the data will live (BAM files, annotations etc)
STORAGE_DIR = "/databases/data"

BLAT_DIR = "/databases/tools/blat"

WORK_DIR=Dir.getwd

# Some key variables:

# -- The user with admin rights to the user DB
user_db = "web_apollo_users_#{GENOME}"
user_db_admin = ENV['WA_USER_DB_ADMIN'] or abort "Environment variable WA_PGUSER not set"
user_db_admin_pw = ENV['WA_WA_USER_DB_ADMIN_PW'] or abort "Environment variable WA_PGPASS not set"

# -- The user with admin rights on the website	<= Default value during installation process
web_apollo_admin = "web_apollo_admin"		
web_apollo_admin_pw = "web_apollo_admin"

# -- The Chado DB stores modified annotations
chado_db = "nbis_chado"
chado_db_user = "chado_admin"
chado_db_user_pw = "chado_admin"

# -- Custom CSS styles we wish to include (e.g. stranded RNAseq data etc)
CSS_STRING = ".plus-cigarM {\nbackground-color: green; /* color for plus matches */\n}\n\n.minus-cigarM {\nbackground-color: blue; /* color for minus matches */\n}\n"

###############################################
#### HERE BE DRAGONS! #########################
####=Things you don't want to mess with########
###############################################

directory 'log'
directory 'scratch'
directory 'split_gff'

# rake targets (files that need to be generated)
# Only way to create targets for PSQL is to 
# write log files. 
user_db_create_log = "log/user_db_create.log"
user_db_populate_log = "log/user_db_populate.log"
user_db_add_web_user_log = "log/user_db_add_web_user.log"
user_db_add_sequence_id_log = "log/user_db_add_sequence_id.log"

extract_sequence_id = "scratch/#{SEQUENCE}.ids.txt"
##

desc 'Create the user database for this annotation project'
namespace :userdb do

	file user_db_create_log => 'log' do |f|	
		system("psql -d template1 -U #{user_db_admin} -c \"DROP DATABASE IF EXISTS #{user_db}\"")
		system("psql -L #{f} -d template1 -U #{user_db_admin} -c \"CREATE DATABASE #{user_db}\"")
	end

	file user_db_populate_log => user_db_create_log do |f|
		system("psql -L #{f} -d #{user_db} -U #{user_db_admin} -f #{WEB_APOLLO_INSTALL_DIR}/tools/user/user_database_postgresql.sql")
	end

	file user_db_add_web_user_log => user_db_populate_log do |f|
		system("perl #{WEB_APOLLO_INSTALL_DIR}/tools/user/add_user.pl -D #{user_db} -U #{user_db_admin} -P #{user_db_admin_pw} -u #{web_apollo_admin} -p #{web_apollo_admin_pw} 2> #{f}")
	end
	
	file extract_sequence_id => [ 'scratch', user_db_add_web_user_log] do |f|
		system("perl #{WEB_APOLLO_INSTALL_DIR}/tools/user/extract_seqids_from_fasta.pl -p Annotations- -i #{SEQUENCE} -o #{f}")
	end
	
	file user_db_add_sequence_id_log => extract_sequence_id do |f|
		system("perl #{WEB_APOLLO_INSTALL_DIR}/tools/user/add_tracks.pl -D #{user_db} -U #{user_db_admin} -P #{user_db_admin_pw} -t #{extract_sequence_id} 2> #{f}")	
		system("perl #{WEB_APOLLO_INSTALL_DIR}/tools/user/set_track_permissions.pl -D #{user_db} -U #{user_db_admin} -P #{user_db_admin_pw} -u #{web_apollo_admin} -t #{extract_sequence_id} -a")	
	end

	task :create => user_db_add_sequence_id_log

end

desc 'Create the Tomcat folder'

namespace :servlet do 

	
	task :make_folder => "userdb:create" do
		system("mkdir -p #{WEB_APOLLO_DIR}")
		FileUtils.cd("#{WEB_APOLLO_DIR}") {
			system("jar -xvf #{WEB_APOLLO_INSTALL_DIR}/war/WebApollo.war 2> /dev/null ")
		}
	end

	task :create => :make_folder	

end

##
config_copy = "scratch/config.xml"
hibernate_copy = "scratch/hibernate.xml"
##

namespace :config do

	file config_copy => 'servlet:create' do |f|
		parsed_content = []
		fil = File.open("#{WEB_APOLLO_DIR}/config/config.xml","r")
		while (line = fil.gets)

			if line.include?("ENTER_DATASTORE_DIRECTORY_HERE")
				line.gsub!(/ENTER_DATASTORE_DIRECTORY_HERE/, "#{STORAGE_DIR}/#{GENOME}/")
			elsif line.include?("ENTER_USER_DATABASE_JDBC_URL")
				line.gsub!(/ENTER_USER_DATABASE_JDBC_URL/, "jdbc:postgresql://localhost/web_apollo_users_#{GENOME}")
			elsif line.include?("ENTER_USER_DATABASE_USERNAME")
				line.gsub!(/ENTER_USER_DATABASE_USERNAME/, "#{user_db_admin}")
			elsif line.include?("ENTER_USER_DATABASE_PASSWORD")
				line.gsub!(/ENTER_USER_DATABASE_PASSWORD/, "#{user_db_admin_pw}")
			elsif line.include?("ENTER_PATH_TO_REFSEQS_JSON_FILE")
				line.gsub!(/ENTER_PATH_TO_REFSEQS_JSON_FILE/, "#{WEB_APOLLO_DIR}/jbrowse/data/seq/refSeqs.json")
			elsif line.include?("ENTER_ORGANISM")
				line.gsub!(/ENTER_ORGANISM/ , "#{GENOME.capitalize.split('_').join(' ')}")
			elsif line.include?("ENTER_CVTERM_FOR_SEQUENCE")
				line.gsub!(/ENTER_CVTERM_FOR_SEQUENCE/, 'sequence:contig')
			end

			parsed_content << line
		end
		fil.close

		parsed_content.compact!

		o = File.new("#{f}", "w+")
		o.puts parsed_content.join("\n")
		o.close

		system("cp #{f} #{WEB_APOLLO_DIR}/config/config.xml")
		
	end

	file hibernate_copy => config_copy do |f|
		
		parsed_content = []
		fil = File.open("#{WEB_APOLLO_DIR}/config/hibernate.xml","r")
		while (line = fil.gets)
			if line.include?("ENTER_DATABASE_CONNECTION_URL")
				line.gsub!(/ENTER_DATABASE_CONNECTION_URL/, "jdbc:postgresql://localhost/#{chado_db}")
			elsif line.include?("ENTER_USERNAME")
				line.gsub!(/ENTER_USERNAME/, "#{chado_db_user}")
			elsif line.include?("ENTER_PASSWORD")
				line.gsub!(/ENTER_PASSWORD/, "#{chado_db_user_pw}")
			end

			parsed_content << line			

		end
		fil.close
		parsed_content.compact!

		o = File.new("#{f}", "w+")
		o.puts parsed_content.join("\n")
		o.close

		system("cp #{f} #{WEB_APOLLO_DIR}/config/hibernate.xml")
	end
	
	task :parse => hibernate_copy

end

namespace :blat do

	blat_db = "#{STORAGE_DIR}/#{GENOME}/blat.2bit"
	blat_config = "#{WEB_APOLLO_DIR}/config/blat_config.xml"
	blat_config_copy = "#{WORK_DIR}/scratch/blat_config.xml"

	file blat_db => 'config:parse' do
		system("mkdir -p #{STORAGE_DIR}/#{GENOME}")
		FileUtils.cd("#{STORAGE_DIR}/#{GENOME}") {
			system("#{BLAT_DIR}/faToTwoBit #{WORK_DIR}/#{SEQUENCE} #{blat_db}")
		}
	end

	file blat_config_copy => blat_db do |f|

		lines = IO.readlines(blat_config)
		parsed_lines = []

		file = File.new("#{WORK_DIR}/scratch/blat_config.xml","w+")

		
		lines.each do |line|
			if line.include?("ENTER_PATH_TO_BLAT_BINARY")
				line.gsub!(/ENTER_PATH_TO_BLAT_BINARY/, "#{BLAT_DIR}/blat")
			elsif line.include?("ENTER_PATH_FOR_TEMPORARY_DATA")
				line.gsub!(/ENTER_PATH_FOR_TEMPORARY_DATA/, '/tmp')
			elsif line.include?("ENTER_PATH_TO_BLAT_DATABASE")
				line.gsub!(/ENTER_PATH_TO_BLAT_DATABASE/, "#{blat_db}")
			elsif line.include?("ENTER_ANY_BLAT_OPTIONS")
				line.gsub!(/ENTER_ANY_BLAT_OPTIONS/, '-minScore=100 -minIdentity=60')
			end

			parsed_lines << line.strip

		end

		file.puts parsed_lines.join("\n")

		file.close

		system("cp #{f} #{blat_config}")

	end

	task :create_config => blat_config_copy

end

namespace :data do


	task :data_dir => 'blat:create_config' do 
		system("mkdir -p #{STORAGE_DIR}/#{GENOME}")
	end

	task :dna_track_setup => :data_dir do
                system("perl #{WEB_APOLLO_DIR}/jbrowse/bin/prepare-refseqs.pl --fasta #{WORK_DIR}/#{SEQUENCE}")
        end

	task :link_data_dir => :dna_track_setup do
		system("ln -sf #{STORAGE_DIR}/#{GENOME} #{WEB_APOLLO_DIR}/jbrowse/data")
	end

	task :copy_json => :link_data_dir do
		system("cp -R #{WORK_DIR}/data/* #{STORAGE_DIR}/#{GENOME}")
	end

	task :add_plugin => :copy_json do
		system("chmod +x #{WEB_APOLLO_DIR}/jbrowse/bin/*.pl")
		system("#{WEB_APOLLO_DIR}/jbrowse/bin/add-webapollo-plugin.pl -i #{WEB_APOLLO_DIR}/jbrowse/data/trackList.json")
	end

	task :split_gff_file => ['split_gff', :add_plugin ] do
		system("perl #{WEB_APOLLO_INSTALL_DIR}/tools/data/split_gff_by_source.pl -i #{ANNOTATION} -d split_gff")
	end
	
end

namespace :gff do


	task :parse => "data:split_gff_file" do

		files = Dir["split_gff/*.gff"]

		warn files.inspect

		maker = files.find{|f| f.include?("maker") } || nil
		protein = files.find{|f| f.include?("protein_coding") } || nil

		system("chmod +x #{WEB_APOLLO_DIR}/jbrowse/bin/flatfile-to-json.pl")

		FileUtils.cd("#{WEB_APOLLO_DIR}/jbrowse") {


			if maker
				files.delete_if{|f| f.include?("maker") }
				system("bin/flatfile-to-json.pl --gff #{WORK_DIR}/split_gff/maker.gff --arrowheadClass trellis-arrowhead --getSubfeatures --subfeatureClasses '{\"wholeCDS\": null, \"CDS\":\"brightgreen-80pct\", \"UTR\": \"darkgreen-60pct\", \"exon\":\"container-100pct\"}' --cssClass container-16px --type mRNA --trackLabel maker")
			elsif protein
				files.delete_if{|f| f.include?("protein_codin") }
                               system("bin/flatfile-to-json.pl --gff #{WORK_DIR}/split_gff/protein_coding.gff --arrowheadClass trellis-arrowhead --getSubfeatures --subfeatureClasses '{\"wholeCDS\": null, \"CDS\":\"brightgreen-80pct\", \"UTR\": \"darkgreen-60pct\", \"exon\":\"container-100pct\"}' --cssClass container-16px --type mRNA --trackLabel EnsEMBLProtein")

			end

			files.each do |gff_file|
				feature = gff_file.split("/")[-1].split(".")[0]
				system("bin/flatfile-to-json.pl --gff #{WORK_DIR}/#{gff_file} --arrowheadClass webapollo-arrowhead --getSubfeatures --subfeatureClasses '{\"match_part\": \"darkblue-80pct\"}' --cssClass container-10px --trackLabel #{feature}")	
			end
		}
	end

	task :add_names => :parse do

		FileUtils.cd("#{WEB_APOLLO_DIR}/jbrowse") {
			system("bin/generate-names.pl")
		}
	end

end

css_file = "#{WEB_APOLLO_DIR}/jbrowse/data/custom.css"
json_file = "#{WEB_APOLLO_DIR}/jbrowse/data/trackList.json"
json_file_copy = "#{STORAGE_DIR}/#{GENOME}/trackList.json.copy"

desc 'Creates the custom CSS file'
namespace :css do

	file css_file => 'gff:parse' do |f|
		file = File.new(css_file,"w+")
		file.puts CSS_STRING
		file.close
	end

	file json_file_copy => css_file do 

		file = File.open(json_file,"r")
		lines = IO.readlines(json_file)
		lines[0] = "\"css\" : \"data/custom.css\",\n"
		lines.unshift("{\n")
		file.close

		copy = File.new(json_file_copy,"w+")
		copy.puts lines.join
		copy.close
				
		system("cp #{json_file_copy} #{json_file}")
	end

	task :generate_custom_css => json_file_copy
end



task :default => 'css:generate_custom_css'

task :clean do 
	
	system("rm -fR #{WEB_APOLLO_DIR}")
	system("rm -fR #{WORK_DIR}/data")
	system("rm -R log")
	system("rm -Rf #{STORAGE_DIR}/#{GENOME}")
	system("rm -R split_gff")
	system("rm -R scratch")
	system("psql -d template1 -U #{user_db_admin} -c \"DROP DATABASE IF EXISTS #{user_db}\"")

end



