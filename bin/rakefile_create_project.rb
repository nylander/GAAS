#!/usr/bin/ruby
# Loads a genome assembly into a new EnsEMBL database
# = PREREQS:
# - contig2chromosome AGP file (.contigs.agp)
# - Contig fasta file 
# = OPTIONAL:
# - Supercontig2contig AGP file (.supercontigs.agp)
# = NOTE:
# Additional analyses need to be added to this script as
# needed! (last part of the setup)
# = USAGE
# 0. Check the settings in this script, especially the analyses/rules

# 3. rake submission:create genome=my_genome_name
# 4. Manually edit configs/pipeline-configs/modules/Bio/EnsEMBL/Pipelines/Config/BatchQueue.pm
# 5. rake sanity_check genome=my_genome_name
# 6. submit with rulemanager.pl

GENOME = ENV['genome'] ||= 'test'
CONTIGS = FileList["*contigs.fa"][0]
CONTIG_AGP_FILE =  FileList['*contigs.agp'][0]

raise "Contig Fasta file is missing" unless CONTIGS
raise "Contig AGP file is missing!" unless CONTIG_AGP_FILE

# The supercontig is optional, let's account for that.
supercontig_file = FileList['*supercontig*.agp']
supercontig_file.length > 0 ? SUPERCONTIG_AGP_FILE = supercontig_file : SUPERCONTIG_AGP_FILE = nil

## TO DO LIST ###

# Run Tests and actual analyses (?)
# Figure out how to set PERL5LIB from script

################################
# CONFIGURATION OF ENVIRONMENT #
################################

MYSQL_ROOT = "ensembl_admin"
MYSQL_ROOT_PW = "annotationadmin"
MYSQL_USER = "ensembl_user"
MYSQL_PW = "annotation"
MYSQL_PORT = 3306
MYSQL_HOST = ENV['host'] || 'bdb'

ENSEMBL_PATH = "/sw/ensembl/ensembl-live"

DUMP_PATH = "dump/"

################################
# PRODUCTION DATABASE ##########
################################

# The Production DB is a fully finished database from which we can steal meta data
# Ideally, we use human for this.

PRODUCTION_DB="homo_sapiens_core_74_37"

################################
# LOCATION OF REFERENCE FILES ##
################################

protein_db = "/projects/references/databases/uniprot/reference_proteomes/caenorhabditis_elegans.fa"
uniprot_db = "/projects/references/databases/uniprot/uniprot_sprot.clean"

################################
# LOCATION OF BINARIES##########
################################

repeat_masker = "/sw/bioinfo/RepeatMasker/RepeatMasker"
pmatch = "pmatch"

def command?(command)
  system("which #{ command} > /dev/null 2>&1")
end
 
raise "Missing RepeatMasker binary at #{repeat_masker}" unless command?(repeat_masker)

################################
# Quick API for short checks ###
################################

require 'active_record'

module EnsemblDB 
  
  include ActiveRecord
  
  class DBConnection < ActiveRecord::Base

		self.abstract_class = true
		self.pluralize_table_names = false
  	
		def self.connect(args={})
			establish_connection(
				:adapter => 'mysql',
				:host => MYSQL_HOST,
				:database => "ensembl_#{GENOME}",
				:username => args[:username] ||= MYSQL_USER,
				:password => args[:password] ||= MYSQL_PW
			)
		end
	end
  
	class InputIdAnalysis < DBConnection
		self.primary_key = 'input_id'	
	end

end

######################################
# CONFIGURATION OF ANALYSES/RULES ####
######################################


analyses = [ 
  
  "[SubmitContig]" ,
	"input_id_type=CONTIG" ,
	"module=Dummy",
  "" ,    		
  "[repeatmask]" ,
  "db=repbase" ,
	"db_version=0129" ,
	"db_file=repbase" , 
	"program=RepeatMasker" ,
	"program_version=3.1.8" ,
	"program_file=#{repeat_masker}" ,
	"parameters=-nolow -species nematoda -s" ,
	"module=RepeatMasker" ,
	"gff_source=RepeatMasker" ,
	"gff_feature=repeat" ,
	"input_id_type=CONTIG" ,
	"" ,
	"[genscan]",
	"db=HumanIso.smat",
	"db_file=/sw/bioinfo/genscan/HumanIso.smat",
	"program_file=/sw/bioinfo/genscan/genscan",
	"module=Genscan",
	"input_id_type=CONTIG",
	"",
	"[trnascan]",
	"db=trna",
	"program_file=/sw/bioinfo/tRNAscan-1.3.1/bin/tRNAscan-SE",
	"module=tRNAscan_SE",
	"parameters=-G",
	"input_id_type=CONTIG",
	"",
	"[uniprot]" ,
	"db=uniprot",
	"db_file=#{uniprot_db}" ,
	"program=blastall" ,
	"program_file=blastall",
	"module=BlastGenscanPep",
	"parameters=",
	"input_id_type=CONTIG",
	"",
	"[blast_wait]",
	"module=Accumulator",
	"input_id_type=ACCUMULATOR",
	"",
	
]

rules = [
  "[repeatmask]" , 
	"condition=SubmitContig",
	"",
	"[trnascan]",
	"condition=SubmitContig",
	"",
	"[genscan]",
	"condition=SubmitContig",
	"",
	"[uniprot]" ,
	"condition=genscan" ,
	"" ,
	"[blast_wait]",
	"condition=uniprot",
	"",

]

###############
# Directories #
###############

directory 'log'
directory 'dump'
directory 'analyses_and_rules'
directory 'configs/pipeline-configs/modules/Bio/EnsEMBL/'
directory 'output'

################
# RAKE targets #
################

# Rake targets for :environment
db_create_log = "log/sql.#{GENOME}.sql"
db_coretable_log = "log/sql.coretables.#{GENOME}.sql"
db_pipelinetable_log = "log/sql.pipelinetables.#{GENOME}.sql"

# Rake targets for :assembly
db_chromosome_level_log = "log/perl.assembly.chromosome.#{GENOME}.log"
db_supercontig_level_log = "log/perl.assembly.supercontig.#{GENOME}.log"
db_contig_level_log = "log/perl.assembly.contig.#{GENOME}.log"
db_chromosome2contig_log = "log/perl.assembly.chromosome2contig.map.#{GENOME}.log"
db_supercontig2contig_log = "log/perl.assembly.supercontig2contig.map.#{GENOME}.log"
db_set_attributes_log = "log/perl.assembly.set_attributes.#{GENOME}.log"
db_set_toplevel_log = "log/perl.assembly.set_toplevel.#{GENOME}.log"

dumped_sequences_log = "log/dumped_seq.toplevel.#{GENOME}.log"
dumped_database = "dump/database.#{GENOME}.dump"
load_genewise_data = "log/mysql.load_genewise.log"


#####################################
# START OF PIPELINE! ################
#####################################

desc 'Creates the blank database and tables'
namespace :environment do

	file db_create_log => ['log', CONTIGS ] do 
		warn "---------------------------------------------------"
		warn "Initializing EnsEMBL database ensembl_#{GENOME}"
		warn "Initializing analysis database genewise_#{GENOME}"
		warn "---------------------------------------------------"
		
		# 1. Drop DBs if exists
		system("mysql -u#{MYSQL_ROOT} -h#{MYSQL_HOST} -e 'DROP DATABASE IF EXISTS ensembl_#{GENOME}'")
		system("mysql -u#{MYSQL_ROOT} -h#{MYSQL_HOST} -e 'DROP DATABASE IF EXISTS genewise_#{GENOME}'")

    # 2. Create DBs
		system("mysql -u#{MYSQL_ROOT} -h#{MYSQL_HOST} -e 'CREATE DATABASE ensembl_#{GENOME}; GRANT SELECT,INSERT,UPDATE,DELETE ON ensembl_#{GENOME}.* TO ensembl_user;' > #{db_create_log}")
		system("mysql -u#{MYSQL_ROOT} -h#{MYSQL_HOST} -e 'CREATE DATABASE genewise_#{GENOME}; GRANT SELECT,INSERT,UPDATE,DELETE ON genewise_#{GENOME}.* TO ensembl_user;' >> #{db_create_log}")
		
	end
	
	file db_coretable_log => db_create_log do
	
		warn "---------------------------------------------------"
		warn "Loading CORE tables into database ensembl_#{GENOME}"
		warn "---------------------------------------------------"
		
		# 1. Load CORE tables
		system("mysql -D ensembl_#{GENOME} -u#{MYSQL_USER} -p#{MYSQL_PW} -h #{MYSQL_HOST} < #{ENSEMBL_PATH}/ensembl/sql/table.sql > #{db_coretable_log}")
		system("mysql -D genewise_#{GENOME} -u#{MYSQL_USER} -p#{MYSQL_PW} -h #{MYSQL_HOST} < #{ENSEMBL_PATH}/ensembl/sql/table.sql >> #{db_coretable_log}")	
		
	end	
	
	file db_pipelinetable_log => db_coretable_log do
		
		warn "-------------------------------------------------------"
		warn "Loading pipeline tables into database ensembl_#{GENOME}"
		warn "-------------------------------------------------------"
		
		# 1. Load Pipeline tables
		system("mysql -D ensembl_#{GENOME} -u#{MYSQL_USER} -p#{MYSQL_PW} -h #{MYSQL_HOST} < #{ENSEMBL_PATH}/ensembl-pipeline/sql/table.sql > #{db_pipelinetable_log}")
		system("mysql -D genewise_#{GENOME} -u#{MYSQL_USER} -p#{MYSQL_PW} -h #{MYSQL_HOST} < #{ENSEMBL_PATH}/ensembl-pipeline/sql/table.sql >> #{db_pipelinetable_log}")
		
		# Fix Engine of job table
		system("mysql -u#{MYSQL_ROOT} -h#{MYSQL_HOST} -D ensembl_#{GENOME} -e 'ALTER TABLE job_status ENGINE=MyISAM;'")
		system("mysql -u#{MYSQL_ROOT} -h#{MYSQL_HOST} -D genewise_#{GENOME} -e 'ALTER TABLE job_status ENGINE=MyISAM;'")
		
		# Fix the table schema to be compatible with Apollo
		system("mysql -u#{MYSQL_USER} -p#{MYSQL_PW} -h#{MYSQL_HOST} -D ensembl_#{GENOME} -e 'create view gene_stable_id AS select gene_id,stable_id,version,created_date,modified_date FROM gene;create view transcript_stable_id AS SELECT transcript_id,stable_id,version,created_date,modified_date FROM transcript;create view exon_stable_id AS SELECT exon_id,stable_id,version,created_date,modified_date FROM exon;create view translation_stable_id AS SELECT translation_id,stable_id,version,created_date,modified_date FROM translation;'")
    
		# 2. Report back the schema version of EnsEMBL used for initialization
		warn "-------------------------------------"
		warn " Schema version used for this Db is: "
		system("mysql -D ensembl_#{GENOME} -u#{MYSQL_USER} -p#{MYSQL_PW} -h #{MYSQL_HOST} -e 'select meta_value from meta where meta_key = \"schema_version\";'")
	
		
	end
	
	desc 'Creates a new EnsEMBL database'
	task :db => db_pipelinetable_log
	
end

########################################
# LOAD ASSEMBLY INTO DB ################
########################################

desc 'Loads the assembly files into the DB'
namespace :assembly do

	file db_chromosome_level_log => db_pipelinetable_log do
		warn "---------------------------------------------------"
		warn "Loading chromosomes into database ensembl_#{GENOME}"
		warn "---------------------------------------------------"
		
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/load_seq_region.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -coord_system_name chromosome -coord_system_version v#{GENOME} -rank 1 -default_version -agp_file #{CONTIG_AGP_FILE} > #{db_chromosome_level_log}")
	end
	
	file db_supercontig_level_log => db_chromosome_level_log do
	
		if SUPERCONTIG_AGP_FILE
			warn "---------------------------------------------------"
			warn "Loading supercontigs into database ensembl_#{GENOME}"
			warn "---------------------------------------------------"
			
			system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/load_seq_region.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -coord_system_name supercontig -coord_system_version v#{GENOME} -rank 2 -default_version -agp_file #{SUPERCONTIG_AGP_FILE} > #{db_supercontig_level_log}")
		else
			system("echo 'Nothing to do' > #{db_supercontig_level_log}")
		end
		
	end

	file db_contig_level_log => db_chromosome_level_log do
		warn "---------------------------------------------------"
		warn "Loading contigs into database ensembl_#{GENOME}"
		warn "---------------------------------------------------"
		
		SUPERCONTIG_AGP_FILE ? seq_region_id = 3 : seq_region_id = 2 
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/load_seq_region.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -coord_system_name contig -coord_system_version v#{GENOME} -rank #{seq_region_id} -sequence_level -default_version -fasta_file #{CONTIGS} > #{db_contig_level_log} 2> /dev/null")
	end
	
	file db_chromosome2contig_log => db_contig_level_log do
		
		warn "-----------------------------------------------------------------"
		warn "Loading chromosome2contig mapping into database ensembl_#{GENOME}"
		warn "-----------------------------------------------------------------"
		
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/load_agp.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -assembled_name chromosome -component_name contig -agp_file  #{CONTIG_AGP_FILE} > #{db_chromosome2contig_log}")
	end
	
	file db_supercontig2contig_log  => db_chromosome2contig_log do
	
		if SUPERCONTIG_AGP_FILE
			
			warn "------------------------------------------------------------------"
			warn "Loading supercontig2contig mapping into database ensembl_#{GENOME}"
			warn "------------------------------------------------------------------"
			
			system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/load_agp.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -assembled_name supercontig -component_name contig -agp_file  #{SUPERCONTIG_AGP_FILE} > #{db_supercontig2contig_log}")
		
		else
			system("echo 'Nothing to do' > #{db_supercontig2contig_log}")
		end
	
	end
		
	file db_set_attributes_log => [ 'dump', db_supercontig2contig_log ] do
		warn "-------------------------------------------------------------"
		warn "Setting attributes for assembly in database ensembl_#{GENOME}"
		warn "-------------------------------------------------------------"
		
		system("perl #{ENSEMBL_PATH}/ensembl-production/scripts/production_database/populate_production_db_tables.pl -d ensembl_#{GENOME} -h #{MYSQL_HOST} -u #{MYSQL_USER} -p #{MYSQL_PW} -dumppath #{DUMP_PATH} -md #{PRODUCTION_DB} -mh #{MYSQL_HOST} -mu #{MYSQL_USER} -mp #{MYSQL_PW}  > #{db_set_attributes_log}")
	
	end
	
	file db_set_toplevel_log => db_set_attributes_log do
		warn "----------------------------------------------------------------------"
		warn "Setting toplevel attributes for assembly in database ensembl_#{GENOME}"
		warn "----------------------------------------------------------------------"
		
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/set_toplevel.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} > #{db_set_toplevel_log}")
		
	end
	
	file dumped_sequences_log => db_set_toplevel_log do
		warn "-------------------------------------------------------------"
		warn "Dumping toplevel sequences from assembly in ensembl_#{GENOME}"
		warn "-------------------------------------------------------------"
		
		system("perl #{ENSEMBL_PATH}/ensembl-analysis/scripts/sequence_dump.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -coord_system_name chromosome -output_dir dump > #{dumped_sequences_log}")
	
	end
	
	file dumped_database => dumped_sequences_log do 
		warn "--------------------------------------------------"
		warn "Dumping database ensembl_#{GENOME} (safe-keeping!)"
		warn "--------------------------------------------------"
	
		system("mysqldump --opt -u#{MYSQL_USER} -p#{MYSQL_PW} -h #{MYSQL_HOST} ensembl_#{GENOME} > #{dumped_database}")
	end
	
	file load_genewise_data => dumped_database do
	  warn "----------------------------------------------------"
	  warn " Loading dumped data into gene wise for replication"
	  warn "----------------------------------------------------"
	  
	  system("mysql -D genewise_#{GENOME} -u#{MYSQL_USER} -p#{MYSQL_PW} -h #{MYSQL_HOST} < #{dumped_database} > #{load_genewise_data}")
	end
	
	desc 'Loads assembly into Database'
	task :load => load_genewise_data

end


###############################
# Generate basic config files #
###############################

analysis_file = "configs/pipeline-confings/modules/Bio/EnsEMBL/Analysis/RunnableDB.pm"
pipeline_file = "configs/pipeline-confings/modules/Bio/EnsEMBL/Pipeline/Analysis.pm"

namespace :configs do
	
	file analysis_file => [ 'configs/pipeline-configs/modules/Bio/EnsEMBL/' , "assembly:load" ] do
		warn "---------------------------------------------------"
		warn "Copying the analysis config files from CVS checkout"
		warn "---------------------------------------------------"
		
		system("cp -R #{ENSEMBL_PATH}/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline configs/pipeline-configs/modules/Bio/EnsEMBL/")
	end
	
	file pipeline_file => analysis_file do
		warn "---------------------------------------------------"
		warn "Copying the pipeline config files from CVS checkout"
		warn "---------------------------------------------------"
	
		system("cp -R #{ENSEMBL_PATH}/ensembl-analysis/modules/Bio/EnsEMBL/Analysis configs/pipeline-configs/modules/Bio/EnsEMBL/")
	end
	
	task :create => pipeline_file do	
		warn "Exporting updated PERL5LIB"
		ENV['PERL5LIB'] = "#{ENV['PERL5LIB']}:/#{Dir.getwd}/configs/pipeline-configs/modules"
	end
	
end
	
#############################################
# INITIALIZING ANALYSES #####################
#############################################

config_file = "analyses_and_rules/contig_ana.conf"
analysis2db_log = "log/analysis.contig2db.#{GENOME}.log"

# Analyses define which software should be run
# and on which sequence type (Contig, Chromosome, etc)

desc 'Analyses to be performed on the assembly'
namespace :analysis do

	file config_file => [ 'analyses_and_rules' , "configs:create" ]  do
		
		f = File.new(config_file, "w+")
		
		f.puts analyses.join("\n")
	  
		f.close
		
	end
	
	file analysis2db_log => config_file do
		warn "-------------------------------------------------------------"
		warn "Adding analyses into DB ensembl_#{GENOME}"
		warn "-------------------------------------------------------------"
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/analysis_setup.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -read -file #{config_file} > #{analysis2db_log}")
	
	end

	task :define => analysis2db_log

end

ana_rule = "analyses_and_rules/rules.conf"
input_ids_log = "log/rules.make_input_ids.log"

###############################
# Generate rules for analyses #
###############################

# Rules define the hierarchy of analyses
# Names correspond to the analysis objects
namespace :rules do

	file ana_rule => "analysis:define" do
		warn "---------------------------"
		warn "Creating rules for analyses"
		warn "---------------------------"
		
		f = File.new(ana_rule,"w+")
		
		f.puts rules.join("\n")
		
		f.close
		
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/rule_setup.pl -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -read -file #{ana_rule}")
	
	end
	
	file input_ids_log => ana_rule do
		
		warn "---------------------------------"
		warn "Generating input ids for analysis"
		warn "---------------------------------"
		
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/make_input_ids -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -logic_name SubmitContig -slice -coord_system contig > #{input_ids_log}")
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/make_input_ids -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -logic_name SubmitChromosome -slice -coord_system chromosome >> #{input_ids_log}")
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/make_input_ids -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -logic_name SubmitGenome -slice -coord_system chromosome >> #{input_ids_log}")
		system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/make_input_ids -dbhost #{MYSQL_HOST} -dbuser #{MYSQL_USER} -dbport #{MYSQL_PORT} -dbname ensembl_#{GENOME} -dbpass #{MYSQL_PW} -logic_name SubmitSlice -slice -coord_system chromosome -slice_size 1000000 >> #{input_ids_log}")
	  
	end

	task :create => input_ids_log

end

analysis_database_pm_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Analysis/Config/Databases.pm"
update_analysis_database_pm_log = "log/submission.update_databasepm.log"
analysis_general_pm_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Analysis/Config/General.pm"
pipeline_batch_queue_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm"
pipeline_general_pm_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Pipeline/Config/General.pm"
analysis_pmatch_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/Pmatch.pm"
analysis_blast_config_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Analysis/Config/Blast.pm"
analysis_minigenewise_config_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/BlastMiniGenewise.pm"
analysis_killlist_config_file = "#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Analysis/Config/GeneBuild/KillListFilter.pm"

##################################
# PREPARE PROJECT FOR SUBMISSION #
##################################

namespace :submission do


	file analysis_database_pm_file => input_ids_log do
		warn "-------------------------------------------------"
		warn "Creating a copy of Databases.pm for this analysis"
		warn "-------------------------------------------------"
		
		system("cp #{ENSEMBL_PATH}/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Config/Databases.pm.example #{analysis_database_pm_file}")
	
	end
	
	# The Database.pm file requires information on our database, needs updating"
	file update_analysis_database_pm_log => analysis_database_pm_file do
		warn "--------------------------------------------"
		warn "Updating Databases.pm file for this analysis"
		warn "--------------------------------------------"
	
		ref_db = false
		gw_db = false
		
		file = File.open(analysis_database_pm_file)
		output = File.new(update_analysis_database_pm_log, "w+")
		while ( line = file.gets)
			
			ref_db = true if line.match(/^(\s+)REFERENCE_DB.*$/)
			ref_db = false if line.match(/^(\s+)\},$/)
			gw_db = true if line.match(/^(\s+)GENEWISE_DB.*$/)
			gw_db = false if line.match(/^(\s+)\},$/)
 					
			if line.include?("-dbname")
				output.puts line.gsub(/\'\'/, "\'ensembl_#{GENOME}\'") if ref_db
				output.puts line.gsub(/\'\'/, "\'genewise_#{GENOME}\'") if gw_db
			elsif line.include?("-host") 
				output.puts line.gsub(/\'\'/, "\'#{MYSQL_HOST}\'") if ref_db or gw_db
			elsif line.include?("-port")
				output.puts line.gsub(/\'\'/, "\'#{MYSQL_PORT}\'") if ref_db or gw_db
			elsif line.include?("-user")
				output.puts line.gsub(/\'\'/, "\'#{MYSQL_USER}\'") if ref_db or gw_db
			elsif line.include?("-pass")  
				output.puts line.gsub(/\'\'/, "\'#{MYSQL_PW}\'") if ref_db or gw_db
			else
				output.puts line
			end
			
		end
		output.close
		file.close
		
		system("cp #{update_analysis_database_pm_log} #{analysis_database_pm_file}")
	end
	
	file analysis_general_pm_file => update_analysis_database_pm_log do
		warn "------------------------------------------"
		warn "Updating General.pm file for this analysis"
		warn "------------------------------------------"
	
		system("cp #{analysis_general_pm_file}.example #{analysis_general_pm_file}")
		
	end
	
	file pipeline_batch_queue_file => analysis_general_pm_file do
		warn "---------------------------------------------"
		warn "Updating BatchQueue.pm file for this Pipeline"
		warn "---------------------------------------------"
		
		system("cp #{pipeline_batch_queue_file}.example #{pipeline_batch_queue_file}")
	end

	file pipeline_general_pm_file => pipeline_batch_queue_file do
		warn "---------------------------------------------"
		warn "Updating General.pm file for this Pipeline"
		warn "---------------------------------------------"
	
		system("cp #{pipeline_general_pm_file}.example #{pipeline_general_pm_file}")

	end
	
	file analysis_pmatch_file => pipeline_general_pm_file do
	  warn "---------------------------------------------"
		warn "Updating Pmatch.pm file for this analysis"
		warn "---------------------------------------------"
		
		system("cp #{analysis_pmatch_file}.example #{analysis_pmatch_file}")
		
	end
	
	file analysis_blast_config_file => analysis_pmatch_file do
	  warn "---------------------------------------------"
		warn "Updating Blast.pm file for this analysis"
		warn "---------------------------------------------" 
	  
	  system("cp #{analysis_blast_config_file}.example #{analysis_blast_config_file}")
	end
	
	file analysis_minigenewise_config_file => analysis_blast_config_file do
		warn "----------------------------------------------------"
		warn "Updating BlastMiniGenewise.pm file for this analysis"
		warn "----------------------------------------------------"
	
		system("cp #{analysis_minigenewise_config_file}.example #{analysis_minigenewise_config_file}")
	
	end
	
	file analysis_killlist_config_file => analysis_minigenewise_config_file do
		warn "----------------------------------------------------"
		warn "Updating KillListFilter.pm file for this analysis"
		warn "----------------------------------------------------"
	
		system("cp #{analysis_killlist_config_file}.example #{analysis_killlist_config_file}")
	end
	
	task :create => [ "output" , analysis_killlist_config_file ] 
		
end


task :synchronise_dbs => "submission:create" do

		warn "--------------------------------------"
		warn " Synchronizing databases"
		warn "--------------------------------------"
		
		tables = [ "analysis" , "assembly", "assembly_exception", "coord_system", "seq_region", "seq_region_attrib", "meta", "attrib_type", "seq_region_synonym" ]	
		
		tables.each do |table|
			system("mysqldump --no-create-db -u#{MYSQL_USER} -p#{MYSQL_PW} -h#{MYSQL_HOST} ensembl_#{GENOME} #{table} > dump/#{table}.sql")
			system("mysql -D genewise_#{GENOME} -u#{MYSQL_USER} -p#{MYSQL_PW} -h#{MYSQL_HOST} < dump/#{table}.sql")
		end
		
		warn "--------------------------------------"
		warn " Done synching tables "
		warn "--------------------------------------"

end



#############################
### STARTS THE PIPELINE #####
#############################

desc 'Start the pipeline'
task :build => "synchronise_dbs" do
  
	  warn "**********************************************************"
	  warn "Everything has been created, now begins the manual labor!"
	  warn "1. Edit BatchQueue.pm to include all necessary analyses and configure environment"
	  warn "2. Edit Pmatch.pm to add protein file and database info"
	  warn "3. Rinse an repeat for all other analyses (need to add to script!)"
	  warn "Run the test_RunnableDB script to check data"
	  warn "Run the pipeline!"
end



task :pathfile do
  path = "PERL5LIB=$PERL5LIB:#{Dir.getwd}/configs/pipeline-configs/modules:#{ENSEMBL_PATH}/ensembl-pipeline/scripts:#{ENSEMBL_PATH}/ensembl-killlist/modules:/usr/bin/tRNAscan-SE"
  f = File.new("update_path.sh","w+")
  f.puts path
  f.close
end


#############################
# Database queries ##########
#############################

task :analysis_ids do
  EnsemblDB::DBConnection.connect
  analyses_input_ids = EnsemblDB::InputIdAnalysis.find(:all).group_by{|i| i.input_id_type }
  analyses_input_ids.each do |typ,input_ids|
    puts "#{typ}"
    input_ids[0..12].each {|i| puts "\t#{i.input_id}"}
  end
end

task :sanity_check do
  system("perl #{ENSEMBL_PATH}/ensembl-pipeline/scripts/pipeline_sanity.pl -dbname ensembl_#{GENOME} -dbuser #{MYSQL_USER} -dbpass #{MYSQL_PW} -dbhost #{MYSQL_HOST}")
end

task :perl_lib do
  # Create a PERL5LIB source
  f = File.new("perl5lib.sh","w+")
  f.puts "PERL5LIB=#{ENSEMBL_PATH}/bioperl-live:#{ENSEMBL_PATH}/ensembl-pipeline/modules:#{ENSEMBL_PATH}/ensembl-pipeline/scripts:#{ENSEMBL_PATH}/ensembl-analysis/modules:#{ENSEMBL_PATH}/ensembl/modules:#{Dir.getwd}/configs/pipeline-configs/modules:/usr/bin/tRNAscan-SE:#{Dir.getwd}/configs/pipeline-configs/modules/Bio/EnsEMBL/Analysis/Runnable:"
  f.puts "export PERL5LIB"
  f.puts "echo $PERL5LIB"
  f.puts "BLASTMAT=/references/software/util/BLOSSUM62"
  f.puts "export BLASTMAT"
  f.close
  
end

#############################
# Cleaning up and resetting #
#############################

desc "Resetting the entire analysis"
task :clean do
	system("rm -R analyses_and_rules")
	system("rm -R log")
	system("rm -R dump")
	system("rm -R configs")
	system("rm -R output*")
end
