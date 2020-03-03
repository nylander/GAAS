#!/usr/bin/ruby
# == NAME
# create_annotation_project
#
# == USAGE
#  ./create_annotation_project [ -h | --help ]
#                    [ -s | --species ] |[ -g | --genome ] |  [ -v | --version ] | [ -f, --force ]
# == DESCRIPTION
# 
#
# == OPTIONS
#  -h,--help::                  Show help
#  -s,--SPECIES=SPECIES::		Species name in snake_case (e.g. homo_sapiens) [Required]
#  -g,--genome=GENOME::			Genome fasta file [Required]
#  -v,--version=VERSION:: 		Assembly version (if not specified, 1.0 is assumed)
#  -f,--force=FORCE::			Force the creation of folder structure if it already exists (everything will be wiped, careful!)
#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, mphoeppner@gmail.com

require 'rdoc/usage'
require 'optparse'
require 'ostruct'
require 'logger'
require 'fileutils'

### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-h","--help","Display the usage information") {RDoc::usage}
opts.on("-s","--species", "=SPECIES","Species") {|argument| options.species = argument }
opts.on("-g","--genome", "=GENOME","Genome sequence") {|argument| options.genome = argument }
opts.on("-v","--version", "=VERSION","Version of genome assembly") {|argument| options.version = argument }
opts.on("-f","--force", "=FORCE","Force project creation if exists") {|argument| options.force = true }

opts.parse! rescue RDoc::usage('usage')

# Check the sanity of the settings

abort "Must specify a species name in snake_case (--species)" unless options.species or options.species.match(/[a-z]*_[a-z]*/) == false
abort "Must specify a genome sequence (--genome)" unless options.genome

unless options.version
	warn "No assembly version specified (--version), assuming 1.0"
	options.version = "1.0"
end

# Set useful, derived variables

WORK_DIR = Dir.getwd
FASTA_EXPLODE = "/sw/bioinfo/exonerate/fastaexplode"
FASTA_FORMATTER = "/sw/bioinfo/fastx-0.0.13/fasta_formatter"

# Stupid hack to check whether genome sequence is here or a path on the system...
options.genome.include?("/") ? genome_location = options.genome : genome_location = "#{WORK_DIR}/#{options.genome}"
genome_name = genome_location.split("/")[-1]

ANNOTATION_ROOT = "/projects/annotation"

PROJECT_PATH = "#{ANNOTATION_ROOT}/#{options.species}"

# Define the folder structure
# Indentation is meaningless, the order matters

folders = [ 
	"#{PROJECT_PATH}/RNAseq" ,
	"#{PROJECT_PATH}/trinity" ,
		"#{PROJECT_PATH}/RNAseq/trimmed" ,
			"#{PROJECT_PATH}/RNAseq/trimmed/normalized" ,
				"#{PROJECT_PATH}/RNAseq/trimmed/normalized/C100/" ,
	"#{PROJECT_PATH}/v#{options.version}" ,
		"#{PROJECT_PATH}/v#{options.version}/cegma" ,
		"#{PROJECT_PATH}/v#{options.version}/genome" ,
		"#{PROJECT_PATH}/v#{options.version}/contigs" ,
		"#{PROJECT_PATH}/ab-initio" ,
			"#{PROJECT_PATH}/ab-initio/snap" ,
			"#{PROJECT_PATH}/ab-initio/augustus" ,
		"#{PROJECT_PATH}/v#{options.version}/maker" ,
		"#{PROJECT_PATH}/v#{options.version}/cufflinks" ,
		"#{PROJECT_PATH}/v#{options.version}/tophat" ,
		"#{PROJECT_PATH}/v#{options.version}/trinity" ,
		"#{PROJECT_PATH}/v#{options.version}/pasa" ,
		"#{PROJECT_PATH}/v#{options.version}/trnascan" ,
		"#{PROJECT_PATH}/v#{options.version}/rfam" ,
	"#{PROJECT_PATH}/refseqs"	
]

# Create project path (or not...)

if File.directory?(PROJECT_PATH) and options.force.nil?
	abort "Can't create project folder - already exists and --force not specified."
elsif File.directory?(PROJECT_PATH) and options.force
	warn "Forced to create project path under #{PROJECT_PATH}, everything will be wiped!"
	#FileUtils.rm_r(PROJECT_PATH)
	#FileUtils.mkdir_p(PROJECT_PATH)
	abort "This function is currently not available. Please remove manually and try again..."
else
	warn "Creating annotation project under #{PROJECT_PATH}"
	FileUtils.mkdir_p(PROJECT_PATH)
end

FileUtils.cd(PROJECT_PATH) 

folders.each do |folder|
    warn "Creating: #{folder}"
	FileUtils.mkdir_p(folder)
	system("touch #{folder}/00README")	
end

FileUtils.cp genome_location , "#{PROJECT_PATH}/v#{options.version}/genome/#{genome_name}"

# Clean up the genome sequence

system("#{FASTA_FORMATTER} -i #{PROJECT_PATH}/v#{options.version}/genome/#{genome_name} -w 80 -o #{PROJECT_PATH}/v#{options.version}/genome/genome.reformatted.fa") 
system("#{FASTA_EXPLODE} -f #{PROJECT_PATH}/v#{options.version}/genome/genome.reformatted.fa -d #{PROJECT_PATH}/v#{options.version}/contigs/")


