#!/usr/bin/ruby
# == NAME
# pasa_create_chunks.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -g | --genome ] |[ -c | --cufflinks ] | [ -t | --trinity ] | [ -s | --chunks ] | [ -r, --ref] | [-a , --analysis ]
# == DESCRIPTION
# Creates several PASA chunks from input data to increase parallelization
#
# == OPTIONS
#  -h,--help::                  Show help
#  -g,--genome=GENOME::         Genome file
#  -c,--cufflinks=CUFFLINKS::   cufflinks file
#  -t,--trinity=TRINITY::       Trinity file
#  -s,--chunks=CHUNKS::         number of chunks to create
#  -r,--ref=REF::               Reference annotation to polish
#  -a,--analysis=ANALYSIS::    A name for this analysis (to later identify the SQL databases)

# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, mphoeppner@gmail.com

require 'rdoc/usage'
require 'optparse'
require 'ostruct'
require 'logger'
require 'bio'
require 'fileutils'
require 'pathname'

### Define modules and classes here

def make_chunk_store(fasta_bin,chunks)
  chunk_store = {}
  chunk_counter = 0
  
  chunks.times do 
      chunk_counter += 1
      chunk_store[chunk_counter] = []
  end
  
  fasta_bin = fasta_bin.sort_by{|k,v| v }
  while fasta_bin.length > 0
    chunk_counter = 0
    chunks.to_i.times do
      next if fasta_bin.empty?
      chunk_counter += 1
      entry = fasta_bin.shift
      chunk_store[chunk_counter] << entry[0]
    end
  end
  return chunk_store
end

def parse_lines(file,sequences)

  answer = []
  infile = File.open(file, "r")
  while (line = infile.gets)
	line.strip!
	e = line.split("\t")[0]
	next unless sequences.include?(e)
	answer << line
  end
  infile.close
  return answer
  
end

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-h","--help","Display the usage information") {RDoc::usage}
opts.on("-g","--genome", "=GENOME","Genome sequence") {|argument| options.genome = argument }
opts.on("-c","--cufflinks", "=CUFFLINKS","Cufflinks GTF file") {|argument| options.cufflinks = argument }
opts.on("-t","--trinity", "=TRINITY","Trinity fasta file") {|argument| options.trinity = argument }
opts.on("-s","--chunks", "=CHUNKS","Number of chunks") {|argument| options.chunks = argument }
opts.on("-r","--ref", "=REF","Reference annotation") {|argument| options.ref = argument }
opts.on("-a","--analysis", "=ANALYSIS","Name for this analysis") {|argument| options.analysis = argument }

opts.parse! rescue RDoc::usage('usage')

options.chunks = 2 unless options.chunks
raise "No Genome sequence provided!" unless options.genome
raise "No Trinity assembly provided!" unless options.trinity
#raise "No reference annotation provided!" unless options.ref
#raise "No cufflinks annotation provuded!" unless options.cufflinks
raise "No analysis name defined!" unless options.analysis
pasa_dir = "/sw/bioinfo/pasa2/r20140414"
pasa_conf_template = "#{pasa_dir}/pasa_conf/pasa.alignAssembly.Template.txt"
wd = Dir.getwd
analysis = "crow27"

genome_file = options.genome
cufflinks_file = options.cufflinks
trinity_file = options.trinity

options.cufflinks ? cufflinks_flag = "--cufflinks_gtf cufflinks.gtf" : cufflinks_flag = ""
models_flag = ""

#Container for fasta
fasta_bin = {}

# Gather stats on fasta sequence lengths
Bio::FastaFormat.open(options.genome).each_entry do |entry|
  fasta_bin[entry.definition] = entry.naseq.length
end

chunk_store = make_chunk_store(fasta_bin,options.chunks.to_i)

chunk_store.each do |chunk,sequences|
  
  # Create the folder for this chunk
  system("mkdir -p chunk_#{chunk}")
  
  # Create a Pasa config file
  system("sed 's/\<__MYSQLDB__\>/pasa2_#{options.analysis}_polish_chunk_#{chunk}/g' #{pasa_conf_template} > chunk_#{chunk}/pasa.conf")
  
  # Grep out the relevant lines from the cufflinks assembly and reference annotation

  if options.cufflinks  
	  cufflinks_lines = parse_lines(options.cufflinks,sequences)
	  cufflinks_file = File.new("#{wd}/chunk_#{chunk}/cufflinks.gtf","w+")
	  cufflinks_lines.each{|l| cufflinks_file.puts l}
	  cufflinks_file.close
  end

  if options.ref
        models_flag = "-A -L --annots_gff3 models.gff"
  	model_lines = parse_lines(options.ref,sequences)
  	model_file = File.new("#{wd}/chunk_#{chunk}/models.gff", "w+")
  	model_lines.each{|l| model_file.puts l }
  	model_file.close
  end


  # Copy over all Trinity files  
  system("cp #{options.trinity}* chunk_#{chunk}/")
  
  outfile = File.new("#{wd}/chunk_#{chunk}/genome_chunk.fa", "w+")
  
  # Extract the contigs for this chunk
  Bio::FastaFormat.open(options.genome).each_entry do |entry|
    next unless sequences.include?(entry.definition)
	outfile.puts entry.to_s
  end
  
  outfile.close
  
  # Write a bsub file  
  bsub_file = File.new("#{wd}/chunk_#{chunk}/bsub.pasa_chunk", "w+")
  
  command_lines = [
    "#BSUB -n 16",
    "#BSUB -R \"span[hosts=1]\"",
	"#BSUB -e err.pasa",
	"#BSUB -J PasaChunk-#{chunk}",
	"Launch_PASA_pipeline.pl -c pasa.conf -C -R --ALIGNER blat,gmap -g genome_chunk.fa -t #{options.trinity}.clean -T -u #{options.trinity} --CPU 16 --TRANSDECODER --transcribed_is_aligned_orient #{cufflinks_flag} #{models_flag}"
  ]
  
  command_lines.each do |cl|
    bsub_file.puts cl
  end
  
  bsub_file.close
  
end

