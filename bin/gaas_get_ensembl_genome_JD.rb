#!/usr/bin/ruby
# A script to download the latest version of an EnsEMBL genome,
# including GTF files.
# = USAGE
# ./get_ensembl_genome.rb gorilla_gorilla 67

def command?(command)
       system("which #{ command} > /dev/null 2>&1")
end

require 'fileutils'
require 'net/ftp'

species = ARGV.shift
version = ARGV.shift
pathVersion=""
nbis_path="/projects/references/genomes"

abort "You should really remember to provide a species_name and optionaly a release!" if species.nil?

raise "This does not look like a valid species name (must be lower case scientific name with underscore)" unless species.match(/[a-z]*_[a-z]*/)
if version.nil? 
	puts "No EnsEMBL release provided, we will take the current release as default"
	pathVersion ="current"
else
	pathVersion="release-"+version
end

$ftp=Net::FTP.new
#connect to Ensembl vertebrates ftp
$ftp.connect("ftp.ensembl.org")
$ftp.passive = true
$ftp.login

# Check for fastx toolkit
fasta_loaded = false
cufflinks_loaded = false

unless command?("fasta_formatter")
	puts "Couldn't find fasta_formatter in your PATH"
#	command("module load fastx/0.0.13")
#	fasta_loaded = true
	raise "You must load the fastx module !"
end

unless command?("gffread -x")
	puts "Couldn't find cufflinks in your PATH, trying to load it"
#	system("module load cufflinks/2.2.1")
#	cufflinks_loaded = true
	raise "You must load the cufflinks module !"
end

found="no"
fetched_path=""
puts "Downloading genome for #{species}, release #{pathVersion} on ensembl FTP server"
#dirExists = ftp.ChangeRemoteDir("/pub/#{pathVersion}/fasta/#{species}/dna/")
begin
	if pathVersion == "current"
		fetched_path="/pub/#{pathVersion}_fasta/#{species}/dna/"
		$ftp.chdir("/pub/#{pathVersion}_fasta/#{species}/dna/")
	else
		fetched_path="/pub/#{pathVersion}/fasta/#{species}/dna/"
		$ftp.chdir("/pub/#{pathVersion}/fasta/#{species}/dna/")
	end
rescue
	puts "Not found on ensembl FTP server !!"
else
	found="yes"
end

if found=="no"
	puts "Downloading genome for #{species}, release #{pathVersion} on ensemblgenomes FTP server"
	$ftp.close
	#connect to the other Ensembl ftp server
	$ftp.connect("ftp.ensemblgenomes.org")
	$ftp.passive = true
	$ftp.login
	begin
		puts "Try on Plants..."
		fetched_path="/pub/plants/#{pathVersion}/fasta/#{species}/dna/"
		$ftp.chdir("/pub/plants/#{pathVersion}/fasta/#{species}/dna/")
	rescue
		puts "Not found on plants !!"
	else
		found="yes"
	end
	if found=="no"
		begin
			puts "Try on metazoa..."
			fetched_path="/pub/metazoa/#{pathVersion}/fasta/#{species}/dna/"
			$ftp.chdir("/pub/metazoa/#{pathVersion}/fasta/#{species}/dna/")
		rescue
			puts "Not found on metazoa !!"
		else
			found="yes"
		end
	end
	if found=="no"
		begin
			puts "Try on fungi..."
			fetched_path="/pub/fungi/#{pathVersion}/fasta/#{species}/dna/"
			$ftp.chdir("/pub/fungi/#{pathVersion}/fasta/#{species}/dna/")
		rescue
			puts "Not found on fungi !!"
		else
			found="yes"
		end
	end
	if found=="no"
		begin
			puts "Try on protists..."
			fetched_path="/pub/protists/#{pathVersion}/fasta/#{species}/dna/"
			$ftp.chdir("/pub/protists/#{pathVersion}/fasta/#{species}/dna/")
		rescue
			puts "Not found on protists !!"
		else
			found="yes"
		end
	end
	if found=="no"
		begin
			puts "Try on bacteria..."
			fetched_path="/pub/bacteria/#{pathVersion}/fasta/#{species}/dna/"
			$ftp.chdir("/pub/bacteria/#{pathVersion}/fasta/#{species}/dna/")
		rescue
			puts "Not found on bacteria !!"
		else
			found="yes"
		end
	end
end

if found == "no"
	raise "Something went wrong. We check all Ensembl database and this species name was not found (Check spelling!)"
end


releaseName=String($ftp.pwd).split('/')[2]
pathToWrite=""
if File.directory?(nbis_path)
	pathToWrite="#{nbis_path}/#{species}/EnsEMBL/#{releaseName}";

else
	puts "apparently we are not on the nbis annotation cluster"
        pathToWrite="#{species}_#{releaseName}"
end

# check presence of directories
if File.directory?(pathToWrite)
        raise "The directory #{pathToWrite} exists ! Apparently you already download that genome release."
end

puts "I will save result here: #{pathToWrite}"

	FileUtils.mkdir_p("#{pathToWrite}")
	FileUtils.cd("#{pathToWrite}")
	file = $ftp.nlst.find{|f| f.include?("dna.primary_assembly.fa.gz") }
	
	file = $ftp.nlst.find{|f| f.include?("dna.toplevel.fa.gz")} if file.nil? or file.length < 1
	
	puts "Downloading genome sequence: #{file}"
	$ftp.get(file,File.basename(file)) # Download the file
	puts "Extracting genomes"
	#system("rm *.fa")unless Dir.entries(Dir.getwd).select{|e| e.include?(".fa")}.empty? #remove .fa file if exists
	system("gunzip #{file}")	
	file = Dir.entries(Dir.getwd).find{|f| f.include?(".fa")}	
	puts "Cleaning genome #{file}"
	system("sed 's/ dna.*//g' #{file} > #{species}.fa.tmp") # Remove the rubbish ensembl encodes in the FASTA headers
	system("rm #{file}")
	puts "Formatting genome"
	system("/sw/bioinfo/fastx-0.0.13/fasta_formatter -i #{species}.fa.tmp -o tmp.fa -w 80") # Ensure that FASTA lines are of equal length
	system("rm *fa.tmp")
	system("mv tmp.fa #{species}.fa")
	system("rm *toplevel*")
	system("rm *.tmp")

begin
	puts "Downloading GTF for #{species}"
	#prepare GTF path
	gtfPath=fetched_path.gsub(/\/fasta\//,"/gtf/")
	gtfPath=gtfPath.gsub(/\/dna\//,"/")
	gtfPath=gtfPath.gsub(/_fasta\//,"_gtf/")
	puts "FTP gtf path ckecked: #{gtfPath}"
	$ftp.chdir("#{gtfPath}")
	file = $ftp.nlst.find{|f| f.include?("gtf.gz")}
	raise "No GTF file was found for this species" if file.nil?
	$ftp.get(file,File.basename(file))
	system("gunzip *.gz")
	system("mv *.gtf #{species}.gtf")
	warn "Reformatting GTF file"
	system("gffread /projects/references/genomes/#{species}/EnsEMBL/#{releaseName}/#{species}.gtf -L -F -g /projects/references/genomes/#{species}/EnsEMBL/#{releaseName}/#{species}.fa -T -o tmp.gtf") # Make EnsEMBL gtf fully compatible with Cufflinks and Tophat
	system("mv tmp.gtf #{species}.gtf")
	system("sed 's/gene_name \"[A-Za-z0-9_]*\"\; //g' #{species}.gtf > #{species}.no_names.gtf") # Create a copy without gene_name attribute (needed for cufflinks package)
rescue
	raise "Something went wrong. Maybe no GTF file exists for that species."
end

$ftp.close

puts "All Done!"


# Unload modules only if we had to specifically load them for this script

#system("module unload cufflinks/2.1.1") if cufflinks_loaded
#system("module unload fastx/0.0.13") if fasta_loaded
