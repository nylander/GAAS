#!/usr/bin/ruby
# == NAME
# add_track.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -g | --gff ] |[ -s | --species ] | 
# == DESCRIPTION
# A script to add a track to an existing WebApollo installation. 
#
# == OPTIONS
#  -h,--help::                  Show help
#  -s,--species=SPECIES::       Name of the species
#  -g,--gff=GFF::       	Annotation file to process (optional)
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
opts.on("-s","--species", "=SPECIES","Name of the species") {|argument| options.species = argument }
opts.on("-c","--category", "=CATEGORY","Track category") {|argument| options.category = argument }
opts.on("-l","--label", "=LABEL","Label") {|argument| options.label = argument }
opts.on("-o","--outfile", "=OUTFILE","Output") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") { 
	puts opts 
	exit 
}

opts.parse! 

### Usernames, passwords and locations

home = ENV['HOME']

config = {
	:web_apollo_build => "/big/webapollo/build/#{options.species}", 	# The location where this WA project is to be build
	:web_apollo_data => "/big/data/#{options.species}/data"	# Location of the data store for this species
}

### The workflow

raise "This species (#{options.species}) does not have a build directory under #{config[:web_apollo_build]})" unless File.directory?(config[:web_apollo_build])
raise "This species (#{options.species}) does not have a build directory under #{config[:web_apollo_data]})" unless File.directory?(config[:web_apollo_data])

tracks = { "protein" => "\"match_part\": \"blue-80pct\"" , 
  "est" => "\"match_part\": \"green-80pct\"" , 
  "tRNA" => "\"exon\": \"green-80pct\"",
  "ncRNA" => "\"exon\": \"green-80pct\"",
  "synteny" => "\"nucleotide_motif\": \"springgreen-80pct\"", 
  "rnaseq" => "\"match_part\": \"orange-80pct\"" , 
  "abinito" => "\"match_part\": \"springgreen-80pct\"",
  "gene" => "\"wholeCDS\": null, \"CDS\": \"blueviolet-80pct\", \"UTR\": \"darkorange-60pct\", \"exon\" : \"container-100pct\"",
  "lift-over" => "\"wholeCDS\": null, \"CDS\": \"green-80pct\", \"UTR\": \"darkorange-60pct\", \"exon\" : \"container-100pct\""
  }

abort "No track provided" if options.infile.nil?
abort "No species provided" unless options.species
abort "Must provide a category for this track" unless options.category
abort "Must provide a label for this track" unless options.label

## Some useful variables

track = nil
track_number = {}
proceed = false
failure = 0

##

# List the current information
puts "Your data:"
puts "----------"
puts "Track to  load: #{options.infile}"
puts "WebApollo portal: #{config[:web_apollo_build]}"
puts "Track label: #{options.label}"
puts "Category: #{options.category}"
puts "Track type selection..."
puts

# Ask about the type of data


while proceed == false
  
  abort "Can't make up your mind, eh? Aborting..." if failure > 3

  puts "###############################################"  
  puts "Choose the track type from the following list:"

  number = 0	

  tracks.sort_by{|t,f| t}.each do |track_type,formatting|
    number += 1 
    track_number[number] = track_type
    track = track_type
    puts "\t#{number}\t#{track_type}"
  end

  puts "##############################################"
  puts "Enter the number corresponding to your choice:"

  selection = gets.chomp

  next unless track_number.has_key?(selection.to_i)

  puts "You selected: #{track_number[selection.to_i]} - is that correct? (Y/N)"
  
  answer = gets.chomp

  if answer.downcase == "y"
	proceed = true
	track = track_number[selection.to_i]
  else
	failure += 1
  end
  
  
end

puts "Loading track as type #{track}, with the label #{options.label} into the category #{options.category} of webapollo directory #{config[:web_apollo_build]}"

puts "Proceed? (Y/N)"

answer = gets.chomp

abort "Aborting..." unless answer.downcase == "y"

if track == "gene" or track == "lift-over"
	system "perl #{config[:web_apollo_build]}/bin/flatfile-to-json.pl --gff #{options.infile} --out #{config[:web_apollo_data]} --arrowheadClass trellis-arrowhead --getSubfeatures --subfeatureClasses '{#{tracks[track]}}'  --cssClass container-16px --config '{ \"category\": \"#{options.category}\" }' --type mRNA --trackLabel #{options.label}"
elsif track == "tRNA"
        system "perl #{config[:web_apollo_build]}/bin/flatfile-to-json.pl --gff #{options.infile} --out #{config[:web_apollo_data]} --arrowheadClass trellis-arrowhead --getSubfeatures --subfeatureClasses '{#{tracks[track]}}'  --cssClass container-16px --config '{ \"category\": \"#{options.category}\" }' --type tRNA --trackLabel #{options.label}"
elsif track == "ncRNA"
	system "perl #{config[:web_apollo_build]}/bin/flatfile-to-json.pl --gff #{options.infile} --out #{config[:web_apollo_data]} --arrowheadClass trellis-arrowhead --getSubfeatures --subfeatureClasses '{#{tracks[track]}}'  --cssClass container-16px --config '{ \"category\": \"#{options.category}\" }' --type ncRNA --trackLabel #{options.label}"
else
	system "perl #{config[:web_apollo_build]}/bin/flatfile-to-json.pl --gff #{options.infile} --out #{config[:web_apollo_data]} --arrowheadClass trellis-arrowhead --getSubfeatures --subfeatureClasses '{#{tracks[track]}}' --cssClass container-16px --config '{ \"category\": \"#{options.category}\" }' --trackLabel #{options.label}"	
end
