#!/usr/bin/ruby
# == NAME
# gff2hints.rb
#
# == USAGE
#  ./this_script.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# Converts Maker-generated evidence alignments from GFF3 format to Augustus-compatbile hint format
#
# == OPTIONS
#  -h,--help::                  Show help
#  -i,--infile=INFILE::         input file
#  -o,--outfile=OUTFILE::       output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'
require 'bio'

### Define modules and classes here


def cufflinks2hints(entries)
  
  bin = entries.records.group_by{|r| r.attributes["transcript_id"].gsub(/\"/, '')}
  
  bin.each do |transcript_id,records|
    
    strand = records[0].strand
    
    strand == "+" ? first_exon = records[0] : first_exon = records[-1]
    strand == "+" ? last_exon = records[-1] : last_exon = records[1]
    
    # Building exon part hints
    records.each do |record|
      
      puts "#{record.seqname}\tc2h\tep\t#{record.start}\t#{record.end}\t.\t#{record.strand}\t.\tgrp=#{transcript_id};pri=4,src=E"
      
      # building exon hints
      unless record == first_exon or record == last_exon
        puts "#{record.seqname}\tc2h\texon\t#{record.start}\t#{record.end}\t.\t#{record.strand}\t.\tgrp=#{transcript_id};pri=4,src=E"
      end
      
      # building TSS hints
      if record == first_exon and records.length > 2
          if strand == "+"
            puts "#{record.seqname}\tc2h\ttss\t#{record.start.to_i-20}\t#{record.start.to_i+20}\t.\t#{record.strand}\t.\tgrp=#{transcript_id};pri=4,src=E"
          else  
            puts "#{record.seqname}\tc2h\ttss\t#{record.end.to_i-20}\t#{record.end.to_i+20}\t.\t#{record.strand}\t.\tgrp=#{transcript_id};pri=4,src=E"
          end
      end
      
    end
    
    # Building the intron hints
    while records.length > 1
      upstream,downstream = records[0..1]
      puts "#{upstream.seqname}\tc2h\tintron\t#{upstream.end.to_i+1}\t#{downstream.start.to_i-1}\t.\t#{upstream.strand}\t.\tgrp=#{transcript_id};pri=4,src=E"     
      records.shift
    end
    
    
  end

  
end

def protein2hints(entries)
    
  bin = entries.records.select{|r| r.feature == "match_part" }.group_by{|r| r.attributes.find{|k,v| k == "Parent" } }

  bin.each do |protein_id,records|
    
    strand = records[0].strand
    
    strand == "+" ? first_exon = records[0] : first_exon = records[-1]
    strand == "+" ? last_exon = records[-1] : last_exon = records[1]
    
    records.each do |record|
      puts "#{record.seqname}\tp2h\tCDSpart\t#{record.start}\t#{record.end}\t.\t#{record.strand}\t.\tgrp=#{protein_id};pri=2,src=P"
     
  
      if record == first_exon
        
        target = record.attributes.find{|k,v| k == "Target" }[1]
                
        # The protein aligns from position one, which we assume is the start codon
        if target.start.to_i == 1
          
          if strand == "+"
            puts "#{record.seqname}\tp2h\tstart\t#{record.start.to_i-20}\t#{record.start.to_i+20}\t.\t#{record.strand}\t.\tgrp=#{protein_id};pri=2,src=P"
          else
            puts "#{record.seqname}\tp2h\tstart\t#{record.end.to_i-20}\t#{record.end.to_i+20}\t.\t#{record.strand}\t.\tgrp=#{protein_id};pri=2,src=P"
          end
          
        end
        
      end
    
    end
    
  end
  
end

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-i","--infile", "=INFILE","Input") {|argument| options.infile = argument }
opts.on("-o","--outfile", "=OUTFILE","Output") {|argument| options.outfile = argument }
opts.on("-s","--source", "=SOURCE","Source of the data") {|argument| options.source = argument }
opts.on("-h","--help","Display the usage information") {
  puts opts
  exit
}

opts.parse! 

sources = [ "cufflinks" , "protein2genome" ]

abort "Must specify a source (#{sources.join(',')})" unless options.source

options.infile ? input_stream = File.open(options.infile) : input_stream = $stdin
options.outfile ? output_stream = File.new(options.outfile,'w') : output_stream = $stdout

if options.infile.include?(".gff") or options.infile.include?(".gtf")
  

  if options.source == "cufflinks"
    
    entries = Bio::GFF.new(input_stream)   
    cufflinks2hints(entries)
    
  elsif options.source == "protein2genome"
    
    entries = Bio::GFF::GFF3.new(input_stream)  
    protein2hints(entries)
    
  end

elsif options.infile.include?(".bed")
  
  
  
end



output_stream.close



