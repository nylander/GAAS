#!/usr/bin/env perl

use strict;
use warnings;
use POSIX;
use File::Basename;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use GAAS::GAAS;

my $header = get_gaas_header();
my $start_run = time();

my $opt_fastafile;
my $opt_output="split_result";
my $opt_help = undef;
my $opt_nb_chunk;
my $opt_leftover="attach";
my $opt_nb_seq_by_chunk;
my $opt_overlap;
my $verbose=0;
my $opt_size_seq;


# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|fa|fasta=s' => \$opt_fastafile,
                  'nb_chunks=i' => \$opt_nb_chunk,
									'nb_seq_by_chunk=i'=> \$opt_nb_seq_by_chunk,
									'size_seq=i'=> \$opt_size_seq,
									'overlap=i'=> \$opt_overlap,
									'leftover=s'=> \$opt_leftover,
									'verbose|v=i'=> \$verbose,
                  'o|output=s'      => \$opt_output,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => "$header\nFailed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

# ---- check options ----
if (! ( $opt_fastafile and ($opt_nb_chunk or $opt_nb_seq_by_chunk)) ) {
    pod2usage( {
           -message => "\nAt least 2 parameters are mandatory:\nInput reference fasta file (-f)".
					 "\nnumber of chuncks (--nb_chunks) or/and number of sequences by chunk (--nb_seq_by_chunk)\n\n".
           "Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}

if ($opt_overlap){
	if ( ! $opt_size_seq){
		print "--overlap option requires the use of --size_seq option !\n";
		exit;
	}
	elsif ($opt_overlap >=  $opt_size_seq){
		print "--overlap option value must be smaller than the --size_seq option value!\n";
		exit;
	}
}

$opt_leftover = lc($opt_leftover);
if (lc($opt_leftover) ne "attach" and lc($opt_leftover) ne "detach" and lc($opt_leftover) ne "remove"){
	print "$opt_leftover is not an accepted value for --leftover parameter.".
	" Accepted value are: attach, detach or remove. Please read the help for more details\n";
	exit;
}

# information
if( $opt_size_seq ){
	print "=> <size_seq> option activated, we cut sequences to get $opt_size_seq nucleotides per sequence\n";
	if ($opt_overlap){
		print "=> <overlap> option activated, we add an overlap of $opt_overlap bp per sequence\n";
	}
	if($opt_leftover eq "attach"){
		print "=> <leftover> option is $opt_leftover. For each sequence the size of the last chunk can be between ".
		($opt_size_seq/2)."<=X<=$opt_size_seq or $opt_size_seq<X<".($opt_size_seq+($opt_size_seq/2))." \n";
	}
	elsif($opt_leftover eq "detach"){
		print "=> <leftover> option is $opt_leftover. For each sequence the size of the last chunk can be <= $opt_size_seq\n";
	}
	else{
		print "=> <leftover> option is $opt_leftover\n";
	}
}

if(! -f $opt_fastafile){
	print "$opt_fastafile is not a file\n."; exit;
}
# ----- output -------

my ($inputname,$path,$ext) = fileparse($opt_fastafile,qr/\.[^.]*/);
if (-d $opt_output) {
	print "$opt_output output directory already exits !\n";exit;
}
mkdir $opt_output;
my $out_prefix = $opt_output."/".$inputname;

##### MAIN ####

my $cpt_seq_out=0;

######### read fasta file #############
my $nb_seq=0;
my $nb_seq_original=0;
my $nb_nucl=0;
my $fasta1  = Bio::SeqIO->new(-file => $opt_fastafile , -format => 'Fasta');
while ( my $seqObj = $fasta1->next_seq() ) {

	# ----------- $opt_size_seq --------
	if( $opt_size_seq ){

		my $nb_seq_here = ceil($seqObj->length() / $opt_size_seq);
		my $total_size_seq = $seqObj->length();

		if ($opt_overlap){

			my $total_overlap = ($nb_seq_here-1) * ($opt_overlap*2); # total bp size of overlap
			my $total_size_with_overlap = $seqObj->length() + $total_overlap;
			$total_size_seq = $total_size_with_overlap;
			$nb_seq_here = ceil($total_size_with_overlap / $opt_size_seq);

		}

		my $leftover = $total_size_seq % $opt_size_seq;
		print "leftover $leftover\n" if ($verbose > 1);

		if($opt_leftover eq "attach"){

			if ($leftover < ceil($opt_size_seq/2) ){
				print "leftover info: last piece shorter than half of split_size ($leftover < ($opt_size_seq / 2) ).".
				"\nWe attach it to the previous chunk that will make one chunck less.\n" if ($verbose > 1);

				$nb_seq_here = ($nb_seq_here - 1 > 0) ? ($nb_seq_here-1) : 1;
				my $last_seq_size = ($seqObj->length() > $opt_size_seq) ? ($leftover+$opt_size_seq) : $seqObj->length();

				my $message = "Sequence ".$seqObj->id()." is length ".$seqObj->length().". I will create:";
				$message .= "\n  * ".($nb_seq_here-1)." sequences of $opt_size_seq nucleotides" if ($nb_seq_here-1 > 0);
				$message .= "\n  * 1 sequence of ".$last_seq_size." nucleotides\n";
				print $message if ($verbose > 0);
			}
			else{
				my $message = "Sequence ".$seqObj->id()." is length ".$seqObj->length().". I will create:";
				$message .= "\n  * ".($nb_seq_here-1)." sequences of $opt_size_seq nucleotides" if ($nb_seq_here-1 > 0);
				$message .= "\n  * 1 sequence of ".$leftover." nucleotides\n";
				print $message if ($verbose > 0);
			}
		}
		elsif($opt_leftover eq "remove"){
			if ($leftover > 0 and $leftover < $opt_size_seq){
				print "opt_leftover = remove last piece is short we will remove it ($leftover)\n" if ($verbose > 1);
				$nb_seq_here--;
			}

				my $message = "Sequence ".$seqObj->id()." is length ".$seqObj->length().". I will create:";
				$message .= "\n  * $nb_seq_here sequences of $opt_size_seq nucleotides" if ($nb_seq_here > 0);
				$message .= "\n  * throw out one sequence of length $leftover\n";
				print $message if ($verbose > 0);
		}
		else{ # ($opt_leftover eq "detach")
				my $message = "Sequence ".$seqObj->id()." is length ".$seqObj->length();
				$message .= "\n  * ".($nb_seq_here-1)." sequences of $opt_size_seq nucleotides" if ($nb_seq_here > 0);
				$message .= "\n  * 1 sequence of ".$leftover." nucleotides\n";
				print $message if ($verbose > 0);
		}
		$nb_seq += $nb_seq_here;
		print "adding $nb_seq_here total sequences is now $nb_seq\n" if ($verbose > 1);
	}
	$nb_seq_original++;
	$nb_nucl += $seqObj->length();
}

# ----- Summerize info ------
print "=> number of sequences in the file $nb_seq_original\n";
if ($opt_size_seq) {
	if ($opt_overlap){
		print "=> number of sequences after cutting at size $opt_size_seq bp and adding an overlap of size $opt_overlap: $nb_seq\n";
	}
	else{
		print "=> number of sequences after cutting at size $opt_size_seq bp: $nb_seq\n";
	}
}
else{
	$nb_seq = $nb_seq_original;
}
print "=> number of nucleotides  $nb_nucl \n";


# ----- compute final number of chunk and seq by chunck depending the situation

my $nb_seq_by_chunk;
my $nb_chunk;
# ----------- $opt_nb_chunk --------
if ($opt_nb_chunk and !$opt_nb_seq_by_chunk){

		$nb_chunk = $opt_nb_chunk;
		$nb_seq_by_chunk = ceil($nb_seq / $opt_nb_chunk); # round up, the last chunck can contains less seqeunces
		my $leftover = $nb_seq % $nb_seq_by_chunk;
		my $nb_chunk_info = $leftover ? $opt_nb_chunk-1 : $opt_nb_chunk;

		my $message = "=> A) opt_nb_chunks $opt_nb_chunk\n";
		$message .= "  * $nb_chunk_info chunks with $nb_seq_by_chunk sequences\n" if ($nb_seq_by_chunk >= 1) ;
		$message .= "  * 1 chunk with $leftover sequences\n" if ($leftover);
		print $message;
}

# ----------- $nb_seq_by_chunk --------
elsif ($opt_nb_seq_by_chunk and !$opt_nb_chunk){

		$nb_chunk = ceil($nb_seq / $opt_nb_seq_by_chunk); # round up, the last chunck can contains less sequences
		$nb_seq_by_chunk = $opt_nb_seq_by_chunk;
		my $leftover = $nb_seq % $nb_seq_by_chunk;
		my $nb_chunk_info = $leftover ? $nb_chunk-1 : $nb_chunk;

		my $message = "=> B) opt_nb_seq_by_chunk $opt_nb_seq_by_chunk\n";
		$message .= "  * $nb_chunk_info chunks with $opt_nb_seq_by_chunk sequences\n" if ($nb_chunk_info >= 1) ;
		$message .= "  * 1 chunk with $leftover sequences\n" if ($leftover);
		print $message;
}

# ----------- $opt_nb_chunk and $nb_seq_by_chunk --------
elsif ($opt_nb_seq_by_chunk and $opt_nb_chunk){

	$nb_chunk = $opt_nb_chunk;
	$nb_seq_by_chunk = $opt_nb_seq_by_chunk;
	my $leftover = $nb_seq % $nb_seq_by_chunk;
	my $chunk_needed_to_fullfill = ceil($nb_seq / $nb_seq_by_chunk); # round up, the last chunck can contains less seqeunces
  my $seq_needed_to_fullfill =  ceil($nb_seq / $nb_chunk);

	print "/!\\ You are subsampling! $chunk_needed_to_fullfill chunks or $seq_needed_to_fullfill sequences per chunk".
	" is needed to store all sequences\n" if ( $chunk_needed_to_fullfill > $nb_chunk );

	my $message = "=> C) Creating $opt_nb_chunk chunks with a maximum of $opt_nb_seq_by_chunk sequences by chunk\n";
	#print $message;
}

split_fasta($opt_fastafile, $nb_chunk, $nb_seq_by_chunk, $opt_size_seq, $opt_overlap);



print "Check: $nb_seq sequences expected <=> $cpt_seq_out sequences written\n";
if ($nb_seq != $cpt_seq_out){
	print "Something went bad!\n";
}

# ---- END ----
print "\nusage: $0 @copyARGV\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

# ----------- print fasta outputs
sub split_fasta{
	my ($fasta_file, $nb_chunk, $nb_seq, $split_size, $overlap) = @_;

	#print "split_fasta in: nb_chunk=$nb_chunk, nb_seq=$nb_seq, split_size=$split_size, overlap=$overlap\n";

	my $stop_signal = 0;
	my $seq_cpt=0;
	my $chunk_cpt=1;

	# create output $ostream
	my $ostream = create_fh($chunk_cpt);

	#read fasta file $seqObj per $seqObj
	my $fasta1  = Bio::SeqIO->new(-file => $fasta_file , -format => 'Fasta');
	while ( my $seqObj = $fasta1->next_seq() ) {
		print "dealing with seqid ".$seqObj->id()."\n" if ($verbose > 1);

		# -------------------------------------------------
		# -------------- Split option ---------------------
		# -------------------------------------------------
		if($split_size){

			# ---- No need to split -----
			if( $seqObj->length < $split_size){
					($ostream, $chunk_cpt, $seq_cpt, $stop_signal) = print_fasta($ostream, $seqObj, $chunk_cpt, $nb_chunk, $nb_seq, $seq_cpt);
				last if $stop_signal;
			}

			# ---- Need to split -----
			else{
				my $part=0;
				my $left=0;
				my $right=0;
				my $size = $split_size;
				my $continue = 1;
				while ($continue){
					$part++;

					#Last piece
					if($left + $split_size >= $seqObj->length){
						print "Last chunck\n";
						if ($opt_leftover eq "remove"){
							print "remove $left - ".$seqObj->length."\n";
							last;
						}
						else{ # here opt_leftover eq "detach" OR [ $opt_leftover eq "attach" AND remaining piece is = to $split_size ]
							$right = $seqObj->length;
							$size = $seqObj->length - $left;
							$continue = undef;
						}
					}
					# Not last piece
					else{
						$right = $left + $split_size;
					}

					if($opt_leftover eq "attach"){

						my $next_left_piece = $overlap ? $right-$overlap : $right;
						my $next_right_piece = $next_left_piece + $split_size;
						#print "next_left_piece $next_left_piece\n";
						#print "next_right_piece $next_right_piece\n";
						#print "seqObj->length ".$seqObj->length."\n";
						if($next_left_piece <  $seqObj->length){
							if($next_right_piece >= $seqObj->length){

								my $next_size = $seqObj->length - $next_left_piece;
								print "next_size = $next_size\n";

								if ($next_size < ceil($split_size/2) ){

									$size = $size+$next_size;
									$continue = undef;
									print "attach: attach because next_size $next_size < $split_size / 2(ceil) \n";

								}
								else{
									print "attach: split because next_size $next_size > $split_size / 2(ceil) \n";
								}
							}
						}
					}

					print "substring $left, $size\n" if ($verbose > 1);
					my $sub_seq = substr $seqObj->seq(), $left, $size;
					my $new_seqObj  = Bio::Seq->new( '-format' => 'fasta' , -id => $seqObj->id()."_part$part", -seq => $sub_seq);

					($ostream, $chunk_cpt, $seq_cpt, $stop_signal) = print_fasta($ostream, $new_seqObj, $chunk_cpt, $nb_chunk, $nb_seq, $seq_cpt);
					last if $stop_signal;

					$left = $overlap ? $right-($overlap) : $right; # make one move according to overlap activated or not

				}
				last if $stop_signal;
			}
		}

		# -------------------------------------------------
		# -------------- No split option ---------------------
		# -------------------------------------------------
		else{
			($ostream, $chunk_cpt, $seq_cpt, $stop_signal) = print_fasta($ostream, $seqObj, $chunk_cpt, $nb_chunk, $nb_seq, $seq_cpt);
			last if $stop_signal;
		}
	}
}

sub print_fasta{
	my ($ostream, $seqObj, $chunk_cpt, $nb_chunk, $nb_seq, $seq_cpt) = @_;

	if ($seq_cpt >= $nb_seq){
		$seq_cpt=0; #re-initialize seq counter
		$chunk_cpt++; # increment chunck counter
		if($chunk_cpt > $nb_chunk){
			print "we reach the number of chuncks asked for.\n" if ($verbose > 1);
			return $ostream, $chunk_cpt, $seq_cpt, 1;
		}
		$ostream->close();
		$ostream = create_fh($chunk_cpt); # create new output with proper chunk info
	}
	print "print one seqObj id=".$seqObj->id()." length ".$seqObj->length."\n" if ($verbose > 1);
	$cpt_seq_out++;
	$ostream->write_seq($seqObj);
	$seq_cpt++;
	return $ostream, $chunk_cpt, $seq_cpt, 0;
}

sub create_fh{
	my ($chunk_cpt) = @_;
	# $out_prefix is available from everywhere

	my $name_out = $out_prefix."_chunk".$chunk_cpt.".fa";
	print "output stream is now $name_out\n" if ($verbose > 1);
  open(my $fh, '>', $name_out) or die "Could not open file '$name_out' $!";
  my $ostream = Bio::SeqIO->new(-fh => $fh, -format => 'Fasta' );

	return $ostream;
}

__END__

=head1 NAME

gaas_fasta_splitter.pl

=head1 DESCRIPTION

This script filter sequences by size. It will remove from the output all sequences under a certain size (1000 bp/aa by default)
We keep all sequences >= --size

=head1 SYNOPSIS

    gaas_fasta_splitter.pl -f infile.fasta [ -o outfile ]
    gaas_fasta_splitter.pl --help

=head1 OPTIONS

=over 8

=item B<-f> or B<--fasta>

Input fasta file.

=item B<--nb_chunks>

Integer - Split the multi fasta in such amount of file. Sequences will be distributed evenly.

=item B<--nb_seq_by_chunk>

Integer - Split the multi fasta in order to get such amount of sequences within each file.

=item B<--size_seq>

Integer - split sequences to be at this maximum size. When split, the suffix
_partX will be added, where X in an incremented integer.

=item B<--overlap>

Integer - overlaping part to keep between the split sequences. Used only when <--size_seq> option activated.

=item B<--leftover>

[value: attach(default),remove,detach]. Used only when <--size_seq> option activated.
When splitting a sequence the last chunk can be very short.
attach: When the last chunk is more than 50% shorter of the --size_seq value we attach it to the previous chunk.
        In such case the last chunk of a split sequence can have as (rounded superior)
				maximum length: size_seq + (size_seq / 2) - 1
				minimum length: (size_seq / 2)
remove: The last chunk will be throw away if shorter than --size_seq value.
detach: The last chunck stay as it is. In the worse case it can be 1 nucleotide long.

=item B<-o> or B<--output>

Output folder. Default split_result

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/GAAS/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/GAAS/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat
