        #HANDLE format
        my $format=3;
        my $problem3=undef;
        my $nbLineChecked=1000; #number line to use to check the formnat
        my $cpt=0;
         open(my $fh, '<', $file) or die "cannot open file $file";
    {
        while(<$fh>){
                $cpt++;
                if($cpt > $nbLineChecked){
                        _printSurrounded("Dosn't look as a GFF file\nLet's see what the Bioperl parser can do with that...",100,"!");  
                        last;
                }
                if($_ =~ /^.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t(.*)/){
                        if(length($1) < 1){next;}
                                if($1 !~ /=/  and $1 !~ /;/ ){  
                                $format=1;
                                _printSurrounded("Problem detected wihtin the 9th colum of the gff file\nYou cannot have space between attributes and between tag and values.\nAll your attributes will be gathered within the GROUP tag",100,"!");  
                                last;
                        }
                        elsif($1 !~ /=/  and $1 =~ /;/ ){       
                                        $format=2;
                                        last;
                        }
                        #my $y;
                        #my $c = () = $y =~ /=/g;
                        #my $d = () = $y =~ /\ /g;
                        #if($c > 1 and $d>1  and $1 !~ /;/ ){
                        #       $problem3=1;
                        #}
                }
        }
    }
    close($fh);
    if($problem3){
        _printSurrounded("Thre is a problem with your GFF format.\nThis format is wrong: tag=value tag=value.\nYou should have: tag=value;tag=value or tag value ; tag value",100,"!");  
    }
        print "GFF version parser used: $format\n";
