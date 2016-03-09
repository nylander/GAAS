#!/usr/bin/perl -w

package BILS::CheckModule;

use strict;
use Exporter qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(module_software_installed module_avail module_unload module_version_avail module_load module_version_is_loaded module_is_loaded);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                         Ok    => [qw(module_software_installed module_avail module_unload module_version_avail module_load module_version_is_loaded module_is_loaded)]);
=head1 SYNOPSIS




=head1 DESCRIPTION

        A library to convert various annotation files
        to GFF3 format. This is the core file in which we
        specify all required resources.

=cut

sub _get_module_avail_hash{
	my $module_avail=`eval "module avail 2>&1"`;
	
	my %modules;
	# cut by line
	my @lines_modules = split/\n/,$module_avail ;
	foreach	my $line (@lines_modules){
        	if($line =~ /^[^-]/){
	        	my @module_list=split /\s/, $line;
                	foreach	my $one_module (@module_list){
                        	my @info=split /\//, $one_module;
                        	if(exists($info[1]) ) {
                        	        $modules{lc($info[0])}{$info[1]}++; 
                        	}
                        	else{
                        	        $modules{lc($info[0])}{-1}++;
                        	}
                	}
        	}
	}
	return \%modules;
}

# Check if the module software is installed
sub module_software_installed{
       
	my $call=`eval "module 2>&1"`;
	if($call =~/command not found/){
		return undef;
	}
	return 1
}

# get last version of a tool's module
sub _get_last_version{
	my ($tool,$modulesH)=@_;
	
	$tool=lc($tool);
	if(! exists($modulesH->{$tool})){print "Module $tool doesn't exists.\n";}

	my $hash_tool=$modulesH->{$tool};
	
	my @versions = ( keys %{$hash_tool} );
	# sort vesions
	@versions = sort {versioncmp($a,$b) } @versions;
	foreach my $key (reverse (@versions)){
		return $key;
	}
}

#check if module given exists, if version is not given or not exists, it return the last version.
# return the version if exits 
# return the last version if a version is not specified
# return undef if the tool/module doesn't exists (doesn't care about version checked)  
sub module_avail {
        my ($moduleName) = @_;

        my $result=undef;
        my $modulesH=_get_module_avail_hash();

	my @tmp = split /\// , $moduleName;
	my $tool=lc($tmp[0]);
	my $version=undef;
	if(exists ($tmp[1])){
		$version=$tmp[1];
	}
#	print "version $version\n";

	if( exists ($modulesH->{$tool}) ) {
                if(! $version){
			print "No version provided for the module $moduleName.\n";
			$result= _get_last_version($tool,$modulesH);
			print "this version is available: $result\n";
		}
		elsif(! exists ($modulesH->{$tool}{$version}) ){
			print "version specified $version of the module $tool not found,\n";
			$result= _get_last_version($tool,$modulesH);
			print "this version is available: $result\n";
		}
		else{		
			$result=$moduleName;
        	        print "module $moduleName exits\n";
		}
        }
	else{
                print "module $moduleName doesn't exits\n";		
        }

return $result;
}

#check if module of the version given exists,
# return 1 if succeed.
sub module_version_avail {
        my ($moduleName) = @_;

        my $result=undef;
        my $modulesH=_get_module_avail_hash();

        my @tmp = split /\// , $moduleName;
        my $tool=lc($tmp[0]);
        my $version=undef;
        if(exists ($tmp[1])){
                $version=$tmp[1];
        }

	if( exists ($modulesH->{$tool}) ) {
                if(! $version){
			return 0;
                }
                elsif(! exists ($modulesH->{$tool}{$version}) ){
			return 0;
                }
                else{
			return 1;
                }
        }
	else{
		return 0;
        }
}

#load a module
sub module_load {
	my ($moduleName) = @_;
	
	my $resu = module_avail($moduleName);
	
	my @tmp = split /\// , $moduleName;
        my $tool=lc($tmp[0]);

	if(! $resu){
            print "module $moduleName doesn't exits\n"; 
            return 0;
        }
	else{
             	exec("module load ".$tool."/".$resu."; perl $0 @ARGV");
                print "module $tool$resu loaded\n";
                return 1;
        }

}

#load a module
sub module_unload {
        my ($moduleName) = @_;



        if(! module_avail($moduleName)){
             print "module $moduleName doesn't exits\n";
            return 0;
        }
	else{
             	exec("module unload ".$moduleName."; perl $0 @ARGV");
                print "module $moduleName loaded\n";
                return 1;
        }

}

#create hash of module list (=loaded)
sub _get_module_list_hash {
	my %modules;

	my $module_list=`eval "module list 2>&1"`;
        my @lines_modules = split/\n/,$module_list;
        my $firstLine=1;
	foreach my $line (@lines_modules){
		if($firstLine){
			$firstLine=0;next;
		}
		else{
			 my @tmp=split/\s/, $line;
			 my @info= split /\//, $tmp[3];
			 if(exists($info[1]) ) {
                                $modules{lc($info[0])}{$info[1]}++;
			}
                         else{
                             	$modules{lc($info[0])}{-1}++;
                         }
			
		}
	}
return \%modules;
}


#check if a module is loaded
sub module_is_loaded{
	 my ($moduleName) = @_;	

	 my $hash_loaded=_get_module_list_hash($moduleName); 

	my @tmp = split /\// , $moduleName;
        my $tool=lc($tmp[0]);


	 if(($hash_loaded->{$tool})){
            return 1;
        }
	else{
                return 0;
        }
}

#check if a specific version of a module is loaded
sub module_version_is_loaded{
         my ($moduleName) = @_;

         my $modulesH=_get_module_list_hash($moduleName);

        my @tmp = split /\// , $moduleName;
        my $tool=lc($tmp[0]);
	my $version=undef;
        if(exists ($tmp[1])){
                $version=$tmp[1];
        }

        if( exists ($modulesH->{$tool}) ) {
                if(! $version){ #no version specified so we cannot check for a specific version (use the module_is_loaded method)
                        return 0;
                }
                elsif(! exists ($modulesH->{$tool}{$version}) ){
                        return 0;
                }
                else{
                        return 1; # version of this module is loaded
                }
        }
        else{
                return 0;
        }
}

#==== internal methods
# compare versions 
sub versioncmp ($$) {
    my @A = ($_[0] =~ /([-.]|\d+|[^-.\d]+)/g);
    my @B = ($_[1] =~ /([-.]|\d+|[^-.\d]+)/g);

    my ($A, $B);
    while (@A and @B) {
	$A = shift @A;
	$B = shift @B;
	if ($A eq '-' and $B eq '-') {
	    next;
	} elsif ( $A eq '-' ) {
	    return -1;
	} elsif ( $B eq '-') {
	    return 1;
	} elsif ($A eq '.' and $B eq '.') {
	    next;
	} elsif ( $A eq '.' ) {
	    return -1;
	} elsif ( $B eq '.' ) {
	    return 1;
	} elsif ($A =~ /^\d+$/ and $B =~ /^\d+$/) {
	    if ($A =~ /^0/ || $B =~ /^0/) {
		return $A cmp $B if $A cmp $B;
	    } else {
		return $A <=> $B if $A <=> $B;
	    }
	} else {
	    $A = uc $A;
	    $B = uc $B;
	    return $A cmp $B if $A cmp $B;
	}	
    }
    @A <=> @B;
}

1;
