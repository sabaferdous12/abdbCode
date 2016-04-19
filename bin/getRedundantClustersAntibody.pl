#!/usr/bin/perl
#*************************************************************************
#
#   Program:    getRedundantAntibodyClusters
#   File:       getRedundantAntibodyClusters.pl
#   
#   Version:    V2.0
#   Date:       09.12.15
#   Function:   For the given directory for antibody PDB structures/files, 
#               computes the clusters of redundant antibodoies
#   
#   Copyright:  (c) UCL, Saba Ferdous, 2014
#   Author:     Miss Saba Ferdous
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      saba@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#   The script finds the redudant clusters of the antibodies with identical 
#   sequences. It uses external script abnr.pl (~martin/scripts/abnr/abnr.pl)
#   that take input 2 antibody numbered files and returns the label SAME, if
#   they are redundant and DIFFERENT otherwise.
#*************************************************************************
#
#   Usage:
#   ======
#  ./getRedundantAntibodyClusters.pl 
#   (Script should be run from the desired directory)
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.1   22.04.15 Original
#
#*************************************************************************

use strict;
#use warnings;
use Data::Dumper;
use Cwd;
use SFPerlVars;

use File::Basename;
use general qw (readDirPDB);


use redund qw (
                  readDir
                  printClusters 
          );

# Reads all the PDB files from current directory and list them in an array
my $dir = getcwd;
my $dirName = basename($dir);
# Output file has same name as directory name with a prefix of Redundant
open(my $OUT, ">", "Redundant_".$dirName.".txt") or 
    die "Can not write file $!";

my @antibodyList = readDirPDB ($dir);
@antibodyList = sort @antibodyList; 

# Variable Declaration
my (%redundantElements, %redundantClusters, %unique) = ();
my @identicals = ();
my ($rd, $tag, $abnrProg);
$abnrProg = $SFPerlVars::abnr;
my $ab2;

# Reads every element of the antibody list and compares with rest
for (my $ab1 = 0; $ab1 <= $#antibodyList ; $ab1++)
{
    # Checks if the antibody has already been compared and present in previous
    # cluster
    next if (exists $redundantElements{$antibodyList[$ab1]});
      
    # Inner loop to compare the rest of element with the outer loop elemenet
    for (my $ab2 = ($ab1+1) ; $ab2 <= $#antibodyList; ++$ab2)
    {
	# Checks if the antibody has already been compared
        next if(exists $redundantElements{$antibodyList[$ab2]});
        
        my $output = `perl $abnrProg $antibodyList[$ab2] $antibodyList[$ab1]`;        

        # if redundant program will return "SAME 1AFV_2.pdb"
        if ( $output =~ /SAME/) 
        {
            
#            ($tag, $rd) = split(" ", $output); 
            push( @identicals, $antibodyList[$ab2] );
            $redundantElements{$antibodyList[$ab2]} = 1; # Flag the redundants
        }
        else {
            next;
        }
    }
  
    if (@identicals)
    {                                         
	# hash with an array of redundant PDBs;     
	$redundantClusters{$antibodyList[$ab1]} = [@identicals];
	@identicals = ();                         
    }
    else 
    {                                                                        
	# hash with PDBs with no redundant PDB
	$unique{$antibodyList[$ab1]} = 1;     
	next;
    }

}

printClusters (\%redundantClusters, \%unique, $OUT);


=cut
=head1 NAME                                                                   


=head1 SYNOPSIS                                                               
use this program, type on command line,                                       
./getRedundantAntibodyClusters.pl

=head1 SYNOPSIS  
Perl Module for generating clusters for redundant antibodies                  

=head1 AUTHOR                                                                 
Saba Ferdous (ucbterd@acrm19)                                                 

=head1 COPYRIGHT AND LICENSE                                                  
Copyright (C) 2014 by Saba Ferdous                                            

This program is free software; you can redistribute it and/or modify          
it under the same terms as Perl itself, either Perl version 5.8.2 or,         
at your option, any later version of Perl 5 you may have available.           

=head1 BUGS                                                                   
                         
None reported... yet.

=cut




