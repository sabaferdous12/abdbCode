#!/usr/bin/perl
#*************************************************************************
#
#   Program:    FreeComplexedAntibody
#   File:       FreeComplexedAntibody.pl
#   
#   Version:    V1.1
#   Date:       22.04.14
#   Function:   This script provides the list of antibody PDBs (file) that are 
#               present in both forms i.e Free and Complexed               
#   Usage:      ./FreeComplexedAntibody.pl
#               
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
use strict;
use warnings;

my ($DIFF, $COMBINED, $OUT);
open ($DIFF, '<', "difference.txt") or die "Can not open file $DIFF\n";
open ($COMBINED, '<', "./Redundant_files/Redundant_Combined_Martin.txt")
    or die "Can not open file $COMBINED\n";
open ( $OUT, '>', "FreeComplexedAntibody.txt")
  or die "Can not open file\n";
print $OUT "Free Antibody:Complex\n";

my (@free_antibody, @complex);
my @combined = <$COMBINED>;

while (my $line = <$DIFF>)
{
    chomp $line;
    my ($records) = grep (m/$line/, @combined);
    my @redundants = split (/,\s+/, $records);
    
    foreach my $elem (@redundants)
    {
	chomp $elem;
	my $elemPdb = $elem.".pdb";
	if ( (-e "./Complex_Martin/$elemPdb") or
	     (-e "./Hapten_Martin/$elemPdb") )
	{
	    push (@complex, $elem);
	}
	else
	{
	    push (@free_antibody, $elem);
	}
	next if ( (!@complex ) or (!@free_antibody) );
    }
    
    print $OUT join(',', @free_antibody), ":", join(',', @complex), "\n";
    @free_antibody = ();
    @complex = ();
}

