#!/usr/bin/perl -s
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

my $differenceFile = $ARGV[0];
my $redundantCombinedFile = $ARGV[1];
my ($antigenDir, $haptenDir);


open ($DIFF, '<', $differenceFile) or
    die "Can not open file $!\n";

open ($COMBINED, '<', $redundantCombinedFile) or
    die "Can not open file $COMBINED\n";
my @combinedClus = <$COMBINED>;

if ( defined ( $::a) ){
    open ( $OUT, '>', "FreeAntibody_AntibodyAntigen.txt")
        or die "Can not open file\n";
    print $OUT "Free Antibody:Complexed Antibody\n";
    $antigenDir = "AntibodyAntigen_Martin";
    $haptenDir = "AntibodyHapten_Martin";
    getFreeVsComplexedList ($DIFF, $OUT, $antigenDir, $haptenDir);    
    
}     
if ( defined ( $::l) ){
    open ( $OUT, '>', "Light_LightAntigen.txt")
        or die "Can not open file\n";
    print $OUT "Free Bence Jones (Light Chains):Complexed Bence Jones\n";
    $antigenDir = "LightAntigen_Martin";
    $haptenDir = "LightHapten_Martin";
    getFreeVsComplexedList ($DIFF, $OUT, $antigenDir, $haptenDir);
}
if  ( defined ( $::h) ) {
    open ( $OUT, '>', "Heavy_HeavyAntigen.txt")
        or die "Can not open file\n";
    print $OUT "Free Camelids (Heavy Chains):Complexed Camelids\n";
    $antigenDir = "HeavyAntigen_Martin";
    $haptenDir = "HeavyHapten_Martin";
    getFreeVsComplexedList ($DIFF, $OUT, $antigenDir, $haptenDir);
}
             

     
sub getFreeVsComplexedList{
    my ($DIFF, $OUT, $antigenDir, $haptenDir) = @_;
    my (@free, @complex);
    my $elemPdb;
    
    while (my $line = <$DIFF>)
        {
            chomp $line;
            my ($records) = grep (m/$line/, @combinedClus);
            my @redundants = split (/,\s+/, $records);
            
            foreach my $elem (@redundants)
                {
                    chomp $elem;
                    $elemPdb = $elem.".pdb";
                    
                        if ( (-e "./$antigenDir/$elemPdb") or
                                 (-e "./$haptenDir/$elemPdb")) 
                            {
                                push (@complex, $elem);
                            }
                    else
                        {
                            push (@free, $elem);
                        }
                }
            # To skip any cluster without free antibody
            if ( (!@complex ) or (!@free)) {
                next;
            }
            else {
                print $OUT join(',', @free), ":", join(',', @complex), "\n";
            }
            
            @free = ();
            @complex = ();
        }
}

