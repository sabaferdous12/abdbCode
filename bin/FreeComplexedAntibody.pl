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
#use warnings;

my ($DIFF, $COMBINED, $OUT);

my $differenceFile = $ARGV[0];
my $redundantCombinedFile = $ARGV[1];
my ($proAntigenDir, $nproAntigenDir);


open ($DIFF, '<', $differenceFile) or
    die "Can not open file for reading... $!\n";

open ($COMBINED, '<', $redundantCombinedFile) or
    die "Can not open file $COMBINED\n";
my @combinedClus = <$COMBINED>;
my $dir = ".";

if ( defined ( $::LH) ){
    open ( $OUT, '>', "$dir/FreeAntibody_AntibodyAntigen.list")
        or die "Can not open file to write\n";
    print $OUT "Free Antibody:Complexed Antibody\n";
    $proAntigenDir = "LH_Protein_Martin";
    $nproAntigenDir = "LH_NonProtein_Martin";
    getFreeVsComplexedList ($DIFF, $OUT, $proAntigenDir, $nproAntigenDir);    
    
}     
if ( defined ( $::L) ){
    open ( $OUT, '>', "$dir/Light_LightAntigen.list")
        or die "Can not open file to write\n";
    print $OUT "Free Bence Jones (Light Chains):Complexed Bence Jones\n";
    $proAntigenDir = "L_Protein_Martin";
    $nproAntigenDir = "L_NonProtein_Martin";
    getFreeVsComplexedList ($DIFF, $OUT, $proAntigenDir, $nproAntigenDir);
}
if  ( defined ( $::H) ) {
    open ( $OUT, '>', "$dir/Heavy_HeavyAntigen.list")
        or die "Can not open file to write\n";
    print $OUT "Free Camelids (Heavy Chains):Complexed Camelids\n";
    $proAntigenDir = "H_Protein_Martin";
    $nproAntigenDir = "H_NonProtein_Martin";
    getFreeVsComplexedList ($DIFF, $OUT, $proAntigenDir, $nproAntigenDir);
}
             

     
sub getFreeVsComplexedList{
    my ($DIFF, $OUT, $proAntigenDir, $nproAntigenDir) = @_;
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

            if ( (-e "./Data/$proAntigenDir/$elemPdb") or
                     (-e "./Data/$nproAntigenDir/$elemPdb")) 
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

