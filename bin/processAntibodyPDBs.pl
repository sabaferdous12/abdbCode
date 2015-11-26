#!/usr/bin/perl -s
# processAntibodyPDBs.pl --- 
# Author: Saba Ferdous <ucbterd@martin-cd01.biochem.ucl.ac.uk>
# Created: 02 Nov 2015
# Version: 0.01

use warnings;
use strict;
use SFPerlVars;
use Data::Dumper;
use File::Copy;
use Cwd;

use general qw (readDirPDB);
use antibodyAntigen qw (processAntibody);
use singleChainAntibody qw (
                               processSingleChainAntibody
                       );

use antibodyProcessing qw (
                              getPDBPath
                              getChains
                              readFileDataInArray
                              dirOperations
                      );
my ($nsch, $numbering);

if ( $::k )
{
     $nsch = '-k';
     $numbering = "Kabat";
 }

if ( $::c )
{
     $nsch = '-c';
     $numbering = "Chothia";
 }

if ( $::a )
{
     $nsch = '-a';
     $numbering = "Martin";
 }

my $inputFile = $ARGV[0];
my @PDBCodes = readFileDataInArray($inputFile);

my $masterDir = getcwd;
my $processDir = $masterDir.'/'.'Dataprep' . $$ ."_". $numbering; # $$ = process id
my $dir;

my ($dataPrep, $proteinAntigensAB, $nonProteinAntigensAB, $freeAntibodiesAB)
    = ("DataPrep", "LH_Protein_".$numbering,
       "LH_NonProtein_".$numbering, "LH_Free_".$numbering);
my ($proteinAntigensLG, $nonProteinAntigensLG, $freeAntibodiesLG)
    = ("L_Protein_".$numbering,
       "L_NonProtein_".$numbering, "L_Free_".$numbering);
my ($proteinAntigensHV, $nonProteinAntigensHV, $freeAntibodiesHV)
    = ("H_Protein_".$numbering,
       "H_NonProtein_".$numbering, "H_Free_".$numbering);

`mkdir -p $proteinAntigensAB $nonProteinAntigensAB $freeAntibodiesAB`;
`mkdir -p $proteinAntigensLG $nonProteinAntigensLG $freeAntibodiesLG`;
`mkdir -p $proteinAntigensHV $nonProteinAntigensHV $freeAntibodiesHV`;

my ($abagCount, $lightCount, $heavyCount, $fcCount,
    $supersededCount, $numErrorCount)
    = (0, 0, 0, 0, 0, 0);
my ($ab, $numberingError);

open ( my $HEADER, ">header.dat") or die "Can not open file $!";
open ( my $AGCHAIN, ">AntigenChains.dat") or die "Can not open file $!";

print {$AGCHAIN} "PDB_ID:Antgen Chains\n";
open ( my $ABCHAIN, ">AntibodyChains.dat") or die "Can not open file $!";
print {$ABCHAIN} "PDB_ID:Antibody Chains\n";

foreach my $pdb ( @PDBCodes )
{
    chomp $pdb;
    print "\nWorking on $pdb\n";
    my $logFile =  "antibody.log";
    
    $dir = dirOperations ($processDir, $pdb);
    my $pdbPath = getPDBPath($pdb);

    open (my $LOG, '>', $logFile) or die "Can not open file\n";
    print {$LOG} "$pdb\n";
    print {$LOG} "The file path for $pdb is obtained: $pdbPath\n";

    # checks for superseded PDB files
    if ( !( -e "$pdbPath" ) ) {
#        push ( @superseded, $pdbId );
        #print {$MASTER_LOG} "$pdbId is found Superseded\n\n";
        $supersededCount++;
        next;
    }
    # Get antibody chain info
    my ($light_ARef, $heavy_ARef, $antigen_ARef) =
        getChains ($pdbPath);
    my @light = @{$light_ARef};
    my @heavy = @{$heavy_ARef};
    my @antigens = @{$antigen_ARef};

    print {$ABCHAIN} $pdb,":",join (",", @light), ",";
    print {$ABCHAIN} join (",", @heavy), "\n";
    print {$AGCHAIN} $pdb,":",join (",", @antigens), "\n";
        
    if ( ( @light ) and (@heavy) )
        {
            $ab = "LH";
            $numberingError =
                processAntibody ($pdb, $pdbPath, $nsch, $ab, $dir, $masterDir, $LOG, $numbering);
            if ( $numberingError)
            {
                $numErrorCount ++;
                next;
            }
            $abagCount++;
        }
   ################# Just light
    elsif ( ( @light ) and ( !@heavy ) )
        {
            $ab = "L";
            $numberingError =
                processSingleChainAntibody($pdb, $pdbPath, $nsch, $ab, $dir, $masterDir, $LOG, $numbering);
            if ( $numberingError)
            {
                $numErrorCount ++;
                next;
            }
            $lightCount++;
        }
    ################# Just heavy
    elsif ( ( !@light ) and (@heavy ) )
        {
            $ab = "H";
            $numberingError =
                processSingleChainAntibody($pdb, $pdbPath, $nsch, $ab, $dir, $masterDir, $LOG, $numbering);
            if ( $numberingError)
            {
                $numErrorCount ++;
                next;
            }
            
            $heavyCount++;
        }
    else {
        $fcCount++;
    }
    
    # Writing the header information into a flat file
    my @headerInfo =  `pdbheader -p -s -m $pdbPath`;
    print {$HEADER} @headerInfo;
    

chdir '..';
#last;

}

print "Processed - ABAG: $abagCount\n";
print "Processed - Light: $lightCount\n";
print "Processed - Heavy: $heavyCount\n";
print "Processed - FC:  $fcCount\n";
print "Processed - Superseded: $supersededCount\n";
print "Kabat Error - : $numErrorCount\n"

__END__

=head1 NAME

processAntibodyPDBs.pl - Describe the usage of script briefly

=head1 SYNOPSIS

processAntibodyPDBs.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for processAntibodyPDBs.pl, 

=head1 AUTHOR

Saba Ferdous, E<lt>ucbterd@martin-cd01.biochem.ucl.ac.ukE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Saba Ferdous

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
