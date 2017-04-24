use strict;
use warnings;
use Data::Dumper;
use File::Copy;
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

use List::Util qw( min max );
use List::MoreUtils qw(uniq);
use antibodyProcessing qw (getResolInfo
                           getPDBPath
                           largestValueInHash);

my $infile = $ARGV[0]; # List of redundant clusters
my $outputDir = $ARGV[1]; # Directory to store non-redundant data

open(my $IN, $infile) or die "Can not open input file$!\n";

my @nrFiles;
my @redundPDBcodes;
my @bestABs;
my %clusterRes;

# Reading each redundant cluster 
while (<$IN>)
    {
        # Split on comma
        my @redundantPDBs = split(/, /, $_);
             
        # Grab PDB code from file names 1AFV_1, make it 1AFV
        foreach my $rAB( @redundantPDBs  )
        {
            chomp $rAB;
            my ($PDB) = split(/_/, $rAB);
            push (@redundPDBcodes, $PDB);
        }
        # Keeping only one of the PDB if there are multiple PDBs
        @redundPDBcodes = uniq @redundPDBcodes;
                
        foreach my $pdb ( @redundPDBcodes )
        {
            # Obtain pdb path to get resolution
            my $pdbPath = getPDBPath ($pdb);
            my %resInfo = getResolInfo ($pdbPath);
            $clusterRes{$pdb} = $resInfo{'Resolution'};
        }
        @redundPDBcodes = ();

        # Select best PDB with minimum resolution
        my $bestPDB = smallestValueInHash(%clusterRes);
        %clusterRes = ();
        # To grep the real file name
        my ($bestABid) = grep (/$bestPDB/, @redundantPDBs);
        push (@bestABs, $bestABid);
#      last;
    }

# Moving NR and low resolution PDBs to NR dataset directory
foreach  my $nrFile (@bestABs)
{
    chomp $nrFile;
    if ( $nrFile =~ /,/)
    {
        $nrFile =~ s/,//;
    }
    my $nrFile = $nrFile.".pdb";
    my $dest = "$outputDir";
    if (-e $nrFile)
    {
        copy ($nrFile, $dest);
    }
    
    else {
        print "$nrFile does not exists\n";
    }
}


sub smallestValueInHash
{
my (%hash) = @_;
my $smallest;
# Selecting smallest value key from the hash
foreach my $name (sort { $hash{$a} <=> $hash{$b} } keys %hash)
{
    $smallest = $name;
    last;
    }
return $smallest;
}
