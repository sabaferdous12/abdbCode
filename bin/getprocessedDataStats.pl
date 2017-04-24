use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Data::Dumper;
use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

use general qw (readDirPDB);
# Only read directories with Martin name
opendir(DIR, ".") or die "Error in opening directory";
my @files = grep(/Martin/,readdir(DIR) );
closedir(DIR);

my @ProcessedFiles;
my $resultantFiles;
my %stats;
foreach my $subDir ( @files )
{
    if ( $subDir =~ /bz2/) {
        next;
    }   
    my @dataDir = readDirPDB($subDir);
    $resultantFiles = scalar @dataDir;
    
    if ($subDir =~ /NR/)
        {
            $stats{"Resultant_$subDir"} = $resultantFiles;
            next;
        }    
   
    foreach my $file (@dataDir)
    {
        my ($pdbId, $ext) =  split ("_", $file);
        push (@ProcessedFiles, $pdbId);
    }
#    print "Original: ", scalar @ProcessedFiles, "\n";
    
    @ProcessedFiles = uniq @ProcessedFiles;

#    print "Unique: ", scalar @ProcessedFiles, "\n";
#    print join ("\n", @ProcessedFiles), "\n";

    
    $stats{"Processed_$subDir"} = scalar @ProcessedFiles;
    $stats{"Resultant_$subDir"} = $resultantFiles;
    
    @ProcessedFiles = ();
}

#print Dumper (\%stats);

system ("bash ~/git/abdbCode/bin/statsProcessed.sh $stats{'Processed_LH_Protein_Martin'} $stats{'Resultant_LH_Protein_Martin'} $stats{'Resultant_NR_LH_Protein_Martin'} $stats{'Processed_LH_NonProtein_Martin'} $stats{'Resultant_LH_NonProtein_Martin'} $stats{'Resultant_NR_LH_NonProtein_Martin'} $stats{'Processed_LH_Free_Martin'} $stats{'Resultant_LH_Free_Martin'} $stats{'Resultant_NR_LH_Free_Martin'} $stats{'Processed_LH_Combined_Martin'} $stats{'Resultant_LH_Combined_Martin'} $stats{'Resultant_NR_LH_Combined_Martin'} $stats{'Processed_L_Protein_Martin'} $stats{'Resultant_L_Protein_Martin'} $stats{'Resultant_NR_L_Protein_Martin'} $stats{'Processed_L_NonProtein_Martin'} $stats{'Resultant_L_NonProtein_Martin'} $stats{'Resultant_NR_L_NonProtein_Martin'} $stats{'Processed_L_Free_Martin'} $stats{'Resultant_L_Free_Martin'} $stats{'Resultant_NR_L_Free_Martin'} $stats{'Processed_L_Combined_Martin'} $stats{'Resultant_L_Combined_Martin'} $stats{'Resultant_NR_L_Combined_Martin'} $stats{'Processed_H_Protein_Martin'} $stats{'Resultant_H_Protein_Martin'} $stats{'Resultant_NR_H_Protein_Martin'} $stats{'Processed_H_NonProtein_Martin'} $stats{'Resultant_H_NonProtein_Martin'} $stats{'Resultant_NR_H_NonProtein_Martin'} $stats{'Processed_H_Free_Martin'} $stats{'Resultant_H_Free_Martin'} $stats{'Resultant_NR_H_Free_Martin'} $stats{'Processed_H_Combined_Martin'} $stats{'Resultant_H_Combined_Martin'} $stats{'Resultant_NR_H_Combined_Martin'} >stats_processed.tt");
