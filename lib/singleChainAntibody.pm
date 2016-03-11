package singleChainAntibody;
use strict;
use Carp;
#use warnings;
use SFPerlVars;
use Data::Dumper;

use Exporter qw (import);
our @EXPORT_OK = qw (
                        processSingleChainAntibody
                );
use antibodyAntigen qw (
                           antibodyNumbering
                           makeFreeAntibodyComplex
                           processAntibodyAntigen
                   );

use antibodyProcessing qw (
                              getChainTypeWithChainIDs
                              splitPdb2Chains
                              checkAntigenChains
                              processHapten
                              movePDBs
                      );

sub processSingleChainAntibody
{
    my ($pdbId, $pdbPath, $nsch, $ab, $dir, $masterDir, $LOG, $numbering) = @_;

    my $destPro = "$masterDir"."/".$ab."_Protein_".$numbering;
    my $destNonPro = "$masterDir"."/".$ab."_NonProtein_".$numbering;
    my $destFreeAb = "$masterDir"."/".$ab."_Free_".$numbering;
    my @singleChainAb;
    my $numberingError = 0;
    
    my ($chainType_HRef, $chainIdChainTpye_HRef) =
        getChainTypeWithChainIDs ($pdbPath);

    print {$LOG} "Antibody chain types with chain IDs: \n";
    print {$LOG} Dumper ($chainType_HRef);
    print {$LOG} "Antibody chain IDs, chain Labels and chain sequences: \n";
    print {$LOG} Dumper ($chainIdChainTpye_HRef);

    my %chainType = %{$chainType_HRef};
    my $count = 1;
    splitPdb2Chains($pdbPath);
    print {$LOG} "The $pdbId PDB has been splitted in to different files".
        " based on number of chains\n";

    if ( $ab eq "L") {
        @singleChainAb = @{$chainType{Light}};
    }
    elsif ( $ab eq "H") {
        @singleChainAb = @{$chainType{Heavy}};
    }

    my @antigenIds =
        checkSingleChainRedundancy ($chainType_HRef, $chainIdChainTpye_HRef, $LOG, $ab);

    eval { antibodyNumbering ( \@singleChainAb, $nsch );
           1;
       };

    if ($@ ) {
        print {$LOG} "Kabat Numbering Failed for $pdbId\n" .
            $@ . "Program exited\n\n";
        $numberingError = 1;
        return $numberingError;
        next;
    }
    else {
        print {$LOG} "All the antibodies in $pdbId has been numbered by".
            " antibody numbering program\n";
    }
    # Check for haptens bound with CDRs
    my $hapten = hasHaptenSingleChain ($pdbPath, \@singleChainAb, $chainIdChainTpye_HRef);
    my $fileType;
    my %fileTypeH;
    
    # Checks for haptens and move them to non-Protein data 
    if ( ($hapten) and (!@antigenIds) )
    {

        %fileTypeH = processHapten($pdbPath, \@singleChainAb, $ab);
        
        makeFreeAntibodyComplex ($pdbId, $pdbPath, \@singleChainAb, $count, $fileType,
                                 $dir, $chainIdChainTpye_HRef, $numbering,
                                 $LOG, $destNonPro, $destFreeAb, %fileTypeH);
#        movePDBs ($dir, $destNonPro, $pdbId);
 #       print {$LOG} "This antibody is bound with hapten -- Moved to non- ".
  #          "protein data\n";
    }    

    # Checks for protein antigens (Further checking is within
    # processAntibodyAntigen subroutine)
    elsif ( (@antigenIds) and (!$hapten) )
    {
        $fileType = "num";
        $numberingError = processAntibodyAntigen($pdbId, $pdbPath, $ab, \@antigenIds,
                               \@singleChainAb, $fileType, $dir, $masterDir,
                               $LOG, $chainIdChainTpye_HRef, $destPro,
                               $destNonPro, $destFreeAb, $numbering, %fileTypeH);
        
    }
    # Antibody bound with antigen and haptens - Moved to protein antigen data
    elsif ( (@antigenIds) and ($hapten) )
    {
#        $fileType = "hap";
        %fileTypeH =  processHapten($pdbPath, \@singleChainAb, $ab);
        $numberingError = processAntibodyAntigen($pdbId, $pdbPath, $ab, \@antigenIds,
                               \@singleChainAb, $fileType, $dir, $masterDir,
                               $LOG, $chainIdChainTpye_HRef, $destPro,
                               $destNonPro, $destFreeAb, $numbering, %fileTypeH);
    }
    # Checks for free light and heavy chains
    else {
        $fileType = "num";
        makeFreeAntibodyComplex($pdbId, $pdbPath, \@singleChainAb, $count,
                                $fileType, $dir, $chainIdChainTpye_HRef, $numbering,
                                $LOG, $destNonPro, $destFreeAb, %fileTypeH);
#        movePDBs ($dir, $destFreeAb, $pdbId);
 #       print {$LOG} "This antibody is free chain (light-heavy) without any type of ".
  #          "bound antigen -- Moved to Free chain Data\n";
     
    }
    return $numberingError;
    
}

sub hasHaptenSingleChain
{
    my ($pdbPath, $singleChainAb_ARef, $chainIdChainTpye_HRef) = @_;
    my $hashapten = $SFPerlVars::hashapten;

    my @singleChains = @{$singleChainAb_ARef};
    my %chainIdChainTpye = %{$chainIdChainTpye_HRef};
    
    my ($hapten, @haptens, @allHaptens);
    my $chainLabel;
    # Mapping chain label for chain ID
    for my $chain (@singleChains )
    {
        foreach my $key ( keys %{$chainIdChainTpye{$chain}} ) {
            $chainLabel = $key;
        }

        if ( $chainLabel eq "L")
        {
            @haptens = `$hashapten -l $chain $pdbPath`;
            push (@allHaptens , @haptens );
        }
        
        elsif ( $chainLabel eq "H")
        {
            @haptens = `$hashapten -h $chain $pdbPath`;
            push (@allHaptens , @haptens );
        }

    }
    
        if (grep (/^HAPTEN/, @allHaptens) )
        {
            $hapten = 1;
        }
        else 
            {
                $hapten = 0;
            }
            return $hapten;

    }
  
sub checkSingleChainRedundancy
{
    my ($chainType_HRef, $chainIdChainTpye_HRef, $LOG, $ab) = @_;

    
    my %chainType = %{$chainType_HRef};
    my  @antigen = checkAntigenChains($LOG, %chainType);

    my  @light =  @{$chainType{Light}};
    my  @heavy =  @{$chainType{Heavy}};
    my $abchainAg;
    
    my @allChains = (@light, @heavy);
    
    my %chainIdChainTpye = %{$chainIdChainTpye_HRef};

    for ( my $i = 0; $i<= $#allChains - 1; $i++)
    {
        for ( my $j = $i+1; $j<= $#allChains; $j++)
        {
            if ( $ab eq "L")
            {
                if ( index ($chainIdChainTpye{$allChains[$i]}{"L"},
                     $chainIdChainTpye{$allChains[$j]}{"L"}) == -1)
                 {
                     $abchainAg = $allChains[$j];
                     push (@antigen, $allChains[$j]);
                      print {$LOG} $allChains[$i]." antibody chain is found as non-redundant ".
                         "to other chains in PDB\n";
                 }
                else {
                    
                    next;
                }
             }
            elsif ($ab eq "Heavy")
            {
                if ( index ($chainIdChainTpye{$allChains[$i]}{"H"},
                            $chainIdChainTpye{$allChains[$j]}{"H"}) == -1)
                    {
                        $abchainAg = $allChains[$j];
                        push (@antigen, $allChains[$j]);
                        print {$LOG} $allChains[$i]." antibody chain is found as non-redundant ".
                            "to other chains in PDB\n";
                    }
                else {
                    
                    next;
                }
            }

        }
    }
    return @antigen;
}
        
