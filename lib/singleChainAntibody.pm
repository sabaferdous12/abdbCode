package singleChainAntibody;
use strict;
use Carp;
#use warnings;
use config;
use Data::Dumper;
use File::Basename;

use Exporter qw (import);
our @EXPORT_OK = qw (
                        processSingleChainAntibody
                );
use antibodyAntigen qw (
                           antibodyNumbering
                           makeFreeAntibodyComplex
                           processAntibodyAntigen
                           getInterchainContacts
                           antibodyAssembly
                           
                   );

use antibodyProcessing qw (
                              getChainTypeWithChainIDs
                              splitPdb2Chains
                              checkAntigenChains
                              processHapten
                              movePDBs
                              dealNumError
                              mergeLH
                              hasHapten
                      );

sub processSingleChainAntibody
{
    my ($pdbId, $pdbPath, $nsch, $ab, $dir, $masterDir, $LOG, $numbering) = @_;

    my $destPro = "$masterDir"."/".$ab."_Protein_".$numbering;
    my $destNonPro = "$masterDir"."/".$ab."_NonProtein_".$numbering;
    my $destFreeAb = "$masterDir"."/".$ab."_Free_".$numbering;
    my @singleChainAb;
    my $numberingError = 0;
    my @dimers;
    
    my ($chainType_HRef, $chainIdChainTpye_HRef) =
        getChainTypeWithChainIDs ($pdbPath, $pdbId, $LOG);

    print {$LOG} "Antibody chain types with chain IDs: \n";
    print {$LOG} Dumper ($chainType_HRef);
    print {$LOG} "Antibody chain IDs, chain Labels and chain sequences: \n";
    print {$LOG} Dumper ($chainIdChainTpye_HRef);

    my %chainType = %{$chainType_HRef};

#   my @antigenIds =
 #       checkSingleChainRedundancy ($chainType_HRef, $chainIdChainTpye_HRef, $LOG, $ab);

    my $count = 1;

    my $pdbsymm = $config::pdbsymm;
    my $pdbcount = $config::pdbcount;
    
    my $symPdb = "temp.sym";

    # Apply Symmetery operation on PDB 
    `$pdbsymm $pdbPath >$symPdb`;
    print {$LOG} "Symmetery Operation has been applied\n";

    # Count number of chains in original and symmetry/translated PDB  
    my $chainsPdb = `$pdbcount $pdbPath | awk -F " " '{print \$2}'`;
    my $chainsSym = `$pdbcount $symPdb | awk -F " " '{print \$2}'`;

    # Compare the number of chains. If symmetry PDB has more chains than original then
    # all the chains from that PDB will be split
   
    my @PDBchains;
    my $dimer_flag = 0;

    # Do symmtery just for L chain 
    if ( $ab eq "L") {

        if ( ( $chainsSym >= $chainsPdb ) and ($chainsPdb != 1) ){
            $pdbPath = $symPdb;
            print {$LOG} "Symmetry PDB file will be used to split the chains\n";
            $dimer_flag = 1;
        }
        
    }
    elsif ( $ab eq "H") {
        $dimer_flag = 0;
    }
    
         @PDBchains = splitPdb2Chains($pdbPath);
        print {$LOG} "The $pdbId PDB has been splitted in to different files".
            " based on number of chains\n";
    
    my $chainCount = scalar @PDBchains;
    my $numError = antibodyNumbering ( \@PDBchains, $nsch );
    my $countFailedchains = dealNumError ($LOG, @PDBchains);
    
    if ( ($numError) and ($chainCount == $countFailedchains))
        {
            $numberingError = 1;
        }
    my $cdrError=0;
    my @antigenIds =
         checkSingleChainRedundancy ($chainType_HRef, $chainIdChainTpye_HRef, $LOG, $ab);
       
##############    
    # If a single light chain is not redundant then it should not be treated as antigen
    # so excluded that from antigen array
    my $index = 0;
    foreach my $ag(@antigenIds ) {
        
     my  @lightChain = @{$chainType{Light}};
               
        if  (grep {$_ eq $ag} @lightChain) {
            splice @antigenIds, $index, 1;
        
        }
     $index++;
     
    }

##############
    
    ##############    
    # check for Antigen Symmetry file
    # For antigen chains, it checks for similar files on basis of number
    # of lines and then updates the antigenIds array

    if ( $dimer_flag) {

        my ($fLines, $AG, @tempArr);
        foreach $AG( @antigenIds) {
            my $f = $AG.".pdb";
            $fLines = `wc -l < $f`;
            
            my ($ch, $chLines);
            foreach $ch ( @PDBchains ) {
                my $chF = $ch.".pdb";
                $chLines = `wc -l < $chF`;
                next if ($AG eq $ch);
                
                if ( $chLines == $fLines) {
                    push (@tempArr, $ch);
                }
            }
        } 
        push (@antigenIds, @tempArr);
    }
    
##################

  # my @antigenIds =
   #     checkSingleChainRedundancy ($chainType_HRef, $chainIdChainTpye_HRef, $LOG, $ab);


    
    if ( $ab eq "L") {
        @singleChainAb = @{$chainType{Light}};
                
        @dimers = getDimers(\@PDBchains, \@antigenIds, $LOG);
        if ( @dimers) {
            mergeLH ($LOG, @dimers);
            print {$LOG} "DIMERS: @dimers\n";
        }
        else {
            @dimers = @singleChainAb;
            $dimer_flag = 0;
            print {$LOG} "NOT DIMERS:\n"; 
        }
    }
    elsif ( $ab eq "H") {
        @singleChainAb = @{$chainType{Heavy}};
    }
        
    my $fileType;
    my %fileTypeH;
    my $hapten;
    
    if ( $dimer_flag ) {
        $hapten = hasHapten ($pdbPath, \@dimers, $LOG);
        %fileTypeH = processHapten($pdbPath, \@dimers, $ab, $LOG, $dimer_flag );
        
    }
    # Check for haptens bound with CDRs
    else {
        $hapten = hasHaptenSingleChain ($pdbPath, \@singleChainAb, $chainIdChainTpye_HRef);
        %fileTypeH = processHapten($pdbPath, \@singleChainAb, $ab, $LOG);
    }

    if ( $dimer_flag ) {
        @singleChainAb = 0;
        @singleChainAb = @dimers;
    }
         
    # Checks for haptens and move them to non-Protein data 
    if ( ($hapten) and (!@antigenIds) )
    {

        makeFreeAntibodyComplex ($pdbId, $pdbPath, \@singleChainAb, $count, $fileType,
                                 $dir, $chainIdChainTpye_HRef, $numbering,
                                 $LOG, $destNonPro, $destFreeAb, $ab, %fileTypeH);
                
    }    

    # Checks for protein antigens (Further checking is within
    # processAntibodyAntigen subroutine)
    elsif ( (@antigenIds) and (!$hapten) )
    {
                
        $fileType = "num";
        $cdrError= processAntibodyAntigen($pdbId, $pdbPath, $ab, \@antigenIds,
                               \@singleChainAb, $fileType, $dir, $masterDir,
                               $LOG, $chainIdChainTpye_HRef, $destPro,
                               $destNonPro, $destFreeAb, $numbering, %fileTypeH);
        
    }
    # Antibody bound with antigen and haptens - Moved to protein antigen data
    elsif ( (@antigenIds) and ($hapten) )
    {
            
        $cdrError = processAntibodyAntigen($pdbId, $pdbPath, $ab, \@antigenIds,
                               \@singleChainAb, $fileType, $dir, $masterDir,
                               $LOG, $chainIdChainTpye_HRef, $destPro,
                               $destNonPro, $destFreeAb, $numbering, %fileTypeH);
    }
    # Checks for free light and heavy chains
    else {
        $fileType = "num";
        makeFreeAntibodyComplex($pdbId, $pdbPath, \@singleChainAb, $count,
                                $fileType, $dir, $chainIdChainTpye_HRef, $numbering,
                                $LOG, $destNonPro, $destFreeAb, $ab,  %fileTypeH);
    }
    return $numberingError;
    
}


sub getDimers
{
    my ($PDBchains_ARef, $antigenIds_ARef, $LOG) = @_;
    my @PDBchains = @{$PDBchains_ARef};
    my @antigenIds = @{$antigenIds_ARef};
    
    @PDBchains = sort @PDBchains;
    
    my @pos = ("36", "87", "40");
    my %dimers;
    my ($d36, $d40, $d87);
    
    for (my $i=0; $i < (scalar @PDBchains); $i++ )
    {
        # ignore antigen chains to check for contacts
        
        next if (grep {$_ eq $PDBchains[$i]} @antigenIds);        
                
        for ( my $j=($i+1); $j < (scalar @PDBchains); $j++ )
        {
            next if (grep {$_ eq $PDBchains[$j]} @antigenIds);
             
         foreach my $p (@pos)
         {
             # Look for numbered files
             my $f1 = $PDBchains[$i]."_num.pdb";
             my $f2 = $PDBchains[$j]."_num.pdb";
             my $fTemp = $PDBchains[$j].".temp";

             # Replace Chain label L with l to compute the distance
             open (my $IF2, '<', $f2) or die "Can not open file: $f2\n"; 
             my @a = <$IF2>;
             @a = map {s/\s+L\s+/ l /g; $_;} @a;
             close $IF2;

             # Writing chain l in temporary file
             open (my $OF2, '>', $fTemp) or die "Can not open file: $fTemp\n";
             print {$OF2} @a;

             
             my @dist = `cat $f1 $fTemp | pdbdist "L"$p "l"$p`;

             
             my @d = split (/ /, $dist[0]);
             if ($p == 40)
                 {
                     $d40=$d[3];
                 }
             elsif ($p == 36)
                 {
                     $d36=$d[3];
                 }
             elsif ($p == 87)
                 {
                     $d87=$d[3];
                 } 
         }

         print {$LOG} "$PDBchains[$i]$PDBchains[$j]\nD36:$d36\nD40:$d40\nD87:$d87\n";
         
             if ( ( $d36 < 20 ) and ($d87 < 22) and ($d40 <= 40) ) # Made D40 flexible
                 {
                     $dimers{"$PDBchains[$i]$PDBchains[$j]"}=1;
                 }
     }
          
    }
    my @dimers = keys (%dimers);
    return @dimers;
    
}

=pod
sub pairDimerChains
{
    my ($pdbPath, $LOG, $singleChainIDs_REF) = @_;
    my %chainContacts = getInterchainContacts($pdbPath);
    print {$LOG} "Antibody inter-chain contacts are: \n";
    print {$LOG} Dumper ( \%chainContacts );
    my %dimers;
    
    my @singleChainIDs = @{$singleChainIDs_REF};
    print "ANTIBODY: @singleChainIDs\n";
    foreach my $ab1 (@singleChainIDs ) {
        foreach my $ab2 (@singleChainIDs) {
            my $pair = $ab1.$ab2;
            if ( exists $chainContacts{$pair}) {
                if ( $chainContacts{$pair} > 80 ) {
                    $dimers{$pair} = $chainContacts{$pair};
                    
                }
            }
        }
    }
     # removing duplicate values
    %dimers = reverse %{ {reverse %dimers} };
    print {$LOG} "The dimerizing pairs are: \n", Dumper (\%dimers);
    return (%dimers);
}
=cut

sub hasHaptenSingleChain
{
    my ($pdbPath, $singleChainAb_ARef, $chainIdChainTpye_HRef) = @_;
    my $hashapten = $config::hashapten;

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
        
