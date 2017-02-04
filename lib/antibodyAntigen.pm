package antibodyAntigen;
use strict;
use Carp;
#use warnings;
use List::MoreUtils qw(uniq);
use config;
use Data::Dumper;
use IO::CaptureOutput qw ( capture capture_exec qxx qxy );

use antibodyProcessing qw (
	getChainTypeWithChainIDs	
	splitPdb2Chains
	hasHapten
	processHapten	
	movePDBs
	largestValueInHash
	extractCDRsAndFrameWorks
	checkAntigenChains
        mapChainsIDs
        printHeader
        mergeLH
                      );
use Exporter qw (import);
our @EXPORT_OK = qw (
	processAntibody
        antibodyNumbering
	makeFreeAntibodyComplex
        processAntibodyAntigen
        getInterchainContacts
        antibodyAssembly
                );
sub processAntibody
{
    my ($pdbId, $pdbPath, $nsch, $ab, $dir, $masterDir, $LOG, $numbering) = @_;
########
    my $destPro = "$masterDir"."/".$ab."_Protein_".$numbering;
    my $destNonPro = "$masterDir"."/".$ab."_NonProtein_".$numbering;
    my $destFreeAb = "$masterDir"."/".$ab."_Free_".$numbering;
    #my $numberingError = 0;  
    my ($chainType_HRef, $chainIdChainTpye_HRef) =
        getChainTypeWithChainIDs ($pdbPath);
   
    print {$LOG} "Antibody chain types with chain IDs: \n";
    print {$LOG} Dumper ($chainType_HRef);

    print {$LOG} "Antibody chain IDs, chain Labels and chain sequences: \n";
    print {$LOG} Dumper ($chainIdChainTpye_HRef);
    
    splitPdb2Chains($pdbPath);

        print {$LOG} "The $pdbId PDB has been splitted in to different files".
            " based on number of chains\n";
    
    my %chainType = %{$chainType_HRef};
    my $count = 1;

    # Obtaining only light and heavy chains labels from %chainType to number chains separately
     my @lChains=@{$chainType{'Light'}};
    my @hChains=@{$chainType{'Heavy'}};
        my @lhHybridChains=@{$chainType{'Light/Heavy'}};
    my @hlHybridChains=@{$chainType{'Heavy/Light'}};

    my @lhChains = (@lChains, @hChains, @lhHybridChains, @hlHybridChains);
    print {$LOG} "Light and heavy chains are @lhChains\n";
    #***


     my $numError = antibodyNumbering ( \@lhChains, $nsch );
    

    
    my ($heavyLightPairContact_HRef, $agContacts_HRef, $antigen_ARef) =
    pairHeavyLightChains ($pdbPath, $chainType_HRef, $chainIdChainTpye_HRef,
                          $LOG);
    my @antibodyPairs = keys % {$heavyLightPairContact_HRef};
    my @antigens = @{$antigen_ARef};# Informs about Free or complexed

    print {$LOG} "$pdbId has ", scalar @antibodyPairs, " pairs of antibodies\n";
    print {$LOG} join ("\n", @antibodyPairs ), "\n";
    print {$LOG} "$pdbId has ",scalar @antigens, " number of antigens\n";
    print {$LOG} join ("\n", @antigens ), "\n";
    
    #****
    # merge each antibody pair
    my $pairCount= scalar (@antibodyPairs);
    my $countFailedPair = mergeLH( $LOG, @antibodyPairs);
    
    my $numberingError=0;

    if ( ($numError) and ($countFailedPair >= $pairCount)) {
        $numberingError = 1;
    }
    
    my $hapten = hasHapten ($pdbPath, \@antibodyPairs);
    my $fileType;
    my %fileType;
    
    my$cdrError=0;
    
        
    # Checks for haptens and move them to non-Protein data
    if ( ( $hapten) and (!@antigens) )
    {
#        $fileType = "hap";
        %fileType = processHapten($pdbPath, \@antibodyPairs, $ab, $LOG);
        makeFreeAntibodyComplex($pdbId, $pdbPath, \@antibodyPairs, $count,
                                $fileType, $dir, $chainIdChainTpye_HRef, $numbering,
                                $LOG, $destNonPro, $destFreeAb, $ab, %fileType);
     
    }
    # Checks for protein antigens (Further checking is within
    # processAntibodyAntigen subroutine)
    elsif ( ( @antigens) and (!$hapten) )
    {
        $fileType = "num";
        $cdrError =
            processAntibodyAntigen($pdbId, $pdbPath, $ab, $antigen_ARef,
                               \@antibodyPairs, $fileType, $dir, $masterDir,
                               $LOG, $chainIdChainTpye_HRef, $destPro,
                               $destNonPro, $destFreeAb, $numbering, %fileType);
    }
    
    elsif ( ( @antigens) and ($hapten))
    {
#        $fileType = "hap";
        my %fileType = processHapten($pdbPath, \@antibodyPairs, $ab, $LOG);
        $cdrError =
            processAntibodyAntigen($pdbId, $pdbPath, $ab, $antigen_ARef,
                               \@antibodyPairs, $fileType, $dir, $masterDir,
                               $LOG, $chainIdChainTpye_HRef, $destPro,
                               $destNonPro, $destFreeAb, $numbering, %fileType);
    }

    # Checks for free antibodies and move then in Free antibodies data
    else
    {
        $fileType = "num";
        makeFreeAntibodyComplex($pdbId, $pdbPath, \@antibodyPairs, $count,
                                $fileType, $dir, $chainIdChainTpye_HRef, $numbering,
                                $LOG, $destNonPro, $destFreeAb, $ab, %fileType);

    }
    return $numberingError;
    
}



sub makeFreeAntibodyComplex
{
    my ($pdbId, $pdbPath, $antibodyPair_ARef, $count, $fileType,
        $dir, $chainIdChainTpye_HRef, $numbering, $LOG,
        $destNonPro, $destFreeAb, $ab, %fileTypeH) = @_;
    my @antibodyPairs = @ {$antibodyPair_ARef};
        
    my ($lookForFile, $newFile);
    
    foreach my $antibodyPair (@antibodyPairs)
    {
        
#        mergeLH($antibodyPair, $LOG);                 
        
        if ( %fileTypeH ) {
            $fileType = $fileTypeH{$antibodyPair};
        }
        else {
            
            $fileType = "num";
        }
        $lookForFile = $antibodyPair."_".$fileType.".pdb";
 #       $newFile = $pdbId."_".$count.".pdb";
        
        my $numLines =`cat $lookForFile | wc -l`;

#        # To check if there is one chain in the $antibodyPair.pdb 

        if ( ($numLines > 1000) and ($ab eq "LH") ){
                        
            $newFile = $pdbId."_".$count.".pdb";
            open (my $ABFILE, '>>',  "$dir/$newFile");
            my %mapedChains = mapChainsIDs ($antibodyPair,$chainIdChainTpye_HRef);
            printHeader($ABFILE, $numbering, $pdbPath, %mapedChains);
                        
            open (my $AB, '<', "$dir/$lookForFile") or die "Can't open File\n";
            while (!eof ($AB))
                {
                    my $freeAB = <$AB>;
                    next if ($freeAB =~ /^MASTER|^END/);
                    print {$ABFILE} $freeAB;
                }
        
            $count++;
        
            if ( $fileType eq "hap") {
                movePDBs ($dir, $destNonPro, $pdbId);
                print {$LOG} "This antibody is bound with hapten -- Moved to non- ".
                    "protein data antigen data\n";
            }
            else {
                movePDBs ($dir, $destFreeAb, $pdbId);
                print {$LOG} "This antibody is Free -- Moved to Free ".
                    "antibody data\n";
            }
            
        } # if numLines

        elsif ( ( ($ab eq "L") or ($ab eq "H") ) and ($numLines > 700) ) 
            {
                $newFile = $pdbId."_".$count.".pdb";    
                open (my $ABFILE, '>>',  "$dir/$newFile");
                my %mapedChains = mapChainsIDs ($antibodyPair,$chainIdChainTpye_HRef);
                printHeader($ABFILE, $numbering, $pdbPath, %mapedChains);
                        
            open (my $AB, '<', "$dir/$lookForFile") or die "Can't open File\n";
            while (!eof ($AB))
                {
                    my $freeAB = <$AB>;
                    next if ($freeAB =~ /^MASTER|^END/);
                    print {$ABFILE} $freeAB;
                }
        
            $count++;
        
            if ( $fileType eq "hap") {
                movePDBs ($dir, $destNonPro, $pdbId);
                print {$LOG} "This antibody is bound with hapten -- Moved to non- ".
                    "protein data antigen data\n";
            }
            else {
                movePDBs ($dir, $destFreeAb, $pdbId);
                print {$LOG} "This antibody is Free -- Moved to Free ".
                    "antibody data\n";
            }


                
            }# elsif
        
        
        
                
    }# $antibody pair
    
    
    return $count;
}
  
sub processAntibodyAntigen
{
    my ($pdbId, $pdbPath, $ab, $antigenIds_ARef, $antibodyPairs_ARef,
        $fileType, $dir, $masterDir, $LOG, $chainIdChainTpye_HRef,
        $destPro, $destNonPro, $destFreeAb, $numbering, %fileType) = @_;
    my @antibodyPairs = @{$antibodyPairs_ARef};

    my $cdrError = 0;    
    
    my @failedNumPair = extractCDRsAndFrameWorks ( \@antibodyPairs, $ab, $LOG );
   # print "FFFF: @failedNumPair \n";
    
    #eval { extractCDRsAndFrameWorks ( \@antibodyPairs, $ab, $LOG );
     #      1;
      # };
    
    #if ( $@ ) {
     #   print {$LOG} "CDR-Error: Can'nt extract CDRs for: $pdbId\n" . $@ .
      #      "Program exited\n\n";
       # $cdrError = 1;
       # return  $cdrError; 
        #next;
   # }
    #else {
     #   print {$LOG} "CDRs and FWs have been extracted from antibody for ".
      #      "contact analysis\n";
    #}
    
    
    my @antigenChains = @{$antigenIds_ARef};
    my %antibodyAntigenContacts =
        assembleCDRsAndFWsWithAntigen(\@antibodyPairs, \@antigenChains, \@failedNumPair);
    print {$LOG} "Antibody-Antigen contacts: \n";
    print {$LOG} Dumper (\%antibodyAntigenContacts);
    
    my %complexInfo = getComplexInfoInHash(%antibodyAntigenContacts);
    print {$LOG} "Antibody-antigen complex information by chain: \n";
    print {$LOG} Dumper (\%complexInfo);
    
    my $biAntigen = makeAntibodyAntigenComplex($pdbId, $pdbPath, $fileType,
                                               $dir, $masterDir, $LOG,
                                               $chainIdChainTpye_HRef,
                                               $destPro, $destNonPro, $destFreeAb,
                                               $numbering, \%complexInfo, \%fileType, $ab);
    return  $cdrError;
    
}


sub makeAntibodyAntigenComplex
    {
    my ( $pdbId, $pdbPath,$fileType, $dir, $masterDir, $LOG,
         $chainIdChainTpye_HRef, $destPro, $destNonPro, $destFreeAb,
         $numbering, $complexInfo_HRef, $fileType_HRef, $ab) = @_;
    #my $dir = '.';
    my $count = 1;
    my $biAntigen = 0;
    my $chaintype = $config::chaintype;
    my %complexInfo = %{$complexInfo_HRef};
    my %fileTypeH = %{$fileType_HRef};
    
    
    my @Freeantibodychains =
        grep {$complexInfo{$_} eq "NULL" } keys %complexInfo;

    if ( @Freeantibodychains ) {
        $fileType = "num";
               
        $count = makeFreeAntibodyComplex($pdbId, $pdbPath,
                                         \@Freeantibodychains, $count,
                                         $fileType, $dir, $chainIdChainTpye_HRef,
                                         $numbering, $LOG, $destNonPro, $destFreeAb, $ab,
                                         %fileTypeH);
    
    }
    
    my @antigen;
    my @antigenChainType;
    
    foreach my $ab_pair( keys %complexInfo )
    {
        $fileType = $fileTypeH{$ab_pair};
        
        my $numberedAntibody;
        
        if ( $fileType eq "num" )
        {
            $numberedAntibody = $ab_pair."_num.pdb";
        }
        elsif ( $fileType eq "hap" )
        {
            $numberedAntibody = $ab_pair."_hap.pdb";
        }
        else {
            
            $numberedAntibody = $ab_pair."_num.pdb";
        }
            my $antigenRef;
        
        next if ( $complexInfo{ $ab_pair } eq  "NULL" );
        if ( $complexInfo{ $ab_pair } ne  "NULL" )
        {
            $antigenRef = $complexInfo{$ab_pair};
            @antigen = @{$antigenRef};
        }
        my $numLines =`cat $numberedAntibody | wc -l`;
        if ( $numLines > 900)
        {
        open ( my $AB_FILE, '<', "$dir/$numberedAntibody" ) or
            die "Could not open file $numberedAntibody";
        
        open ( my $AG_AB_FILE, '>>',"$pdbId"."_"."$count".".pdb" ) or
            die "Can not write complex";

        # Print headers
        my %mapedChains = mapChainsIDs ($ab_pair,$chainIdChainTpye_HRef, %complexInfo );     
        printHeader($AG_AB_FILE, $numbering, $pdbPath, %mapedChains);
                
        # In case of multiple antigens making contacts with antibody
        
        while (!eof ( $AB_FILE ) )
        {
            my $antibody = <$AB_FILE>;
            next if ($antibody =~ /^MASTER|^END/);
            print $AG_AB_FILE $antibody;
        }
        # In case of multiple antigens making contacts with antibody 
        foreach my $agChain(@antigen)
        {
            $agChain = $agChain.".pdb";
            open ( my $AG_FILE , '<', "$dir/$agChain" ) or
                die "Could no open file $agChain";
            
            while ( !eof ( $AG_FILE ))
            {
                my $antigen = <$AG_FILE>;
                # To ignore the Footer of PDB
                next if ($antigen =~ /^MASTER|^END/);
                print $AG_AB_FILE $antigen;
            }
            # To obtain the antigen chain label and chain type (N Protein)
            @antigenChainType = `$chaintype -c $agChain $pdbPath`;
        }

        
            if ( grep {/DNA|RNA/} @antigenChainType)
            {
                # Move to Non-Protein Dir
                movePDBs ($dir, $destNonPro, $pdbId);
                print {$LOG} "The $pdbId has non-protein antigen -- ".
                    "Moved to Non-protein antigen data\n";
            }
            else {
                # Move to protein Dir
                movePDBs ($dir, $destPro, $pdbId);
                print {$LOG} "The $pdbId has protein antigen -- ".
                    "Moved to Protein antigen data\n";
            }
            
        $count++;    
    }

    }
    
    
    if ( (scalar @antigen) >= 2 )
    {
        $biAntigen = 1;
        print {$LOG} "$pdbId: Antibody is bound with 2 or more antigens (multi-chain antigen)\n";
    }

    return $biAntigen;
}

# ************* getComplexInfoInHash *****************
# Description: Computes antibody and antigen complex on basis of contact info
#              of antigen (based on CDRs) with antibody. There are certain PDBs
#              where 2 antigens binds with single antibody. Here such antigens
#              are paired with antibodies by calculating highest and second
#              highest number of contacts
# Inputs: A hash containing antibody and antigen contact information
# Outputs: A hash with antibody and antigen chain IDs information for correct
#          pair
# Subroutine call/Testing: my %complex =
#                               getComplexInfoInHash(%antibodyAntigenContacts);
# Other subroutine calls: None
# Date: 26 June 2014
# Author: Saba 
sub getComplexInfoInHash
{
    my ( %antibodyAntigenContacts ) = @_;
    my %complex;
    my %temp;
    
    # Contacts information is used to form antibody and antigen complexes
    foreach my $ab_pair (keys %antibodyAntigenContacts)
    {
        my %antigenChains = %{ $antibodyAntigenContacts{$ab_pair} };
        %temp = ();      
        foreach my $key (keys %antigenChains)
        {
            # Storing non zero antigen contacts keys (chain Ids)
            if ( $antigenChains{$key} > 0 )
            {
                $temp{$key} = $antigenChains{$key};
            }

        }
        if ( !%temp ) {
            $complex{$ab_pair} = "NULL";
        }
        else {
            $complex{$ab_pair} = [keys %temp];
        }
        
    }
    
    return %complex;
}
   
# ************* pairHeavyLightChains *****************
# Description: This sub-routines does multiple jobs:
#              1. returns light and heavy pairs with inter-chain contacts
#                 in a hash
#              2. returns antigen inter-chain contacts in a hash
#              3. returns antigen chain IDs
# Inputs: PDB file path
# Outputs: 2 hashes - 1) antibody pair (light and heavy chains) with contacts
#                     2) antigen-antigen pair with contacts
#          1 array - containg antigen chain IDs
# Subroutine call/Testing: my ( $heavyLightPairContact_HRef, $agContacts_HRef,
#                          $antigenIds_ARef) =
#                       &pairHeavyLightChains ( $pdbPath, % {$chainType_HRef});
# Other subroutine calls: 1) checkAntigenChains
#                         2) getInterchainContacts
#                         3) largestValueInHash
# Date: 26 June 2014
# Author: Saba
sub pairHeavyLightChains
{
    my ($pdbPath, $chainType_HRef, $chainIdChainTpye_HRef, $LOG) = @_;
    my %chainType = %{$chainType_HRef};
    
    my  @light =  @{$chainType{Light}};
    my  @heavy =  @{$chainType{Heavy}};
    my  @antigen = checkAntigenChains($LOG, %chainType);
    
    my (%heavyLightPairContact, %agContacts);

    # using contact hash (hash with inter-chain contact information)
    my %chainContacts = getInterchainContacts($pdbPath);
    print {$LOG} "Antibody inter-chain contacts are: \n";
    print {$LOG} Dumper ( \%chainContacts );
    
    my %antibodyChains = ();
    
    # find light and heavy pair with maximum contacts
    foreach my $lg ( @light )
    {
        foreach my $hv ( @heavy )
        {
            my $abPair = $lg.$hv;
            if ( exists $chainContacts{$abPair} )
            {
                $antibodyChains{$abPair} = $chainContacts{$abPair};
            }
            elsif ( !exists $chainContacts{$abPair} )
            {
                $antibodyChains{$abPair} = 0;
            }
        }
        # Select key with largest value
        my ( $key, $val ) = largestValueInHash( \%antibodyChains);
        $heavyLightPairContact{$key}= $val;
        %antibodyChains = (); # Empty hash
    }
    # Find antigen-antigen chain contact pairs
    foreach my $ag1( @antigen )
    {
        foreach my $ag2( @antigen )
        {
            my $agPair = $ag1.$ag2;
            if ( exists $chainContacts{$agPair} )
                {
                    $agContacts{$agPair} = $chainContacts{$agPair};
                }
        }
    }
    # removing duplicate values
    %agContacts = reverse %{ {reverse %agContacts} };

    my ($heavyLightPairContact_HRef2, $antigenIds_ARef2) =
        checkChainRedundancy (\%heavyLightPairContact, $chainIdChainTpye_HRef,
                              $LOG, \@antigen);
    
    
    return ( $heavyLightPairContact_HRef2, \%agContacts, $antigenIds_ARef2);
}

# ************* getInterchainContacts *****************
# Description: Calculates the inter-chain contacts and return them into a hash
# Inputs: PDB file path
# Outputs: Returns a hash with chain pair as key and number of contacts as val
# Subroutine call/Testing: my %chainContacts = getInterchainContacts($pdbPath);
# Date: 26 June 2014
# Author: Saba

sub getInterchainContacts
{
    my ( $pdbPath) = @_;
    # The program chaincontacts calculates inter chain contacts
    my $chaincontacts = "$config::chaincontacts -r 4.00";
    my @chainContacts = `$chaincontacts $pdbPath`;

    # To remove the headers from the output of chaincontacts program
    splice @chainContacts, 0, 8;
    my %chainTable=();
    my $ncontacts=0;
    
    foreach my $line ( @chainContacts )
    {
        chomp ($line);
        # Chain: A Res:  79  - Chain: H Res: 105  Contacts: 13
        
        if ($line =~ /^Chain*/)
        {
            my @cont= split /:/, $line;
            my $n = $cont[5];
            chomp($n);
            my $ch1 = $cont[1];
            my $ch2 = $cont[3];
            chomp($ch1);
            chomp($ch2);
            $ch1 =~ s/Res//; # Convert A Res to A 
            $ch2 =~ s/Res//;
            $ch1 =~ s/^\s*(.*?)\s*$/$1/; # Removing leading and trailing space
            $ch2 =~ s/^\s*(.*?)\s*$/$1/;
            # Create Chain's pair e.g: LH
            my $chkey=$ch1.$ch2;
            # Counting number of contacts between 2 chains and storing them
            # in a hash like LH 150
            if ( exists $chainTable{$chkey} )
                {
                    $ncontacts = $chainTable{$chkey} + $n;
                    $chainTable{$chkey}=$ncontacts;
                }
            else
                {
                    $ncontacts=$n;
                    $chainTable{$chkey}=$ncontacts;
                }
        }
    }

    return %chainTable;

}

sub checkChainRedundancy
{
    my ($heavyLightPairContact_HRef, $chainIdChainTpye_HRef, $LOG,
        $antigen_ARef) = @_;

    my %chainIdChainTpye = %{$chainIdChainTpye_HRef};
    my @antigens = @{$antigen_ARef};
    my %heavyLightPairContact = %{$heavyLightPairContact_HRef};    
    my @abPairs = sort keys (%{$heavyLightPairContact_HRef});

    
    for (my $i = 0 ; $i<= $#abPairs-1; $i++ )
    {
        my ($l1, $h1) = split ("", $abPairs[$i]);
        my $abConts = $heavyLightPairContact{$abPairs[$i]};
        
        for (my $j = $i+1; $j<= $#abPairs; ++$j )
        {
            my ($l2, $h2) =split ("", $abPairs[$j]);
            # Check pairwise redundancy by using index function
            # which check if 2nd string is substring of 1st string
            # If an antibody is not redundant to the other then treat that
            # antibody chain (L and H) as antigens

            if ( ( index ($chainIdChainTpye{$l1}{"L"},
                          $chainIdChainTpye{$l2}{"L"}) == -1 )
                     and ( index ($chainIdChainTpye{$h1}{"H"},
                                  $chainIdChainTpye{$h2}{"H"}) == -1 ) )
                {
                    push (@antigens, $l2, $h2, $l1, $h1);
                    # Treat antibody chains as antigens
                    # Also keep antibody pair in hash
                    $heavyLightPairContact{$abPairs[$i]} = $abConts;
                    print {$LOG} $abPairs[$j]." is ". 
                        "found as non-redundant to other antibodies in PDB ".
                            "and will be treated as idiotypic antigen\n";
                }
            else {
                next;
            }

        }
    }
    @antigens = uniq (@antigens);
    
    return (\%heavyLightPairContact, \@antigens);
}

# ************* antibodyAssembly *****************
# Description: It assembles two files (L and H) of an antibody together for a
#              given chains pair.
#              In our case it assembles L and H chain of an antibody in 1 PDB
# Inputs: An array contaning chain label pairs
# Outputs: It returns number of antibodies assembled and writes separate file
#          for each antibody
# Subroutine call/Testing: my $count = antibodyAssembly ( @antibodyPairs );
# Other subroutine calls: 
# Date: 26 June 2014
# Author: Saba
sub antibodyAssembly
{
    my ( @antibodyPairs ) = @_;
    my $dir = '.';
    my $count = 0;
    my ($L_READ, $H_READ, $OUT);
    foreach my $pair ( @antibodyPairs )
    {
        my @chain = ();
        my ( $L, $H ) = split ( '', $pair );
        my $l_pdb = $L.".pdb";
        my $h_pdb = $H.".pdb";
        open ( $L_READ, '<', "$dir/$l_pdb" ) or
            die "Can not open file $l_pdb\n";
        open ( $H_READ, '<', "$dir/$h_pdb" ) or
            die "Can not open file $h_pdb\n";
        open ( $OUT, '>', "$dir/$pair.pdb" ) or
            die "Can not write file $pair.pdb\n";

        while ( !eof ($L_READ ) )
        {
            my $light = <$L_READ>;
            print {$OUT} $light;
        }

        while ( !eof ($H_READ ) )
        {
            my $heavy = <$H_READ>;
            print {$OUT} $heavy;
        }
        $count++;
    }
    close $L_READ;
    close $H_READ;
    close $OUT;

    return $count;
}
# ************* antibodyNumbering *****************
# Description: Numbers the antibody according to the user defined numbering
#              scheme
# Inputs: A reference to array containing antibody (LH) chain labels and
#         numbering scheme
# Outputs: It returns number of antibodies numbered and writes separate file
#          for each numbered antibody
# Subroutine call/Testing: my $count1 = antibody_number (\@keys, $nsch);
# Date: 26 June 2014

sub antibodyNumbering
{
    my ( $antibodyChains, $nsch ) = @_;
    my ( $out, $error, $numberedPDB );
    open (my $OUT, '>', "numberingFailedChain.dat");
    my $fail;    
    foreach my $c( @$antibodyChains )
    {
        my $antibody = $c.".pdb";
        my $kabatnum = $config::kabatnum;
        $numberedPDB = $c."_num.pdb";

        # Capturing error from Stderr
        # In case of error (kabat program failure), output file ($numbered_pdb)
        # would be empty
               
        ($out, $error) = qxx ( "$kabatnum $nsch $antibody $numberedPDB" );

        if ( -z $numberedPDB )
            { # Checks if file is empty
                print {$OUT} "$numberedPDB failed to number\n";
                print {$OUT} "$error\n";
                $fail = 1;
                
                #croak $error;
            }
    }
    return $fail;
}


# ************* assembleCDRsAndFWsWithAntigen *****************
# Description: 1) Assembles CDR and frame work regions with antigen and writes
#              them into a separate PDB file in order to calculate the number
#              of contacts of each antibody (CRDs) with antigen
#              2) Number of CRD and FW contacts with antigen are computed and
#              used to find the pairing of antigen with antibody
# Inputs: 2 array references containing antibody (LH) chain labels and
#         antigen chain labels
# Outputs:  Returns a hash of hash containg antibody contacts with
#           correnponding antigen (computed based on number of contacts)
#           between CDRs/FWs and antigen
# Subroutine call/Testing: my %antibodyAntigenContacts =
#                             assembleCDRsAndFWsWithAntigen
#                                  (\@antibodyPairs, \@antigenChains);
# Other subroutine calls: antigenAntibodyContacts 
# Date: 26 June 2014
# Author: Saba
sub assembleCDRsAndFWsWithAntigen
{
    my ( $antibodyPairs, $antigenChains, $failedNumPair) = @_;
    my $dir = '.';
    my ($AG_FILE, $CDR_FILE, $AG_CDR_FILE, $cdr_ag_conts );
    my ($FW_FILE, $AG_FW_FILE, $fw_ag_conts);
    
    my %antibodyAntigenContacts = ();
    my $nonAgCount;
   
    
    foreach my $antibody ( @$antibodyPairs )
    {

        if ( $antibody ~~ @{$failedNumPair}) {
            next;
        }
        my $cdrs = $antibody."_CDR.pdb";
        my $fw = $antibody."_FW.pdb";
        
        $antibodyAntigenContacts{$antibody} = {};
        
        foreach my $agChain ( @$antigenChains )
        {
            # If non-redundant antibody is treated as antigen and it has L and H
            # chain labels then ignore them and do not assemble them with
            # antibody chain
            next if ( ( $agChain eq "L") or ( $agChain eq "H") );
            
            my $agPDB = $agChain.".pdb";
            #  next if ( $agPDB eq '0.pdb');# to deal exception in 3J30 (split prob)
            # Open antigen PDB file
            open ( $AG_FILE, '<', "$dir/$agPDB" ) or
                die "Could no open file $agPDB\n";
            # Open CDRs PDB file
            open ( $CDR_FILE, '<', "$dir/$cdrs" ) or
                die "Could not open file $cdrs\n";
            # Open new file to write complex of antigen and CDR co-ordinates
            open ( $AG_CDR_FILE, '>',"$agChain"."_"."$cdrs" ) or
                die "Can not write file $cdrs\n";

            # Open FWs PDB file
            open ( $FW_FILE, '<', "$dir/$fw" ) or
                die "Could not open file $fw\n";
            # Open new file to write complex of antigen and FW co-ordinates 
            open ( $AG_FW_FILE, '>',"$agChain"."_"."$fw" ) or
                    die "Can not write file $fw\n";
            
            while ( !eof ( $AG_FILE ) )
            {
                my $antigen_pdb = <$AG_FILE>;
                print $AG_CDR_FILE $antigen_pdb;
            }

            while ( !eof ($CDR_FILE ) )
            {
                my $cdrs_pdb = <$CDR_FILE>;
                print $AG_CDR_FILE $cdrs_pdb;
            }

            close $AG_FILE;
            # Re-open antigen file
            open ( $AG_FILE, '<', "$dir/$agPDB" ) or
                    die "Could no open file $agPDB\n";
            while ( !eof ( $AG_FILE ) )
            {
                my $antigen_pdb = <$AG_FILE>;
                print $AG_FW_FILE $antigen_pdb;
            }
            
            while ( !eof ($FW_FILE ) )
            {
                my $fw_pdb = <$FW_FILE>;
                print $AG_FW_FILE $fw_pdb;
            }
            my $ag_cdr_filename = "$agChain"."_"."$cdrs";
            my $ag_fw_filename = "$agChain"."_"."$fw";
                
            # To calculate the contacts of antigen with CDRs and FWs regions
            # of an antibody
            ($cdr_ag_conts, $fw_ag_conts) =
                antigenAntibodyContacts ( $ag_cdr_filename, $ag_fw_filename,
                                          $agChain );
            # Checks if CDR contacts with antigen are more than FW contacts
            # then forms complex of antigen with that antibody chains pair
            if ( $cdr_ag_conts > $fw_ag_conts)
            {
                # To ensure that antibody chain (as antigen ) does not check
                # contacts to itself 
                if ( $antibody =~ /$agChain/) {
                    $antibodyAntigenContacts{$antibody}->{$agChain} = 0;
                }
                else {
                    # To drop antigens making 2-15 contacts with CDRs
                    # They are most likely X-ray crytal of antigen bound
                    # with another antibody
                    if ( $cdr_ag_conts >= 15 ) {
                        $antibodyAntigenContacts{$antibody}->{$agChain}
                            = $cdr_ag_conts;
                    }
                    else {
                        $antibodyAntigenContacts{$antibody}->{$agChain} = 0;
                    }
                }
                
            }    
            else {
                $antibodyAntigenContacts{$antibody}->{$agChain} = 0;
            }
        }
    }
    close ($AG_FILE);
    close ($CDR_FILE);
    close ($FW_FILE);
    close ($AG_CDR_FILE);
    close ($AG_FW_FILE);
    return  %antibodyAntigenContacts;
}


# ************* antigenAntibodyContacts *****************
# Description: Computes number of contacts between CDRs/FWs and antigen
# Inputs: 2 PDB files (Ag and CDs complex, Ag and FWs complex) and Ag chain ID
# Outputs: Number of contacts CDRs and FWs making with antigens
# Subroutine call/Testing: antigenAntibodyContacts
#                                ( $ag_cdr_filename, $ag_fw_filename,$agChain);
# Other subroutine calls: countCDRContacts
# Date: 26 June 2014
# Author: Saba

sub antigenAntibodyContacts
{
    my ($AgCDRsPDB, $AgFWsPDB, $agChain ) = @_;
    my $CDRsContacts =
        "$config::chaincontacts -r 4.00 -x LH -y $agChain";
    my @CDRsContacts = `$CDRsContacts $AgCDRsPDB`;

    my $FWsContacts =
      "$config::chaincontacts -r 4.00 -x LH -y $agChain";
    my @FWsContacts = `$FWsContacts $AgFWsPDB`;

    my $n;
    splice @CDRsContacts, 0, 8;
    splice @FWsContacts, 0, 8;
    
    my $nCDRs = countCDRContacts (@CDRsContacts);
    my $nFWs = countCDRContacts (@FWsContacts);

    return ($nCDRs, $nFWs);
}

# ************* countCDRContacts *****************
# Description: Counts CDR contacts from given array (contains output from
#              external program - chaincontacts)
# Inputs: An array containing contacts information 
# Outputs: Total number of contacts
# Subroutine call/Testing: countCDRContacts (@CDRsContacts);
# Other subroutine calls: None
# Date: 26 June 2014
# Author: Saba

sub countCDRContacts
{
    my (@cdr_conts) = @_;
    my $n = 0;
    foreach my $line ( @cdr_conts )
    {
        chomp ($line);
        if ( $line =~ /^Chain*/ )
        {
            my @conts = split /:/, $line;
            $n += $conts[5];
            chomp($n);
        }
    }
    return $n;
}


