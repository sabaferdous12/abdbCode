package antibodyProcessing;

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw (uniq);
use SFPerlVars;
use Carp;
use File::Copy;
use general qw (readDirPDB);
use IO::CaptureOutput qw ( capture capture_exec qxx qxy );
use Exporter qw (import);

our @EXPORT_OK = qw (
                        getPDBPath 
                        readDirPDB
                        getChains
			getChainTypeWithChainIDs
			splitPdb2Chains
                        extractCDRsAndFrameWorks
                        readFileDataInArray
                        getSEQRESFromPDB
			Aa3to1 
                        dirOperations
			hasHapten
			processHapten
			movePDBs
			largestValueInHash
			checkAntigenChains	
                        mapChainsIDs
                        printHeader
                        getProcessedABchains
                );
# ************* getPDBPath *****************
# Description: For given pdb it return local pdb file path
# Inputs: PDB code
# Outputs: Complete file path for local pdb e.g: /acrm/data/pdb/pdb3u7a.ent
# Subroutine call/Testing: $filePath = getPDBPath($pdb_id)
# Date: 26 June 2014
# Author: Saba 
sub getPDBPath
{
    my ( $pdbId ) = @_;
    chomp $pdbId;
    my $filePath;
    
    if ( $pdbId=~m/[0-9A-Z]{4}/ )
    {
        my $pdbprep = $SFPerlVars::pdbprep;
        my $pdbext = $SFPerlVars::pdbext;
        $filePath = $pdbprep . lc ( $pdbId ) . $pdbext;
    }
    else
    {
        croak "This is not a valid PDB Code\n";
    }

    return $filePath;
}

# **************** getChainTypeWithChainIDs *************
# Description: For given antibody, it determines chain types and chain IDs and
#              sequence from SEQRES records of a PDB and place them in a hash. 
# Inputs: PDB Path
# Outputs: 2 hashes; 1) hash containing array reference where hash key is chain
#          type and value is array reference, array contains chains IDs for the
#          chain type. 2) hash of hashes containing chain IDs as keys and chain#          label (L for light, H for heavy and A for antigen) and sequence for
#          each of the chain as values.
# Subroutine call/Testing: my ($chainType_HRef, $chainIdChainTpye_HRef) =
#                             getChainTypeWithChainIDs($pdbPath);
# Other subroutine Calls: 1) getSEQRESFromPDB
# Date: 19 Nov 2015
# Author: Saba
sub getChainTypeWithChainIDs
{
    my ($pdbPath) = @_;
    my $idabchain = $SFPerlVars::idabchain;
    # Gets chain IDs and chain Types from program idabchain
    my @chainInfo = `$idabchain $pdbPath`;
    
    my %chainType = ();
    my %chainIdChainTpye;
    my $chainResSeq;

    $chainType{'Antigen'} = [];
    $chainType{'Light'} = [];
    $chainType{'Heavy'} = [];
    
    foreach my $line ( @chainInfo )
        {
            #Chain 1, A: Antigen
            my @a = split(' ', $line);
            my $chainid = $a[2];
            # Retaing only chain ID e.g 'A' removing any leading char
            ($chainid) = $chainid =~ m/([A-Za-z0-9]{1})/;
            my $chainType = $a[3];
            
            if ( $chainType eq "Antigen" )
                {    
                    push ( @{ $chainType{'Antigen'}}, $chainid );
                    # Get SEQRES sequence from PDB file (SEQRES records)
                    $chainResSeq = getSEQRESFromPDB($pdbPath, $chainid);
                    $chainIdChainTpye{$chainid}{"A"} = $chainResSeq; 
                }
            elsif ( $chainType eq "Light" )
                {
                    push ( @{ $chainType{'Light'}}, $chainid );
                    # Get SEQRES sequence from PDB file (SEQRES records)
                    $chainResSeq = getSEQRESFromPDB($pdbPath, $chainid);
                    $chainIdChainTpye{$chainid}{"L"}= $chainResSeq;
                }
            elsif ( $chainType eq "Heavy" )
                {
                    push ( @{ $chainType{'Heavy'}}, $chainid );
                    # Get SEQRES sequence from PDB file (SEQRES records)   
                    $chainResSeq = getSEQRESFromPDB($pdbPath, $chainid);
                    $chainIdChainTpye{$chainid}{"H"}= $chainResSeq;
                }
        }
    
    return (\%chainType, \%chainIdChainTpye);
}



sub getChains
{
    my ($pdbPath) = @_;
    my ($chainType_HRef, $chainIdChainTpye_HRef) =
        getChainTypeWithChainIDs ($pdbPath);
    my %chainType = %{$chainType_HRef};

    my (@antigen, @light, @heavy);
    @antigen = @{$chainType{'Antigen'}};
    @light = @{$chainType{'Light'}};
    @heavy = @{$chainType{'Heavy'}};

    return (\@light, \@heavy, \@antigen);
}

# ************* splitPdb2Chains *****************
# Description: Splits PDB into chains and keep in separate files
# Inputs: PDB file path
# Outputs: None
# Subroutine call/Testing: $chainno = splitPdb2Chains($pdbPath);
# Other subroutine calls: 1) readDirPDB($dir)
# Date: 26 June 2014
# Author: Saba

sub splitPdb2Chains
{
    my $dir = ".";
    my ( $pdbPath ) = @_;
    open ( my $PDBFILE, $pdbPath) or die "Could not open file \"$pdbPath\"" ;
    # This program splits chains in a single PDB file into multiple PDBs
    my $pdbSplitChains = $SFPerlVars::pdbsplitchains;
    my $pdbatoms = $SFPerlVars::pdbatoms;
    my $pdbgetchain = $SFPerlVars::pdbgetchain;
        
    # To strip header and keeping only ATOMs and then spliting PDB into chains
    `$pdbatoms $pdbPath | $pdbSplitChains`;
    my @PDBchains = readDirPDB($dir); # Reads all PDB files in directory
    
    foreach my $chain (@PDBchains)
        {
            my ($chainID, $ext) = split (/\./, $chain);
            # To discard HETATM by selecting each chain and keeping
            # only ATOM lines                                       
            `$pdbgetchain -a $chainID $chain >$chainID`;
            `mv $chainID $chainID.pdb`; # Rename A as A.pdb
        }
    close $PDBFILE;
}

# ************* checkAntigenChains *****************
# Description: checks antigen array for L and H chain IDs and rename them as l
#              and h
# Some of antigens have chain labels of H and L (Same as antibody chains).
# These need to be renamed. The reason is to avoid the problem that arises
#  when grouping light and heavy LH with antigen to compute number of contacts
# in CDRs of antibody with antigen.
# e.g: In cases with antigen labelled as H will give wrong results. LH->H
# Therefore, this bit of code is renaming the antigen by converting the
# chain label to lower case and renaming the PDB file as %h or %l
# and updating the antigen labels array

# Inputs: chain Type hash
# Outputs: Array with chain Ids (new chain IDs if antigen has L or H chain IDs)
# Subroutine call/Testing: my  @antigen = checkAntigenChains(%chainType);
# Other subroutine calls: None
# Date: 26 June 2014
# Author: Saba

sub checkAntigenChains
{
    my ($LOG, %chainType) = @_;
    # Extract antigen chains from chain Type hash
    my @antigen = @{$chainType{Antigen}};
    my $renumpdb = $SFPerlVars::renumpdb;
    my $index=0;
    my $flag = 0;
    
    foreach my $antigenChain( @antigen )
    {
        chomp $antigenChain;
        # If antigen has L or H chain IDs, rename them as l or h and
        # save it as %l and %h 
        if( ( $antigenChain eq "L" ) or
                ( $antigenChain eq "H" ) )
        {
            my $antigenPDB = $antigenChain.".pdb"; # L as L.pdb
            my $newChainLabel = lc ($antigenChain);
            # Rename L chain ID as l and save as %l.pdb 
        `$renumpdb -c $newChainLabel -n -d $antigenPDB "%"$newChainLabel.pdb`;
            my $filename = "%".$newChainLabel;
            rename ( $antigenPDB, "Backup_$antigenPDB" );#L.pdb as BackUp_L.pdb
            # Replacing H or L chain IDs by %h or %l in antigen array  
            splice (@antigen, $index, 1, $filename);
            print {$LOG} "Antigen chain $antigenChain has been renumbered ".
                "for chain ID as $newChainLabel and file PDB file for chain ".
                    "has been renamed as %".$newChainLabel.".pdb\n";
            $flag++;
        }
        $index++;

}
    if ( $flag ) {
        print {$LOG} "New antigens ID data is: \n";
        print {$LOG} join ("\n", @antigen), "\n";
    }
    
    return @antigen;                
}




# ************* checkChainRedundancy *****************
# Description: Redundancy among antibody pairs is checked and if one antibody
#              is redundant to the other then that antibody will be treated as
#              antigen and antigen chain IDs array will be updated
# Inputs: Requires 2 heavyLightPair hash and chain type hash and antigen IDs
#         array - all of these are passed as references
# Outputs: Returns hash with antibody chain pairs and antigen IDs array 
# Subroutine call/Testing: ($heavyLightPairContact_HRef, $antigenIds_ARef )
#                         = checkChainRedundancy($heavyLightPairContact_HRef,
#                            $chainIdChainTpye_HRef, $antigenIds_ARef);
# Other subroutine calls: None
# Date: 26 June 2014
# Author: Saba



# ************* largestValueInHash *****************
# Description: Finds the largest value in a hash
# Inputs: Takes hash as input
# Outputs: returns key with largest value
# Subroutine call/Testing: my ($key, $val) = largest_value_hash(\%hash);
# Date: 26 June 2014
# Author: Saba 
sub largestValueInHash
{
    my $hashRef   = shift;
    my ( $key, @keys ) = keys   %$hashRef;
    my ( $large, @vals ) = values %$hashRef;
    
    for ( 0 .. $#keys )
    {
        if ( $vals[$_] > $large )
            {
                $large = $vals[$_];
                $key = $keys[$_];
            }
    }
    return $key, $large;
}


# ************* extractCDRsAndFrameWorks *****************
# Description: Extracts CDR and framework regions from antibody according
#              to kabat definition
# Inputs: A reference to array containing antibody (LH) chain labels and
#         complex type
# Outputs: Write 2 PDB files containing CDR and frame work co-ordinates
# Subroutine call/Testing: extractCDRsAndFrameWorks ( \@antibodyPairs, $ab )
# Other subroutine calls: None
# Date: 26 June 2014
# Author: Saba

sub extractCDRsAndFrameWorks
{
    my ( $antibodyPairs, $ab) = @_;
    my ( $out, $error, $success, $CDR );
    
    foreach my $pair ( @$antibodyPairs )
    {
        my $numberedAntibody = $pair."_num.pdb";
        my $getpdb= $SFPerlVars::getpdb;
        my $CDRsPDB = $pair."_CDR.pdb";
        my $FWsPDB = $pair."_FW.pdb";
        
        my @numbering;
        
        if ( $ab eq "LH")
        {
            @numbering =  (["L24", "L34"], ["L50", "L56"], ["L89", "L97"],
                           ["H31", "H35"], ["H50", "H65"], ["H95", "H102"] );
                    # Extract Frame works
`$getpdb -v L24 L34 $numberedAntibody |
$getpdb -v L50 L56 |
$getpdb -v L89 L97 |
$getpdb -v H31 H35 |
$getpdb -v H50 H65 |
$getpdb -v H95 H102 >>$FWsPDB`;
        }
        
        elsif ( $ab eq "L")
        {
            @numbering =  (["L24", "L34"], ["L50", "L56"], ["L89", "L97"]);
            `$getpdb -v L24 L34 $numberedAntibody |
$getpdb -v L50 L56 |
$getpdb -v L89 L97 >>$FWsPDB`;
        }
        
        elsif ( $ab eq "H")
        {
            @numbering =  (["H31", "H35"], ["H50", "H65"], ["H95", "H102"] );
            `$getpdb -v H31 H35 $numberedAntibody | 
$getpdb -v H50 H65 |
$getpdb -v H95 H102 >>$FWsPDB`;
        }
        my $all_error = '';
        # to avoid append and clearing the file contents if already present
        if (-e $CDRsPDB)
        {
            open ($CDR, '>', $CDRsPDB) or die "can not open $CDRsPDB";
        }

        # Extract CDRs
        foreach my $cdr (@numbering)
        {
            my ( $start, $end ) = @{$cdr};
            ($out, $error) =
                qxx ( "$getpdb $start $end $numberedAntibody >>$CDRsPDB" );
            if ( $error )
            {
                $all_error .= $error;
            }
            
        }
        # Extract Frame works
#`$getpdb -v L24 L34 $numberedAntibody |
#$getpdb -v L50 L56 |
#$getpdb -v L89 L97 |
#$getpdb -v H31 H35 |
#$getpdb -v H50 H65 |
#$getpdb -v H95 H102 >>$FWsPDB`;

        if ( $all_error )
            {
                croak $all_error;
            }
    }
}



# ************* getSEQRESFromPDB *****************
# Description: For given PDB file (path) and chain ID, it extracts SEQRES
#              records from original PDB file
# Inputs: PDB file path and chain ID
# Outputs: String of sequence
# Subroutine call/Testing:  $chainResSeq =
#                                    getSEQRESFromPDB($pdbPath, $chainid);
# Other subroutine calls: 1) readFileDataInArray
#                         2) aa3to1
# Date: 26 June 2014
# Author: Saba

sub getSEQRESFromPDB
{
    my ($pdbPath, $chainID) = @_;
    my @pdbFile = readFileDataInArray ($pdbPath);
    my @seqResLines;
    my @seqRes3;
    foreach my $line ( @pdbFile )
    {
        if ( ($line =~ /^SEQRES/) and ( $line =~ m/\s+$chainID\s+/) )
        {
            @seqResLines = split (/\s+/, $line);
            splice ( @seqResLines, 0, 4);
            push (@seqRes3, @seqResLines);
        }
    }
    
    # Convert 3 letter AA in to 1 letter AA
    my @seqRes1;
    foreach my $aa3 (@seqRes3 )
    {
        chomp $aa3;
        my $aa1 = aa3to1 ($aa3);
        push (@seqRes1, $aa1);
    }

    # Converting array to string of sequence
    my $chainResSeq = join ("", @seqRes1);
    
    return $chainResSeq;    
}

# ************* readFileDataInArray *****************
# Description: Reads a file data into an array line by line
# Inputs: A file path 
# Outputs: An array with file data
# Subroutine call/Testing: my @pdbFile = readFileDataInArray ($pdbPath);
# Date: 26 June 2014
# Author: Saba

sub readFileDataInArray
{
    my ( $filename ) = @_;
    chomp $filename;
    my @filedata = ( );
    my $FILE_DATA;
    unless ( open ( $FILE_DATA, $filename ) )
    {
        croak "Cannot open file \"$filename\"\n\n";
    }
    @filedata = <$FILE_DATA>;
    close $FILE_DATA;
    return @filedata;
}

# ************* aa3to1 *****************                               
# Description: Returns a single letter code for given 3 letter code of A.A
# Inputs: 3 letter code of A.A
# Outputs: 1 letter code of A.A
# Subroutine call/Testing: my $lg = aa3to1($pro[3])
# Date: 22 April 2015       
# Author: Saba
sub aa3to1 
{
    my ($input) = @_;
    my %three2one = (
	'ALA' => 'A',
	'VAL' => 'V',
	'LEU' => 'L',
	'ILE' => 'I',
	'PRO' => 'P',
	'TRP' => 'W',
	'PHE' => 'F',
	'MET' => 'M',
	'GLY' => 'G',
	'SER' => 'S',
	'THR' => 'T',
	'TYR' => 'Y',
	'CYS' => 'C',
	'ASN' => 'N',
	'GLN' => 'Q',
	'LYS' => 'K',
	'ARG' => 'R',
	'HIS' => 'H',
	'ASP' => 'D',
	'GLU' => 'E',
	);

    # clean up the input
    $input =~ s/\n/ /g;
    
    my $seq = '';
    # This use of split separates on any contiguous whitespace
    my @code3 = split(' ', $input);
    
    foreach my $code (@code3) 
    {
        # A little error checking
        if(not defined $three2one{$code}) {
            print "Code $code not defined\n";
            next;
        }
        $seq .= $three2one{$code};
    }
    return $seq;
}




sub hasHapten
{
    my ($pdbPath, $antibodyPair_ARef) = @_;
    my $hashapten = $SFPerlVars::hashapten;
    
    my @antibodyPairs = @{$antibodyPair_ARef};
    my ($hapten, @haptens, @allHaptens);
    
    for my $antibody (@antibodyPairs  )
    {
            my ($L, $H) = split ("", $antibody);
            @haptens = `$hashapten -l $L, -h $H $pdbPath`;
            push (@allHaptens , @haptens );
    }

    if (grep (/^HAPTEN/, @allHaptens) )
    {
        $hapten = 1;
    }
    else {
        $hapten = 0;
    }
    
    return $hapten;
    
}

sub processHapten
{
    my ($pdbPath, $antibodyPair_ARef ) =@_;

    my @antibodyPairs = @ {$antibodyPair_ARef};

    foreach  my $antibodyPair (@antibodyPairs)
    {
        my $antibodyFile = $antibodyPair."_num.pdb";
        my $outputFile = $antibodyPair."_hap.pdb";
        `pdbaddhet $pdbPath $antibodyFile $outputFile`;
        
    }
}


sub movePDBs
{
    my ($dir, $dest, $pdb_id) = @_;
    opendir ( my $DIR, $dir ) or
        die "Can not open directory $dir";
    foreach my $fl ( grep m/$pdb_id/, readdir ( $DIR ) )
    {
        my $from = $dir."/".$fl;
        move $from, $dest;
    }
}

sub dirOperations
{
    my ($process_dir, $pdb_id) = @_;
    mkdir $process_dir;
    chdir $process_dir;
    mkdir $pdb_id;
    chdir $pdb_id;
    my $dir = $process_dir.'/'.$pdb_id;
    return $dir;
}

sub mapChainsIDs
{
    my ($abPair, $chainIdChainTpye_HRef, %complexInfo) = @_;
    my %mapedChains;
    my %chainIdChainTpye = %{$chainIdChainTpye_HRef};
    
    if ( length ($abPair) == 2 )
    {
        my ($L, $H);
        ($L, $H) = split ("", $abPair);
        $mapedChains{'L'} = $L;
        $mapedChains{'H'} = $H;
        if ( !%complexInfo) {
            $mapedChains{'A'} = [];
        }
        else {
            $mapedChains{'A'} = $complexInfo{$abPair};
        }
        
    }
    else {
        my $chainLabel;
        my $key;
        
        foreach $key ( keys %{$chainIdChainTpye{$abPair}} )
        {
            $chainLabel = $key;
        }
        
        if ( $chainLabel eq "L")
        {
            $mapedChains{'L'} = $abPair;
        }
        elsif ( $chainLabel eq "H")
        {
            $mapedChains{'H'} = $abPair;
        }
        
        if ( $complexInfo{$abPair}) {
            $mapedChains{'A'} = $complexInfo{$abPair};
        }
        else {
            $mapedChains{'A'} = [];
        }
        
    }
    
    return %mapedChains;
}
    
sub printHeader
{
    my ($INFILE, $numbering, $pdbPath, %mapedChains ) = @_;
    my %resInfo = getResolInfo($pdbPath);
    print Dumper (\%mapedChains );
    
    my $L = $mapedChains{L};
    my $H = $mapedChains{H};
    my $AgRef = $mapedChains{A};
    
    my @ag = @{$AgRef};

    headerBasic($INFILE, $numbering, %resInfo);    
    headerLHchains ($INFILE, $pdbPath, $L, $H, \@ag);
}

sub headerLHchains
{
    my ($INFILE, $pdbPath, $L, $H, $ag_ARef) = @_;
    my @ag = @{$ag_ARef};
    my ($sym, $AgID);
    
    if ( ($L) and ($H) and (!@ag ) ) 
    {
        print $INFILE "REMARK 950 CHAIN L    L    $L\n";
        print $INFILE "REMARK 950 CHAIN H    H    $H\n";
        print $INFILE "REMARK 950 ", `pdbheader -c $L -m $pdbPath`;
        print $INFILE "REMARK 950 ", `pdbheader -c $L -s $pdbPath`;
        print $INFILE "REMARK 950 ", `pdbheader -c $H -m $pdbPath`;
        print $INFILE "REMARK 950 ", `pdbheader -c $H -s $pdbPath`;
        
    }
    elsif ( ($L) and ($H) and ( @ag ) )
    {
        print $INFILE "REMARK 950 CHAIN L    L    $L\n";
        print $INFILE "REMARK 950 CHAIN H    H    $H\n";
        foreach my $Ag ( @ag ) {
            if ( ( $Ag eq "%l") or ($Ag eq "%h") )
            {
                ($sym, $AgID) = split ("", $Ag);
                $Ag = uc ($AgID);
                print $INFILE "REMARK 950 CHAIN A    $AgID    $Ag\n";
                
            }
            else {
                print $INFILE "REMARK 950 CHAIN A    $Ag    $Ag\n";
            }
            
        }
        print $INFILE "REMARK 950 ", `pdbheader -c $L -m $pdbPath`;
        print $INFILE "REMARK 950 ", `pdbheader -c $L -s $pdbPath`;
        print $INFILE "REMARK 950 ", `pdbheader -c $H -m $pdbPath`;
        print $INFILE "REMARK 950 ", `pdbheader -c $H -s $pdbPath`;

        foreach my $Ag ( @ag ) {
            if ( ( $Ag eq "%l") or ($Ag eq "%h") )
            {
                ($sym, $AgID) = split ("", $Ag);
                $Ag = uc ($AgID);
            }
            print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
        }
    }

    elsif ( ($L) and (!$H) and (!@ag) )
        {
            print $INFILE "REMARK 950 CHAIN L    L    $L\n";
            print $INFILE "REMARK 950 ", `pdbheader -c $L -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $L -s $pdbPath`;
        }
    elsif ( ($L) and (!$H) and (@ag) )
        {
            print "I AM LIGHT_ANTIGEN\n";
            
            print $INFILE "REMARK 950 CHAIN L    L    $L\n";
            foreach my $Ag ( @ag ) {
                print $INFILE "REMARK 950 CHAIN A    $Ag    $Ag\n";
            }

            print $INFILE "REMARK 950 ", `pdbheader -c $L -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $L -s $pdbPath`;
            foreach my $Ag ( @ag) {
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
            }
            
        }

    elsif ( (!$L) and ($H) and (!@ag) )
        {
            print $INFILE "REMARK 950 CHAIN H    H    $H\n";
            print $INFILE "REMARK 950 ", `pdbheader -c $H -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $H -s $pdbPath`;
        }
    elsif ( (!$L) and ($H) and (@ag) )
        {
            print "I AM HEAVY_PROTEIN\n";
            
            print $INFILE "REMARK 950 CHAIN H    H    $H\n";
            foreach my $Ag ( @ag ) {
                print $INFILE "REMARK 950 CHAIN A    $Ag    $Ag\n";
            }

            print $INFILE "REMARK 950 ", `pdbheader -c $H -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $H -s $pdbPath`;
            foreach my $Ag ( @ag) {
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
            }
            
        }
}

sub headerBasic
{
    my ($INFILE, $numbering, %resInfo) = @_;
    print $INFILE "REMARK 950 NUMBERING  $numbering\n";
    print $INFILE "REMARK 950 METHOD     $resInfo{Type}\n";
    print $INFILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
    print $INFILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
    print $INFILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
    print $INFILE "REMARK 950 CHAIN-TYPE LABEL ORIGINAL\n";
    
}
    

sub getResolInfo
    {
        my ($pdbFile) = @_;
        my %resInfo;
        my ($tags, $val);

        my @resolInfo = `pdbheader -r -n $pdbFile`;
        foreach my $line (@resolInfo) {
            chomp $line;
            ($tags, $val) = split (/:\s+/, $line);
            $resInfo{$tags} = $val;

        }
        return %resInfo;

    }

sub getProcessedABchains
{
    my (@fileData) = @_;
    my ($L, $H);
           
    foreach my $line (@fileData)
    {
        if ( $line =~ m/CHAIN L/ )
        {
            my @lineColmn = split (/\s+/, $line);
            $L = $lineColmn[5];
        }
        elsif ( $line =~ m/CHAIN H/ )
        {
            my @lineColmn = split (/\s+/, $line);
            $H = $lineColmn[5];
        }
    }
    return ($L, $H);
}

#  LocalWords:  pdbPath
