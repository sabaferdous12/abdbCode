package antibodyProcessing;

use strict;
#use warnings;
use Data::Dumper;
use List::MoreUtils qw (uniq);
use config;
use Carp;
use Cwd;
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
                        getResolInfo
                        mergeLH
                        dealNumError
                        getChainLabels
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
        my $pdbprep = $config::pdbprep;
        my $pdbext = $config::pdbext;
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
    my ($pdbPath, $pdb, $LOG) = @_;
    my $idabchain = $config::idabchain;
    my $idabFiledir = $config::idabFiledir;    
    my @chainInfo;
    
    # Check if pdb code exists in idabMisLabel.dat file
    # if yes then get info from file otherwise prompt idabachain program
    my $nPdb = `grep $pdb $idabFiledir/idabMisLabel.dat`;
    chomp $nPdb;
    
    if ( $nPdb ) {
        my $idabInfo = $nPdb.".dat";
        open (my $IN, '<', "$idabFiledir/$idabInfo") or die "Can not open file: $idabInfo";
        @chainInfo = <$IN>;
        print {$LOG} "Chain Type information is taken from manual fixed idabchain cases\n";
    }
    # Gets chain IDs and chain Types from program idabchain
    else {
        @chainInfo = `$idabchain $pdbPath`;
        print {$LOG} "Chain Type information is taken from idabchain program\n";
    }
    
    my %chainType = ();
    my %chainIdChainTpye;
    my $chainResSeq;

    $chainType{'Antigen'} = [];
    $chainType{'Light'} = [];
    $chainType{'Heavy'} = [];
    $chainType{'Light/Heavy'} = [];
    $chainType{'Heavy/Light'} = [];
    
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
            elsif ( $chainType eq "Light/Heavy" )
                {    
                    push ( @{ $chainType{'Light/Heavy'}}, $chainid );
                    # Get SEQRES sequence from PDB file (SEQRES records)
                    $chainResSeq = getSEQRESFromPDB($pdbPath, $chainid);
                    $chainIdChainTpye{$chainid}{"LH"} = $chainResSeq; 
                }

            elsif ( $chainType eq "Heavy/Light" )
                {    
                    push ( @{ $chainType{'Heavy/Light'}}, $chainid );
                    # Get SEQRES sequence from PDB file (SEQRES records)
                    $chainResSeq = getSEQRESFromPDB($pdbPath, $chainid);
                    $chainIdChainTpye{$chainid}{"HL"} = $chainResSeq; 
                }

            
        }
    
    return (\%chainType, \%chainIdChainTpye);
}



sub getChains
{
    my ($pdbPath, $pdb, $LOG) = @_;
    my ($chainType_HRef, $chainIdChainTpye_HRef) =
        getChainTypeWithChainIDs ($pdbPath, $pdb, $LOG);
    my %chainType = %{$chainType_HRef};

#    print Dumper (\%chainType);
    
    my (@antigen, @light, @heavy, @LHhybrid, @HLhybrid);
    @antigen = @{$chainType{'Antigen'}};
    @light = @{$chainType{'Light'}};
    @heavy = @{$chainType{'Heavy'}};
    @LHhybrid = @{$chainType{'Light/Heavy'}};
    @HLhybrid = @{$chainType{'Heavy/Light'}};
    
    return (\@light, \@heavy, \@antigen, \@LHhybrid, \@HLhybrid);
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
    my $pdbSplitChains = $config::pdbsplitchains;
    my $pdbatoms = $config::pdbatoms;
    my $pdbgetchain = $config::pdbgetchain;
        
    # To strip header and keeping only ATOMs and then spliting PDB into chains
    `$pdbatoms $pdbPath | $pdbSplitChains`;
    my @PDBchainsFiles = readDirPDB($dir); # Reads all PDB files in directory
    my @PDBchains;
    
    foreach my $chain (@PDBchainsFiles)
        {
            my ($chainID, $ext) = split (/\./, $chain);
            # To discard HETATM by selecting each chain and keeping
            # only ATOM lines                                       
            `$pdbgetchain -a $chainID $chain >$chainID`;
            `mv $chainID $chainID.pdb`; # Rename A as A.pdb
            push (@PDBchains, $chainID);           
        }
    close $PDBFILE;
    return @PDBchains;
    
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
    my $renumpdb = $config::renumpdb;
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
sub getChainLabels
    {
        my ($file) = @_;
        my @cols;
        my %chains;
        open (my $IN, '<', $file) or die "Can not find file $file\n";         
        while (my $line = <$IN>)
            {

                chomp $line;
                if ($line =~ /^ATOM/) {
                
                @cols = split(' ', $line);
                $chains{$cols[4]} = 1;
            }
                
            }
        my @ch = keys %chains;
        
        return @ch;
        
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
    my ( $antibodyPairs, $ab, $LOG) = @_;
    my ( $out, $error, $success, $CDR );
    my $numberedAntibody;
    my @numFailedPair;    
    foreach my $pair ( @$antibodyPairs )
    {
       $numberedAntibody = $pair."_num.pdb";         
      my @chains = getChainLabels($numberedAntibody);
             
     #  my $numLines =`cat $numberedAntibody | wc -l`;
       
        my $getpdb= $config::getpdb;
        my $CDRsPDB = $pair."_CDR.pdb";
        my $FWsPDB = $pair."_FW.pdb";
         
        my @numbering;
        
        if ( $ab eq "LH")
        {
            if ( (scalar @chains ) < 2) # Checks if there is any chain missing
                {
                print {$LOG} "One of the chain is missing: FWs'/CDRs' can not extracted\n";
                push (@numFailedPair, $pair);
                next; 
            }

            else {
                
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
        }
       
        
        elsif ( $ab eq "L")
        {
            if ( (scalar @chains ) < 1)
                {
                print {$LOG} "One of the chain is missing: FW'/CDRs' can not extracted\n";
                push (@numFailedPair, $pair);
                next;
            }
            else
                {
                    @numbering =  (["L24", "L34"], ["L50", "L56"], ["L89", "L97"]);
                    `$getpdb -v L24 L34 $numberedAntibody |
$getpdb -v L50 L56 |
$getpdb -v L89 L97 >>$FWsPDB`;
                }
        }
       
        
        elsif ( $ab eq "H")
        {
            if ((scalar @chains ) < 1) {
                print {$LOG} "One of the chain is missing: FWs'/CDRs' can not extracted\n";
                push (@numFailedPair, $pair);
                next;
            }
            else {
                @numbering =  (["H31", "H35"], ["H50", "H65"], ["H95", "H102"] );
                `$getpdb -v H31 H35 $numberedAntibody | 
$getpdb -v H50 H65 |
$getpdb -v H95 H102 >>$FWsPDB`;
            }
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
    return @numFailedPair;
    
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
    my $hashapten = $config::hashapten;
    
    
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
    my ($pdbPath, $antibodyPair_ARef,$ab, $LOG ) =@_;
    my $hashapten = $config::hashapten;
    my $haptenFlag=0;
    my %fileType;
    
    my @antibodyPairs = @ {$antibodyPair_ARef};
    my (@hapMole, @hapInter);
    
    foreach  my $antibodyPair (@antibodyPairs)
    {
        
        #mergeLH($LOG, @antibodyPairs);

        @hapMole = ();
        @hapInter = ();
        my @haptens;
        
        if ( $ab eq "L") {
            @haptens = `$hashapten -l $antibodyPair $pdbPath`;
        }
        elsif (  $ab eq "H" ) {
            @haptens = `$hashapten -h $antibodyPair $pdbPath`;
        }
        elsif ( $ab eq "LH" ){
            my ($L, $H) = split ("", $antibodyPair);
            @haptens = `$hashapten -l $L -h $H $pdbPath`;
        }
        
        if ( ! @haptens) {
            $fileType{$antibodyPair} = "num";
            next;
        }
        
        # This For Loop Is To Store hapten molecule name from
        # hasHapten program output.
        foreach my $hapten( @haptens ) {
            my @hap = split (" ", $hapten);
            my $mol = $hap[1]; # molecule name
            my $molInteraction = $hap[2]; # chain label and resSeq

            push (@hapInter, $molInteraction);
            push (@hapMole, $mol); 
        }
        # To remove duplicates
        @hapMole = uniq (@hapMole);
        @hapInter = uniq (@hapInter);
        
        #print "Test: @hapMole:   @hapInter\n";
        
        my $antibodyFile = $antibodyPair."_num.pdb";
        my $outputFile_temp = $antibodyPair."_hap_temp.pdb";
        my $outputFile = $antibodyPair."_hap.pdb";
        open (my $OUT, '>', $outputFile);
                
        # Writing all HETATM records to antibody file
        if (-z !$antibodyFile) {
            `pdbaddhet $pdbPath $antibodyFile $outputFile_temp`;
           }     
        # Parsing HETATMs for only identified haptens and discarding the rest
        open (my $IN, '<', $outputFile_temp);
        my @HET = <$IN>;
                
        my @atomRec = grep { $_ =~ m/^ATOM|^TER/} @HET;
        print {$OUT} @atomRec;
                
        # Obtaining HETATMS for Haptens
        foreach my $hapMol( @hapMole ) {
            # Obtaining Haptens for multiple resSeq and Chains
            foreach my $molInter( @hapInter) {
                my $chainId = substr ($molInter, 0, 1);
                my $resSeq = substr ($molInter, 1);

                # To grep this pattern SO4 L 214 from HETATM records
                my @hapHetRec = grep { $_ =~ m/$hapMol $chainId $resSeq/ |
                                           $_ =~ m/$hapMol $chainId  $resSeq/ |
                                               $_ =~ m/$hapMol $chainId   $resSeq/ |
                                                   $_ =~ m/$hapMol $chainId$resSeq/ } @HET;
                # This check is to eliminate any small solvent molecule like: SO4, PO4
                                
                if ( scalar @hapHetRec > 5) {
                    print {$OUT} @hapHetRec;
                    $haptenFlag = 1;
                    $fileType{$antibodyPair} = "hap";
                                        
                }
                else {
                    $fileType{$antibodyPair} = "num";
                    next;
                }
                                
            }
                        
        }
                
    }#abtibody pair
        
    return %fileType;
    
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
    my ($ab, $abPair, $chainIdChainTpye_HRef, %complexInfo) = @_;
    my %mapedChains;
    my %chainIdChainTpye = %{$chainIdChainTpye_HRef};
    
    if ( ( length ($abPair) == 2 ) and ($ab eq "LH") )
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

    if ( ( length ($abPair) == 2 ) and ( ($ab eq "L") or ($ab eq "H") ) )
    {
        my ($c1, $c2);
        ($c1, $c2) = split ("", $abPair);
        if ( $ab eq "L") {
            $mapedChains{'L1'} = $c1;
            $mapedChains{'L2'} = $c2;
        }
        elsif ( $ab eq "H" ) {
            $mapedChains{'H1'} = $c1;
            $mapedChains{'H2'} = $c2;
        }
        
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
    my ($ab, $INFILE, $numbering, $pdbPath, %mapedChains ) = @_;
    my %resInfo = getResolInfo($pdbPath);
    my ($L, $H);

#    print "GGG\n";
#    print Dumper (\%mapedChains);
    
    my $AgRef = $mapedChains{A};
    my @ag = @{$AgRef};

    if ( $ab eq "L") {
        if ( exists $mapedChains{L}) {
            $L = $mapedChains{L};
        }
        else {
            $L = $mapedChains{L1};
            $H = $mapedChains{L2};
        }
        
     }
    elsif ( $ab eq "H")
        {
            if ( exists $mapedChains{H}) {
            $H = $mapedChains{H};
        }
        else {
            $L = $mapedChains{H1};
            $H = $mapedChains{H2};
        }
        }
    
    elsif ( $ab eq "LH") {
         $L = $mapedChains{L};
         $H = $mapedChains{H};
         }
    headerBasic($INFILE, $numbering, %resInfo);
    headerLHchains ($ab, $INFILE, $pdbPath, $L, $H, \@ag, \%mapedChains);
    
}
# The subroutine extracts SEQRES from PDB file for given chain label...
sub getSEQRES
{
    my ($pdbPath, $chainID, $chainLabel) = @_;
    open (my $IN, $pdbPath) or die "Error in opening file\n";
    my @fileContent = <$IN>;
    my @seqRes;
    
    foreach my $line (@fileContent ) {
        if ( ( $line =~ m/SEQRES/ ) and ($line =~ m/\s+$chainID\s+/) )
        {
            push (@seqRes, $line);
        }
    }
    # Replace the Light and heavy chain labels with L or H
    # This is to have consistency in the headers of  PDB output file 
    if ( $chainLabel ) {
        @seqRes = map {s/\s+$chainID\s+/ $chainLabel /g; $_; } @seqRes; 
    }
    
   return @seqRes;
}     
          

sub headerLHchains
{
    my ($ab, $INFILE, $pdbPath, $L, $H, $ag_ARef, $mappedChain_HRef) = @_;
    my @ag = @{$ag_ARef};
    my %mappedChains = %{$mappedChain_HRef};
#   sort {$mappedChains{$a} <=> $mappedChains{$b} or $a cmp $b} values %mappedChains ;
    
    my ($sym, $AgID);
    my ($LL, $HH);

        $LL = $L;
        $HH = $H;
    
    # if antibody chain as real labels as 0
    if (  $L eq 0 ){
        $L = 1;
    }
    if ( $H eq 0 ) {
        $H = 1;
    }

# print Headers for complete Antibody
    if ( ($L) and ($H) and (!@ag ) ) 
    {
        # Check if its L and H but not single chain dimer
        if ( (exists $mappedChains{L}) and (exists $mappedChains{H}) ) {
            print $INFILE "REMARK 950 CHAIN L    L    $LL\n";
            print $INFILE "REMARK 950 CHAIN H    H    $HH\n";
            print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $LL -s $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $HH -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $HH -s $pdbPath`;
            my @seqRes = getSEQRES($pdbPath, $LL, "L");
            print $INFILE @seqRes;

            @seqRes = getSEQRES($pdbPath, $HH, "H");
            print $INFILE @seqRes;
        }
        else {
            # if its a light dimer
            if ( $ab eq "L") {
                print $INFILE "REMARK 950 CHAIN L    L    $LL\n";
                print $INFILE "REMARK 950 CHAIN L    L    $HH\n";
                print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $LL -s $pdbPath`;
                # Same $LL for both chain because one chain is symmetry dimer and
                # is not present in original PDB
                #print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
                #`sed -i '0,/MOLECULE : C :/! s/MOLECULE : C :/MOLECULE : $HH :/' $dir/$newFile`;
                                
                my @seqRes = getSEQRES($pdbPath, $LL, "L");
                print $INFILE @seqRes;

                @seqRes = getSEQRES($pdbPath, $LL, "L");
                print $INFILE @seqRes;
                                 
            }

            elsif ( $ab eq "H") {
                print $INFILE "REMARK 950 CHAIN H    H    $LL\n";
                print $INFILE "REMARK 950 CHAIN H    H    $HH\n";
                print $INFILE "REMARK 950 ", `pdbheader -c $HH -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $HH -s $pdbPath`;
          
                my @seqRes = getSEQRES($pdbPath, $HH, "H");
                print $INFILE @seqRes;

                @seqRes = getSEQRES($pdbPath, $HH, "H");
                print $INFILE @seqRes;
            }
        }
        
    }
    
  
    elsif ( ($L) and ($H) and ( @ag ) )
    {
        if ( (exists $mappedChains{L}) and (exists $mappedChains{H}) ) 
            {
                print $INFILE "REMARK 950 CHAIN L    L    $LL\n";
                print $INFILE "REMARK 950 CHAIN H    H    $HH\n";
# To print Antigen Header information 
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
                                
                print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $LL -s $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $HH -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $HH -s $pdbPath`;
                
                foreach my $Ag ( @ag ) {
                    if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                        {
                            ($sym, $AgID) = split ("", $Ag);
                            $Ag = uc ($AgID);
                        }
                    print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
                    print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
                }

                # Printing SeqRes
                # Printing SEQRES for antibody
                my @seqRes = getSEQRES($pdbPath, $LL, "L");
                print $INFILE @seqRes;
                
                @seqRes = getSEQRES($pdbPath, $HH, "H");
                print $INFILE @seqRes;
                # Printing SEQRES for antigen
                foreach my $Ag ( @ag ) {
                    if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                        {
                            ($sym, $AgID) = split ("", $Ag);
                            $Ag = uc ($AgID);
                        }
                    @seqRes = getSEQRES($pdbPath, $Ag);
                    print $INFILE @seqRes;
                    
                }
                
            }
        else {
            if ($ab eq "L")
                {
                    print $INFILE "REMARK 950 CHAIN L    L    $LL\n";
                    print $INFILE "REMARK 950 CHAIN L    L    $HH\n";

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
                # For Dimer print header information just for one chain                    
                    print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
                    print $INFILE "REMARK 950 ", `pdbheader -c $LL -s $pdbPath`;

                    foreach my $Ag ( @ag ) {
                        if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                        {
                            ($sym, $AgID) = split ("", $Ag);
                            $Ag = uc ($AgID);
                        }
                        
                        print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
                        print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
                    }
                    
                        # Printing SeqRes
                        # Printing SEQRES for antibody
                        my @seqRes = getSEQRES($pdbPath, $LL, "L");
                        print $INFILE @seqRes;
                        
                        @seqRes = getSEQRES($pdbPath, $HH, "L");
                        print $INFILE @seqRes;
                        # Printing SEQRES for antigen
                        foreach my $Ag ( @ag ) {
                            if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                                {
                                    ($sym, $AgID) = split ("", $Ag);
                                    $Ag = uc ($AgID);
                                }
                            @seqRes = getSEQRES($pdbPath, $Ag);
                            print $INFILE @seqRes;
            
                        }
                                        
                }
            elsif ($ab eq "H")
                {
                    print $INFILE "REMARK 950 CHAIN H    H    $LL\n";
                    print $INFILE "REMARK 950 CHAIN H    H    $HH\n";

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
                                        
                    print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
                    print $INFILE "REMARK 950 ", `pdbheader -c $LL -s $pdbPath`;

                    foreach my $Ag ( @ag ) {
                        if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                        {
                            ($sym, $AgID) = split ("", $Ag);
                            $Ag = uc ($AgID);
                        }
                        
                        print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
                        print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
                    }

                    # Printing SEQRES for antibody
                    my @seqRes = getSEQRES($pdbPath, $LL, "H");
                    print $INFILE @seqRes;
                    
                    @seqRes = getSEQRES($pdbPath, $HH, "H");
                    print $INFILE @seqRes;
                    # Printing SEQRES for antigen
                    foreach my $Ag ( @ag ) {
                        if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                            {
                                ($sym, $AgID) = split ("", $Ag);
                                $Ag = uc ($AgID);
                            }
                        @seqRes = getSEQRES($pdbPath, $Ag);
                        print $INFILE @seqRes;
            
                    }

                }
        } 
        
            
    }
    
# Print headers for light chain only 
    elsif ( ($L) and (!$H) and (!@ag) )
        {
            print $INFILE "REMARK 950 CHAIN L    L    $LL\n";
            print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $LL -s $pdbPath`;

            my @seqRes = getSEQRES($pdbPath, $LL, "L");
            print $INFILE @seqRes;
           
        }
    elsif ( ($L) and (!$H) and (@ag) )
        {
            
            print $INFILE "REMARK 950 CHAIN L    L    $LL\n";
            foreach my $Ag ( @ag ) {
            if ( ( $Ag eq "%l") or ($Ag eq "%h") )
            {
                ($sym, $AgID) = split ("", $Ag);
                $Ag = uc ($AgID);
            }
                print $INFILE "REMARK 950 CHAIN A    $Ag    $Ag\n";
            }

            print $INFILE "REMARK 950 ", `pdbheader -c $LL -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $LL -s $pdbPath`;
            
            foreach my $Ag ( @ag) {
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
            }

            my  @seqRes = getSEQRES($pdbPath, $LL, "L");
            print $INFILE @seqRes;
            
            foreach  my $Ag ( @ag ) {
                if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                 {
                     ($sym, $AgID) = split ("", $Ag);
                     $Ag = uc ($AgID);
                 }
                @seqRes = getSEQRES($pdbPath, $Ag);
                print $INFILE @seqRes;
                
            }
                        
        }

# Print header for Heavy chain only     
    elsif ( (!$L) and ($H) and (!@ag) )
        {
            print $INFILE "REMARK 950 CHAIN H    H    $HH\n";
            print $INFILE "REMARK 950 ", `pdbheader -c $HH -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $HH -s $pdbPath`;
            my @seqRes = getSEQRES($pdbPath, $HH, "H");
            print $INFILE @seqRes;
            
        }
    elsif ( (!$L) and ($H) and (@ag) )
        {
            print $INFILE "REMARK 950 CHAIN H    H    $HH\n";
            foreach my $Ag ( @ag ) {
                if ( ( $Ag eq "%l") or ($Ag eq "%h") )
                {
                    ($sym, $AgID) = split ("", $Ag);
                    $Ag = uc ($AgID);
                }   
                print $INFILE "REMARK 950 CHAIN A    $Ag    $Ag\n";
            }

            print $INFILE "REMARK 950 ", `pdbheader -c $HH -m $pdbPath`;
            print $INFILE "REMARK 950 ", `pdbheader -c $HH -s $pdbPath`;
            foreach my $Ag ( @ag) {
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -m $pdbPath`;
                print $INFILE "REMARK 950 ", `pdbheader -c $Ag -s $pdbPath`;
            }
            my @seqRes = getSEQRES($pdbPath, $HH, "H");
            print $INFILE @seqRes;
            
            foreach  my $Ag ( @ag ) {
                if  ( ( $Ag eq "%l") or ($Ag eq "%h") )
                {
                    ($sym, $AgID) = split ("", $Ag);
                    $Ag = uc ($AgID);
                }
                @seqRes = getSEQRES($pdbPath, $Ag);
                print $INFILE @seqRes;
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

sub mergeLH
    {
        my ( $LOG, @antibodyPairs) = @_;

        my $numberedAntibody;
        my $countFailedPairs=0;
        foreach my $antibodyPair (@antibodyPairs)
        {
            chomp $antibodyPair;
                     
        my ($l, $h) = split ("", $antibodyPair);

            
       #     if ( ($l) and ($h)) {
        #        print "ASSSEM: $l and $h\n";
                $numberedAntibody = $antibodyPair."_num.pdb";
                
                    `cat $l"_num.pdb" $h"_num.pdb" >$numberedAntibody`;

                    print {$LOG} "Each of the antibody in the PDB has been assembled with".
                        " light and heavy chain in one file\n"; 
                    my @chains = getChainLabels ($numberedAntibody);
                    
                    #my $numLines =`cat $numberedAntibody | wc -l`;
                    # To check if there is one chain in the $antibodyPair.pdb
                    if ( (scalar @chains) < 2 )
                        {
                            # Another check to confirm if kabat numbering has failed
                        my $errorCheck=`grep "pdbpatchnumbering" numberingFailedChain.dat`;
                        #my $chainCheck=`grep $lchain|$hchain numberingFailedChain.dat`;
                        
                        if ( $errorCheck ){
                            print {$LOG} "Kabatnum Error-Numbering failed on this pair:$antibodyPair\n";
                            $countFailedPairs++;
                        }
                    }
         #       }
    }
        
        return $countFailedPairs++;
        
    }

sub dealNumError
{
    my ($LOG, @singleChainAb)=@_;
    my $countFailedPairs=0;
    foreach my $chain (@singleChainAb )
        {
            chomp $chain;
            
            my $nchain=$chain."_num.pdb";
            my $errorCheck=`grep "pdbpatchnumbering" numberingFailedChain.dat`;      
            my $chainCheck = `grep $nchain numberingFailedChain.dat`;
            
            if ( ($errorCheck) and ($chainCheck))
            {
                print {$LOG} "Kabatnum Error-Numbering failed on this chain:$chain\n";
                $countFailedPairs++;
            }
        }
    return $countFailedPairs;

}
