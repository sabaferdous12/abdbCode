#!/acrm/usr/local/bin/perl -s

#*************************************************************************
#
#   Program:    processAntibodyPDBs
#   File:       processAntibodyPDBs.pl
#   
#   Version:    V1.1
#   Date:       22.04.14
#   Function:   This script reads a list of antibody PDBs (obtained from SACs)
#               and process it. The procssing includes 1) antibody numbering
#               according to given numbering scheme 2) splitting PDB file into
#               multiple PDB files for antibodies present in individual PDB
#               3) Moves the processed antibody/antigen complexes, antibody 
#               only and antibody/hapten complxes into separate directories
#
#   Modules:    pdb.pm, pdbWrap
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
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#   This script reads a list of antibody PDBs (obtained from SACs)  
#   and process it. The procssing includes 1) antibody numbering               
#   according to given numbering scheme 2) splitting PDB file into            
#   multiple PDB files for antibodies present in individual PDB              
#   3) Moves the processed antibody/antigen complexes, antibody              
#   only and antibody/hapten complxes into separate directories 
#   
#   
#
#*************************************************************************
#
#   Usage:
#   ======
#  ./processAntibodyPDBs.pl [-k -c -a] [input_file]
#   (Script should be run from the directory where you want to generate data)
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.1   22.04.15 Original
#
#*************************************************************************

use strict;
use IO::CaptureOutput qw(capture qxx qxy);
use Carp;
#use warnings;
use Data::Dumper;
use File::Copy;
use Cwd;
#use lib "./scripts";
use pdb qw(
    get_pdb_path                                                              
    get_file_data                                                             
    get_pdbcode_from_xml                                                      
    check_chain                                                               
    split_pdb_to_chains                                                       
    get_pdbchains_contacts                                                    
    pair_heavy_light_antigen                                                  
    get_hash_key
    antibody_assembly                                                         
    antibody_number                                                           
    extract_CDR                                                               
    assemble_CDR_antigen                                                      
    antigen_CDR_conts                                                         
    get_complex                                                               
    get_antibody_antigen_complex                                              
    make_antibody_complex
    hasHapten
    processHapten
    getResolInfo
    getsingleChainAntibody
      );

use pdbWrap qw(
dirOperations 
processAntibody 
processAntibodyAntigen 
movePDBs
processSingleChainAntibody
processSingleChainAntibodyAntigen);

my $Usage = <<'EOF';
Usage:

    program_name [-k -c -a] [input_file]
 
Options:

numbering scheme

    k : number the antibodies by Kabat numbering scheme
    c : number the antibodies by Chothia numbering scheme
    a : number the antibodies by Martin numbering scheme
    
EOF

if ( defined ( $::help ) ) 
{
    print "$Usage\n";
    exit 0;
}
my ($FreeAntibodies, $AntibodyAntigen, $AntibodyHapten);
my ($LightChain, $LightAntigen,  $LightHapten);
my ($HeavyChain, $HeavyAntigen, $HeavyHapten);
my $numbering;

# Initial numbering scheme variable depending upon user choice on command line
my $nsch;

if ( $::k )
{
    $nsch = '-k';
    $FreeAntibodies  = "FreeAntibody_Kabat";
    $AntibodyAntigen= "AntibodyAntigen_Kabat";
    $AntibodyHapten = "AntibodyHapten_Kabat";

    $LightChain = "LightChain_Kabat";
    $LightAntigen = "LightAntigen_Kabat";
    $LightHapten = "LightHapten_Kabat";
    
    $HeavyChain = "HeavyChain_Kabat";
    $HeavyAntigen = "HeavyAntigen_Kabat";
    $HeavyHapten = "HeavyHapten_Kabat";

    $numbering = "KABAT";
}
elsif ( $::c )
{
    $nsch = '-c';
    $FreeAntibodies  = "FreeAntibody_Chothia";
    $AntibodyAntigen = "AntibodyAntigen_Chothia";
    $AntibodyHapten = "AntibodyHapten_Chothia";

    $LightChain = "LightChain_Chothia";
    $LightAntigen = "LightAntigen_Chothia";
    $LightHapten = "LightHapten_Chothia";
        
    $HeavyChain = "HeavyChain_Chothia";
    $HeavyAntigen = "HeavyAntigen_Chothia";
    $HeavyHapten = "HeavyHapten_Chothia";
    $numbering = "CHOTHIA";
}

elsif ( $::a )
{
    $nsch = '-a';
    $FreeAntibodies  = "FreeAntibody_Martin";
    $AntibodyAntigen = "AntibodyAntigen_Martin";
    $AntibodyHapten = "AntibodyHapten_Martin";
    
    $LightChain = "LightChain_Martin";
    $LightAntigen = "LightAntigen_Martin";
    $LightHapten = "LightHapten_Martin";
    
    $HeavyChain = "HeavyChain_Martin";
    $HeavyAntigen = "HeavyAntigen_Martin";
    $HeavyHapten = "HeavyHapten_Martin";

    $numbering = "MARTIN";
}

# Define and create process directory 
my $master_dir = getcwd;
my $process_dir = $master_dir.'/'.'Dataprep' . $$ . $nsch; # $$ = process id
mkdir $FreeAntibodies;
mkdir $AntibodyAntigen;
mkdir $AntibodyHapten;

mkdir $LightChain;
mkdir $LightAntigen;
mkdir $LightHapten;

mkdir $HeavyChain;
mkdir $HeavyAntigen;
mkdir $HeavyHapten;

    
my $input_file = "@ARGV";

# Here are 2 functions,they can be used depending upon type of input file 
# i.e, .xml or .txt
# Use function get_pdbcode_from_xml, if input file is .xml file
# Use get_file_data, if input file is .txt file

#my @all_pdbs = get_pdbcode_from_xml ( $input_file ) or 
#    die "Can not open $input_file\n";
my @all_pdbs = get_file_data ( $input_file ) or 
    die "Can not open $input_file\n"; 

# ****** Global Variable Declaration ******
my $dir;
my ($count_superseded, $count_kabatnum_error, $count_CDR_error,
    $count_idabError) = (0, 0, 0, 0);
my ($count_FreeAntibody, $count_AbAg, $count_Abhapten, $count_AbDNA_RNA,
    $countNonAg, $countBiAg)
    = (0, 0, 0, 0, 0, 0 );
my ($count_Lg, $count_LgAnt, $count_Lhapten, $count_LDNA_RNA)
    = (0, 0, 0, 0);
my ($count_Hv, $count_HvAnt, $count_Hhapten, $count_HDNA_RNA, $count_Ant)
    = (0, 0, 0, 0, 0); 
my (@Lg, @LgAnt, @LgHapten, @LgDRNA,
    @Hv, @HvAnt, @HvHapten, @HvDRNA,
    @FreeAntibody, @AbAg_complex, @AbAghapten, @AbAgDRNA, @nonAg, @biAg,
    @superseded, @kabat_failed, @CDR_failed, %contact_hash, @Ant, @idab_failed);
my ( $MASTER_LOG, $LOG, $SUMMARY );
my $destABAG =  "$master_dir"."/".$AntibodyAntigen;
my $destAB = "$master_dir"."/".$FreeAntibodies;
my $destLg = "$master_dir"."/".$LightChain;
my $destHv = "$master_dir"."/".$HeavyChain;
my $destLgAG = "$master_dir"."/".$LightAntigen;
my $destHvAG = "$master_dir"."/".$HeavyAntigen;

#my (@LgAnt, @HvAnt, @Lg, @Hv, @Ant);
#my ($LgAnt, $HvAnt, $Lg, $Hv, $Ant) = (0, 0, 0, 0, 0);
my($light_c, $heavy_c, $antigen_c);

# Open Master file for recording information about Superceded, failed
# Kabat numbering and details of processed PDBs
open ( $MASTER_LOG, ">masterlog".$nsch.".log" );
open ( my $HEADER, ">header.dat") or die "Can not open file $!"; 
open ( my $AGCHAIN, ">AntigenChains.dat") or die "Can not open file $!";
print {$AGCHAIN} "PDB_ID:Antgen Chains\n";
open ( my $ABCHAIN, ">AntibodyChains.dat") or die "Can not open file $!";
print {$ABCHAIN} "PDB_ID:Antibody Chains\n";

# Reading the list of antibodies
foreach my $pdb_id ( @all_pdbs )
{
    print "\nWorking on $pdb_id\n"; 
    chomp $pdb_id;

    $dir = dirOperations ($process_dir, $pdb_id);
    # Open log files                                                          
    open ( $LOG, ">$dir/antibody.log" ) or
        die "Error in creating log file";
    open ( $SUMMARY, ">$dir/cont_summary.txt" ) or
        die "Error in creating log file";
    
    my $file_path = get_pdb_path ( $pdb_id ); # Reading the local pdb on acrm

    # Obtain a hash with resolution information of structure
    my %resoluInfo = getResolInfo($file_path);
        
    # Writing the header information into a flat file
    my @headerInfo =  `pdbheader -p -s -m $file_path`;
    print {$HEADER} @headerInfo;

    # Check if file exists, if it doesn't exist (i.e Superceded) then writes the
    # protein ID in Master file.  
    if ( !( -e "$file_path" ) ) {
	push ( @superseded, $pdb_id );
	print {$MASTER_LOG} "$pdb_id is found Superseded\n\n";
	$count_superseded++;
	next;
    }

    # To find Antigen Chain IDs and chain types
    my $chaintype = $SFPerlVars::chaintype;
    # To obtain the antigen chain label and chain type (N Protein)
    my (@antigenChainType) =
        `$chaintype $file_path`;
        
    # Check antibody for types of chains either light, heavy or antigen
    # This script only processes 3 types of complexes, 
    # 1) Antibody/antigen complexes (light, heavy and antigen)
    # 2) Antibody complexes (light and heavy)
    # 3) Antibodies bound with Hapten  
    # It ignores all other types of complexes e.g. Onlyantigen, Only Light,
    # Only Heavy and Light or Heavy with antigen. However the details about these 
    # are present in master log file
    
    ( $light_c, $heavy_c, $antigen_c ) = check_chain ( $file_path );

    if ( ( !$light_c) and (!$heavy_c) and (!$antigen_c) )
        {
            push ( @idab_failed, $pdb_id );
            $count_idabError++;
            chdir '..';
            next;
        }

     
    my $flag = 0;               # Flag if antigen is a Hapten molecule

    # This subroutine returns true if Hapten molecule
    my $hapten = hasHapten ($file_path);

    if ($hapten) {
        #        print "TST: $pdb_id: $hapten\n";
 	if ( ( $antigen_c eq "" ) and  (( $heavy_c eq 'Heavy' ) or ($heavy_c eq ""))
                 and ( ( $light_c eq 'Light' ) or ($light_c eq "") ) ) {
	    $antigen_c = 'Antigen';
	    $flag = 1; 
	}

    }

    # Checks to count unprocessed Antibody
    if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq "" ) and 
             ( $light_c eq 'Light') ) {
        my ($chainsHashRef, $ABChainsRef, $AGChainsRef) =
            processSingleChainAntibody ($file_path, $LOG, $SUMMARY, $pdb_id);
        my $ab = "Light";
            
        my @antibodyChain = map {split ("", $_)} @{$ABChainsRef};
        print {$ABCHAIN} $pdb_id,":",join (",", @antibodyChain), "\n";
        print {$AGCHAIN} $pdb_id,":",join (",", (map {split ("", $_)} @{$AGChainsRef}) ), "\n";        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eval { antibody_number ( $ABChainsRef, $nsch ); 1; };
	if ( $@ ) {
            foreach my $mLOG ($LOG, $MASTER_LOG) {
                print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
                    $@ . "Program exited\n\n";
            }
                
            push ( @kabat_failed, $pdb_id );
            $count_kabatnum_error++;
            chdir '..';
            next;
        } else {
            print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
                " numbering scheme\n\n";
        }
                
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        eval { extract_CDR ( $ABChainsRef, $ab ); 1; };
        if ( $@ ) {
            foreach my $mLOG ( $LOG, $MASTER_LOG ) {
                print {$mLOG} "Can'nt extract CDRs for: $pdb_id\n" . $@ .
                    "Program exited\n\n";
            }
            push ( @CDR_failed, $pdb_id );
            $count_CDR_error++;
            chdir '..';
            next;
        } else {
            print {$LOG} "CDRs for $pdb_id has been extracted from each light/Heavy".
                "chain\n\n";
        }
        #~~~~~~~~~~~~~~~~~~~~~~~
                
        if ( ($hapten) and ($flag)) {
            my $dest = "$master_dir"."/".$LightHapten;
            processHapten($file_path, $ABChainsRef);
            getsingleChainAntibody ( $ABChainsRef, \$pdb_id, $hapten, 1, $numbering,
                                     $file_path, $chainsHashRef);
            movePDBs ($dir, $dest, $pdb_id);
            push (@LgHapten, $pdb_id); 
            $count_Lhapten++;
            chdir '..';
            next; 
	}

        processSingleChainAntibodyAntigen
            ($ABChainsRef, $AGChainsRef, $pdb_id, $destLgAG, $destLg, $LOG, $SUMMARY,
             $numbering, $file_path, $chainsHashRef);

        
        # To deal with DNA/RNA antigens and move them in Hapten folder (Non-protein)
        if ( grep {/DNA|RNA/} @antigenChainType) {
            my $dest = "$master_dir"."/".$LightHapten;
            movePDBs ($dir, $dest, $pdb_id);
            push (@LgDRNA, $pdb_id);
            $count_LDNA_RNA++;
        } else {
            # Moving complexes from individual directory to a separate folder to 
            # store all antibody/antigen complexes 
            my $dest = "$master_dir"."/".$LightAntigen;
            movePDBs ($dir, $dest, $pdb_id);
            push ( @LgAnt, $pdb_id );
            $count_LgAnt++;
        }
    }
    ##############################################################################
    if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq 'Heavy' ) and 
             ( $light_c eq "" ) ) {
        my ($chainsHashRef, $ABChainsRef, $AGChainsRef) =
            processSingleChainAntibody ($file_path, $LOG, $SUMMARY, $pdb_id);
        my $ab = "Heavy";

        my @antibodyChain = map {split ("", $_)} @{$ABChainsRef};
        print {$ABCHAIN} $pdb_id,":",join (",", @antibodyChain), "\n";
        print {$AGCHAIN} $pdb_id,":",join (",", (map {split ("", $_)} @{$AGChainsRef}) ), "\n";

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eval { antibody_number ( $ABChainsRef, $nsch ); 1; };
	if ( $@ ) {
            foreach my $mLOG ($LOG, $MASTER_LOG) {
                print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
                    $@ . "Program exited\n\n";
            }
                
            push ( @kabat_failed, $pdb_id );
            $count_kabatnum_error++;
            chdir '..';
            next;
        } else {
            print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
                " numbering scheme\n\n";
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        eval { extract_CDR ( $ABChainsRef, $ab ); 1; };
        if ( $@ ) {
            foreach my $mLOG ( $LOG, $MASTER_LOG ) {
                print {$mLOG} "Can'nt extract CDRs for: $pdb_id\n" . $@ .
                    "Program exited\n\n";
           }
            push ( @CDR_failed, $pdb_id );
            $count_CDR_error++;
            chdir '..';
            next;
        } else {
            print {$LOG} "CDRs for $pdb_id has been extracted from each light/Heavy".
                "chain\n\n";
        }
      

        if ( ($hapten) and ($flag)) {
            my $dest = "$master_dir"."/".$HeavyHapten;
            processHapten($file_path, $ABChainsRef);
            getsingleChainAntibody ( $ABChainsRef, \$pdb_id, $hapten, 1, $numbering,
                                     $file_path, $chainsHashRef);
            movePDBs ($dir, $dest, $pdb_id);
            push (@HvHapten, $pdb_id); 
            $count_Hhapten++;
            chdir '..';
            next; 
	}

                
        processSingleChainAntibodyAntigen
            ($ABChainsRef, $AGChainsRef, $pdb_id, $destHvAG, $destHv,  $LOG,
             $SUMMARY, $numbering, $file_path, $chainsHashRef);

        
        # To deal with DNA/RNA antigens and move them in Hapten folder (Non-protein)
        if ( grep {/DNA|RNA/} @antigenChainType) {
            my $dest = "$master_dir"."/".$HeavyHapten;
            movePDBs ($dir, $dest, $pdb_id);
            push (@HvDRNA, $pdb_id);
            $count_HDNA_RNA++;
        } else {
            # Moving complexes from individual directory to a separate folder to 
            # store all antibody/antigen complexes 
            my $dest = "$master_dir"."/".$HeavyAntigen;
            movePDBs ($dir, $dest, $pdb_id);
            push ( @HvAnt, $pdb_id );
            $count_HvAnt++;
        }
    }
    ############################################################################
    if ( ( $antigen_c eq "" ) and ( $heavy_c eq 'Heavy' ) and 
             ( $light_c eq "" ) ) {
        my ($chainsHashRef, $ABChainsRef) =
            processSingleChainAntibody ($file_path, $LOG, $SUMMARY, $pdb_id);
        my @antibodyChain = map {split ("", $_)} @{$ABChainsRef};
        print {$ABCHAIN} $pdb_id,":",join (",", @antibodyChain), "\n";
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eval { antibody_number ( $ABChainsRef, $nsch ); 1; };
	if ( $@ ) {
            foreach my $mLOG ($LOG, $MASTER_LOG) {
                print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
                    $@ . "Program exited\n\n";
            }
                
            push ( @kabat_failed, $pdb_id );
            $count_kabatnum_error++;
            chdir '..';
            next;
        } else {
            print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
                " numbering scheme\n\n";
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        my $dest = "$master_dir"."/".$HeavyChain;
        getsingleChainAntibody ( $ABChainsRef, \$pdb_id, $hapten, 1, $numbering,
                                 $file_path, $chainsHashRef);
        movePDBs ($dir, $dest, $pdb_id);
	push ( @Hv, $pdb_id );
	$count_Hv++;
	#next;
    }
    #########################################################################
    if ( ( $antigen_c eq "" ) and ( $heavy_c eq "" ) and 
             ( $light_c eq 'Light' ) ) {
        my ($chainsHashRef, $ABChainsRef) =
            processSingleChainAntibody ($file_path, $LOG, $SUMMARY, $pdb_id);
               
        my @antibodyChain = map {split ("", $_)} @{$ABChainsRef};
        print {$ABCHAIN} $pdb_id,":",join (",", @antibodyChain), "\n";

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        eval { antibody_number ( $ABChainsRef, $nsch ); 1; };
	if ( $@ ) {
	    foreach my $mLOG ($LOG, $MASTER_LOG) {
		print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
		    $@ . "Program exited\n\n";
	    }
	    
	    push ( @kabat_failed, $pdb_id );
	    $count_kabatnum_error++;
	    chdir '..';
	    next;
	} else {
	    print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
		" numbering scheme\n\n";
	}
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        my $dest = "$master_dir"."/".$LightChain;
        getsingleChainAntibody ( $ABChainsRef, \$pdb_id, $hapten, 1, $numbering,
                                 $file_path, $chainsHashRef);
        movePDBs ($dir, $dest, $pdb_id);
	push ( @Lg, $pdb_id );
	$count_Lg++;
        next;
    }
    #########################################################################
    if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq "" ) and 
             ( $light_c eq "" ) ) {
        push ( @Ant, $pdb_id );
        $count_Ant++;
        next;
    }
    #########################################################################    
    # 1st If-check start
    # This if block processes the antibodies with only Light and Heavy chains
    if ( ( $antigen_c eq "" ) and ( $heavy_c eq 'Heavy' ) and 
             ( $light_c eq 'Light' ) ) {
	# This subroutine wraps up several other subroutines
	my ($antigen, $antigen_IDs, $hash_keysRef, $chainsHashRef) = 
	    processAntibody($file_path, $LOG, $SUMMARY, $nsch, $pdb_id);
        # Writing antibody Chain IDs on a file for future record
        # To split the paired chain IDs into a single array
        my @antibodyChain = map {split ("", $_)} @{$hash_keysRef};
        print {$ABCHAIN} $pdb_id,":",join (",", @antibodyChain), "\n";

        # Kabatnum is an external program to number antibodies according to 
	# standard numbering scheme - For some antibodies, it does not work and 
	# gives an error. On Failure, the subroutine, antibody_number return an 
	# error by Croak which is being printed on Master log and protein log 
	# file.    
	eval { antibody_number ( $hash_keysRef, $nsch ); 1; };
	if ( $@ ) {
	    foreach my $mLOG ($LOG, $MASTER_LOG) {
		print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
		    $@ . "Program exited\n\n";
	    }
            push ( @kabat_failed, $pdb_id );
	    $count_kabatnum_error++;
	    chdir '..';
	    next;
	} else {
	    print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
		" numbering scheme\n\n";
	}
	
	# Forming antibody complexes and renaming them depending on number of 
	# complexes. Like, 1AFV having 2 complexes will be named as 
	# 1AFV_1 and 1AFV_2
	make_antibody_complex ( $hash_keysRef, \$pdb_id, $hapten, 1, $numbering,
                                $file_path, $chainsHashRef);

        # Moving complexes from individual directory to a separate folder to 
	# store all antibody complexes
	my $dest = $master_dir."/". $FreeAntibodies;
	movePDBs ($dir, $dest, $pdb_id); 
	print {$LOG} "Antibody has been moved to the folder of antibody" . 
	    "complexes\n\n";
    
	push ( @FreeAntibody, $pdb_id );
	$count_FreeAntibody++;

    }                           # 1st If-check ends here

    ###########################################################################################
    # 2nd If-check start                                                          
    # This if block processes the antibodies with Light, Heavy and antigen chains 
    if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq 'Heavy' ) and 
             ( $light_c eq 'Light' ) ) {
        my $ab = "LightHeavy";
        
	my $dest = "$master_dir"."/".$AntibodyAntigen;

	# This subroutine wraps up several other subroutines 
	my ($antigen, $antigen_IDs, $hash_keysRef, $chainsHashRef) =
	    processAntibody($file_path, $LOG, $SUMMARY, $nsch, $pdb_id);

        #        print "TEST: ", Dumper (\%{$chainsHashRef});
        #        exit;
        
        
        # Writing Antigen and antibody Chain IDs on a file for future record
        print {$AGCHAIN} $pdb_id,":",join (",", @{$antigen_IDs}), "\n";
        # To split the paired chain IDs into a single array
        my @antibodyChain = map{split ("", $_)} @{$hash_keysRef};
        print {$ABCHAIN} $pdb_id,":",join (",", @antibodyChain), "\n";
           
        # Kabatnum is an external program to number antibodies according to 
        # standard numbering scheme - For some antibodies, it does not work and 
        # gives an error. On Failure, the subroutine, antibody_number return an 
        # error by Croak which is being printed on Master log and protein log file. 
   
	eval { antibody_number ( $hash_keysRef, $nsch ); 1; };
	if ( $@ ) {
	    foreach my $mLOG ($LOG, $MASTER_LOG) {
		print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
		    $@ . "Program exited\n\n";
	    }
	    
	    push ( @kabat_failed, $pdb_id );
	    $count_kabatnum_error++;
	    chdir '..';
	    next;
	} else {
	    print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
		" numbering scheme\n\n";
	}
	
	# The antibody with Hapten is a complex of numbered antibody and hapten 
	# The molecule is hapten and non protein antigen
	if ( ($hapten) and ($flag)) {
	    my $dest = "$master_dir"."/".$AntibodyHapten;
	    processHapten($file_path, $hash_keysRef);
	    make_antibody_complex ( $hash_keysRef, \$pdb_id, $hapten, 1, $numbering,
                                    $file_path, $chainsHashRef);
	    movePDBs ($dir, $dest, $pdb_id);
	    push (@AbAghapten, $pdb_id); 
	    $count_Abhapten++;
	    next; 
	}

        eval { extract_CDR ( $hash_keysRef, $ab ); 1; };
        if ( $@ ) {
            foreach my $mLOG ( $LOG, $MASTER_LOG ) {
                print {$mLOG} "Can'nt extract CDRs for: $pdb_id\n" . $@ .
                    "Program exited\n\n";
            }
            push ( @CDR_failed, $pdb_id );
            $count_CDR_error++;
            chdir '..';
            next;
        } else {
            print {$LOG} "CDRs for $pdb_id has been extracted from each light/Heavy".
                "chain\n\n";
        }
      
    
	# Wrapper subroutine with several sub routines from pdbWrap
        my ($nonAg, $biAntigen);
        ($nonAg, $biAntigen) = processAntibodyAntigen ($antigen, $antigen_IDs, $hash_keysRef,
				$LOG, $SUMMARY, $pdb_id, $destABAG, $destAB, $dir, $chainsHashRef,
                                $file_path, $numbering);
                
        if ( $nonAg) {
            push (@nonAg, $pdb_id );
            $countNonAg++;
            chdir '..';
            next;
             
        }

        if ( $biAntigen ) {
            push (@biAg, $pdb_id );
            $countBiAg++;
        }
        
        # To deal with DNA/RNA antigens and move them in Hapten folder (Non-protein)
        if ( grep {/DNA|RNA/} @antigenChainType) {
            my $dest = "$master_dir"."/".$AntibodyHapten;
            movePDBs ($dir, $dest, $pdb_id);
            push (@AbAgDRNA, $pdb_id);
            $count_AbDNA_RNA++;
        } else {
            # Moving complexes from individual directory to a separate folder to 
            # store all antibody/antigen complexes 
            #            my $dest = "$master_dir"."/".$AntibodyAntigen;
            #            movePDBs ($dir, $dest, $pdb_id);
            push ( @AbAg_complex, $pdb_id );
            $count_AbAg++;
        }
            
    }                           # 2nd If-check ends here
    
    chdir '..';
    #  last; 
}   # Main For loop ends here
#chdir '..';
print "$countBiAg proteins were found bound with 2 antigens\n\n\n";

print "$count_FreeAntibody proteins were found as only antibody complexes\n";
print "$count_AbAg proteins has been successfully processed " .
    "for antibody/antigen complexes\n";
print "$countNonAg proteins has been successfully processed " .
    "for non antigen complexes\n";
print "$count_Abhapten free antibodies were found complexed with Hapten\n";
print "$count_AbDNA_RNA free antibodies were found complexed with RNA and DNA\n";

print "$count_Lg proteins were found as Light chains\n";
print "$count_LgAnt proteins were found as Light-Antigen complexes\n";
print "$count_Lhapten light chain antibodies were found as Light-Hapten\n";
print "$count_LDNA_RNA light chain antibodies were found complexed with RNA and DNA\n";

print "$count_Hv proteins were found as Heavy chains\n";
print "$count_HvAnt proteins were found as Heavy-Antigen complexes\n";
print "$count_Hhapten heavy chain antibodies were found as Heavy-Hapten\n";
print "$count_HDNA_RNA heavy chain antibodies were found complexed with RNA and DNA\n";

print "$count_Ant proteins were found as Antigen chains (Fc Fragments)\n";
print "$count_kabatnum_error proteins failed to number by program Kabatnum\n";
print "$count_superseded proteins were found as Superseded\n";
print "$count_CDR_error proteins were failed for CDR extraction\n";
print "$count_idabError proteins were failed to identify by idabchain program\n";

#######################

chdir '..';

# open files for recording the PDB codes for processed and unprocessed Data 
open (my $NONAGCOMPLEX, '>NonAg_complexes.list') or
    die "Can not open file\n";
open (my $BIAGCOM, '>BiAntigens_complexes.list') or
    die "Can not open file\n";

open (my $ABAGCOMPLEX, '>AbAg_complexes.list') or
    die "Can not open file\n";
open (my $FREE_ANTIBODY, '>Antibody_complexes.list') or
    die "Can not open file\n";
open (my $ABAG_HAPTEN, '>Antibody_Hapten.list') or
    die "Can not open file\n";
open (my $ABAG_DRNA, '>Antibody_DRNA.list') or
    die "Can not open file\n";


open (my $LG, '>Light_Only.list') or
    die "Can not open file\n";
open (my $LG_AG, '>Light_Antigen_complexes.list') or
    die "Can not open file\n";
open (my $LG_HAPTEN, '>Light_Hapten.list') or
    die "Can not open file\n";
open (my $LG_DRNA, '>Light_DRNA.list') or
    die "Can not open file\n";


open (my $HV, '>Heavy_Only.list') or
    die "Can not open file\n";
open (my $HV_AG, '>Heavy_Antigen_complexes.list') or
    die "Can not open file\n";
open (my $HV_HAPTEN, '>Heavy_Hapten.list') or
    die"Can not open file\n";
open (my $HV_DRNA, '>Heavy_DRNA.list') or
    die"Can not open file\n";


open (my $FC, '>FC_Fragments.list') or
    die "Can not open file\n";
open (my $FAILED, '>Kabat_Failed.list') or
    die "Can not open file\n";
open (my $SUPERSEDED, '>Superseded.list') or
    die "Can not open file\n";
open (my $CDR_ERROR, '>CDR_Error.list') or
    die "Can not open file\n";
open (my $IDAB_ERROR, '>IDAB_Error.list') or
    die "Can not open file\n";

print {$MASTER_LOG} ">$count_AbAg proteins has been successfully " . 
    "processed for antibody/antigen complexes\n";
print {$MASTER_LOG} "@AbAg_complex\n";
print {$ABAGCOMPLEX} join("\n", @AbAg_complex);


print {$MASTER_LOG} ">$countNonAg proteins has been successfully " .
    "processed for non antigen complexes\n";
print {$MASTER_LOG} "@nonAg\n";
print {$NONAGCOMPLEX} join("\n", @nonAg);

print {$MASTER_LOG} ">$countBiAg proteins has been successfully " .
    "processed for bi-antigen complexes\n";
print {$MASTER_LOG} "@biAg\n";
print {$BIAGCOM} join("\n", @biAg);


print {$MASTER_LOG} ">$count_FreeAntibody proteins were found as only antibody " .
    "complaxes\n";
print {$MASTER_LOG} "@FreeAntibody\n";
print {$FREE_ANTIBODY} join ("\n", @FreeAntibody), "\n";

print {$MASTER_LOG} "$count_Abhapten free antibodies were found complexed with Hapten\n";
print {$MASTER_LOG} "@AbAghapten";
print {$ABAG_HAPTEN} join ("\n", @AbAghapten), "\n";

print {$MASTER_LOG} "$count_AbDNA_RNA free antibodies were found complexed with RNA and DNA\n";
print {$MASTER_LOG} "@AbAgDRNA";
print {$ABAG_DRNA} join ("\n", @AbAgDRNA), "\n";



print {$MASTER_LOG} ">$count_Lg proteins were found as Light chains\n";
print {$MASTER_LOG} "@Lg\n";
print {$LG} join ("\n", @Lg), "\n";

print {$MASTER_LOG} ">$count_LgAnt proteins were found as Light-Antigen complexes\n";
print {$MASTER_LOG} "@LgAnt\n";
print {$LG_AG} join ("\n", @LgAnt), "\n";

print {$MASTER_LOG} "$count_Lhapten light chain antibodies were found complexed with Hapten\n";
print {$MASTER_LOG} "@LgHapten";
print {$LG_HAPTEN} join ("\n", @LgHapten), "\n";

print {$MASTER_LOG} "$count_LDNA_RNA light chain antibodies were found complexed with RNA and DNA\n";
print {$MASTER_LOG} "@LgDRNA";
print {$LG_DRNA} join ("\n", @LgDRNA), "\n";


print {$MASTER_LOG} ">$count_Hv proteins were found as Heavy chains\n";
print {$MASTER_LOG} "@Hv\n";
print {$HV} join ("\n", @Hv), "\n";

print {$MASTER_LOG} ">$count_HvAnt proteins were found as Heavy-Antigen complexes\n";
print {$MASTER_LOG} "@HvAnt\n";
print {$HV_AG} join("\n", @HvAnt), "\n";

print {$MASTER_LOG} "$count_Hhapten heavy chain antibodies were found complexed with Hapten\n";
print {$MASTER_LOG} "@HvHapten";
print {$HV_HAPTEN} join ("\n", @HvHapten), "\n";

print {$MASTER_LOG} "$count_HDNA_RNA heavy chain antibodies were found complexed with RNA and DNA\
n";
print {$MASTER_LOG} "@HvDRNA";
print {$HV_DRNA} join ("\n", @HvDRNA), "\n";



print {$MASTER_LOG} ">$count_Ant proteins were found as Antigen chains (Fc Fragments)\n";
print {$MASTER_LOG} "@Ant\n";
print {$FC} join ("\n", @Ant), "\n";

print {$MASTER_LOG} ">$count_kabatnum_error proteins failed to number by" .
    "program Kabatnum\n";
print {$MASTER_LOG} "@kabat_failed\n"; 
print {$FAILED} join ("\n", @kabat_failed), "\n"; 

print {$MASTER_LOG} ">$count_superseded proteins were found as superseded\n";
print {$MASTER_LOG} "@superseded\n";
print {$SUPERSEDED} join ("\n", @superseded), "\n";

print {$MASTER_LOG} ">$count_CDR_error proteins were failed for CDR extraction".
    "\n";
print {$MASTER_LOG} "@CDR_failed\n";
print {$CDR_ERROR} join ("\n", @CDR_failed), "\n";

print {$MASTER_LOG} ">$count_idabError proteins were failed to identify by idabchain program".
    "\n";
print {$MASTER_LOG} "@idab_failed\n";
print {$IDAB_ERROR} join ("\n", @idab_failed), "\n";


my $count_AbHapten = $count_Abhapten+$count_AbDNA_RNA;
my $totalAB = $count_AbAg+$count_FreeAntibody+$count_Abhapten+$count_AbDNA_RNA, "\n";

my $count_LgHapten = $count_Lhapten+$count_LDNA_RNA; 
my $totalLG = $count_Lg+$count_LgAnt+$count_LgHapten;

my $count_HvHapten = $count_Hhapten+$count_HDNA_RNA;
my $totalHV = $count_Hv+$count_HvAnt+$count_HvHapten;
print {$MASTER_LOG} "Bi-Antigen=$countBiAg\n";
print {$MASTER_LOG} "NonAntigen=$countNonAg\n";
print {$MASTER_LOG} "AntibodyAntigen=$count_AbAg\n";
print {$MASTER_LOG} "Free_Antibody=$count_FreeAntibody\n";
print {$MASTER_LOG} "AB-Hapten=$count_AbHapten\n";
print {$MASTER_LOG} "CompleteAntibody_Dataset=$totalAB\n";

print {$MASTER_LOG} "Bence-Jones=$count_Lg\n";
print {$MASTER_LOG} "Light-Antigen=$count_LgAnt\n";
print {$MASTER_LOG} "Light-Hapten=$count_LgHapten\n";
print {$MASTER_LOG} "CompleteLight_Dataset=$totalLG\n";

print {$MASTER_LOG} "Camelids=$count_Hv\n";
print {$MASTER_LOG} "Heavy-Antigen=$count_HvAnt\n";
print {$MASTER_LOG} "Heavy-Hapten=$count_HvHapten\n";
print {$MASTER_LOG} "CompleteHeavy_Dataset=$totalHV\n";

print {$MASTER_LOG} "Fc=$count_Ant\n";
print {$MASTER_LOG} "Failed=$count_kabatnum_error\n";
print {$MASTER_LOG} "Superseded=$count_superseded\n";
print {$MASTER_LOG} "CDR-Error=$count_CDR_error\n";

print {$MASTER_LOG} "IDAB-Error=$count_idabError\n";


print "See masterlog.log for details...\n";

=head1 NAME

get_antibody_complex.pl 

=head1 SYNOPSIS

To use this program, type on command line, 
./processAntibodyPDBs.pl [-k -c -a] input file

=head1 DESCRIPTION

This script reads a list of antibody PDBs (obtained from SACs)  
and process it. The procssing includes 1) antibody numbering    
according to given numbering scheme 2) splitting PDB file into  
multiple PDB files for antibodies present in individual PDB     
3) Moves the processed antibody/antigen complexes, antibody     
only and antibody/hapten complxes into separate directories 
 
=head1 AUTHOR

Saba Ferdous (ucbterd@acrm19)

=head1 COPYRIGHT AND LICENSE
Copyright (C) 2014 by Saba Ferdous

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS
None reported... yet.

=cut
