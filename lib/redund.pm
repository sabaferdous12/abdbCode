package redund;
use Carp;
use strict;
#use warnings;
use Data::Dumper;

use Exporter qw (import);
our @EXPORT_OK = qw (getAntibodySeq readDir getFileData aa3to1 splitString 
getEndsSeq compareStarts compareEnds checkSubstring printClusters 
isIdenticalEnds);
# ************* getAntibodySeq *****************                               
# Description: Reads a PDB (Antibody-Antigen Complex) file and extracts 
#              sequence from structure
# Inputs: A PDB file
# Outputs: 2 strings each with sequence of light and heavy chain
# Subroutine call/Testing: $ab1Light, $ab1Heavy) = 
#                          getAntibodySeq($antibodyList[$ab1]);
# Date: 22 April 2015       
# Author: Saba
sub getAntibodySeq
{
    my ($inputPro) = @_;
    my @protein = getFileData ($inputPro);
    my (@light,@heavy);
    # Reads each line of PDB file
    foreach my $line (@protein)
    {
	if ($line =~ /^ATOM*/)
	{	
	    my @pro = split (' ', $line);
	    # Reads CAs of light chain
	    if ($pro[4] eq 'L' and $pro[2] eq 'CA')
	    {
		my $lg = aa3to1($pro[3]);
		push (@light, $lg);
	    }
	    # Reads CAs of heavy chain
	    elsif ($pro[4] eq 'H' and $pro[2] eq 'CA')
	    {
		my $hv = aa3to1($pro[3]);
		push (@heavy, $hv);
	    }
	}
    }
    my $seqLight = join('',@light);
    my $seqHeavy = join('',@heavy);
    return ($seqLight, $seqHeavy);
}

# ************* readDir *****************                               
# Description: Reads all the PDB files in the directory and return as a list
# Inputs: Directory name
# Outputs: An array with the name of all PDB files
# Subroutine call/Testing: my @antibodyList = readDir ($dir);
# Date: 22 April 2015       
# Author: Saba
sub readDir
{
    my ($dir) = @_;
    my @list;
    opendir (DIR, "$dir") or die "Can't open directory: $!\n";
    while ( my $file = readdir(DIR))
    {
	next unless (-f "$dir/$file"); ## Reads only files
	next unless ($file =~ m/\.pdb$/); ## Reads files with .pdb extension
	push(@list, $file);
    }
    closedir (DIR);
    return @list;
}

# ************* getFileData *****************                               
# Description: Reads the contents of a file into an array
# Inputs: A file name
# Outputs: An array with each line of file as an element of an array
# Subroutine call/Testing: my @fileData = getFileData ($fileName); 
# Date: 22 April 2015       
# Author: Saba
sub getFileData
{
    my($filename) = @_;
    chomp $filename;
    my @filedata = ();
    unless( open(FILE_DATA, $filename) )
    {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    @filedata = <FILE_DATA>;
    close FILE_DATA;
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

# ************* splitString *****************                               
# Description: Split the given string in to 3 sections
# Inputs: A string of sequence
# Outputs: 3 strings, one for each section
# Subroutine call/Testing: ($ab1Lstart, $ab1Lmid, $ab1Lend)
#                             = splitString($ab1Light);
# Date: 22 April 2015       
# Author: Saba
sub splitString
{
    my ($inputStr) = @_;
    my ($startSec, $midSec, $endSec, $tempLen, $tempStr);

    $startSec = substr($inputStr, 0, 5); # First section with length 5
    $tempStr = substr($inputStr, 5); # Rest of string excluding first 5 A.A
    $midSec = substr($tempStr, 0, length($tempStr) - 5); # Middle Sec Sequence
    $tempLen = length($startSec) + length($midSec);
    $endSec = substr($inputStr, $tempLen); # Last section sequence
    
    return ($startSec, $midSec, $endSec);
}

# ************* getEndsSeq *****************                               
# Description: Finds sequence of end sections
# Inputs: 2 strings 1) a reference string 2) sub string to be compared
# Outputs: 2 strings with the sequence of end sections and an index of 
#          substring
# Subroutine call/Testing: my ($ab2Lstart, $ab2Lend, $index_s) = 
#                             getEndsSeq($ab2Light, $ab1Lmid);
# Date: 22 April 2015       
# Author: Saba
sub getEndsSeq
{
    my ($fullStr, $subStr) = @_;
    my ($startSec, $endSec, $indexE, $indexS);

    $indexS = index($fullStr, $subStr); # string index where substring starts  
    $startSec = substr($fullStr, 0, $indexS); # sequence of start section 
    $indexE = length($startSec) + length($subStr); # index where substring ends
    $endSec = substr($fullStr, $indexE); # sequence of end section

    return ($startSec, $endSec, $indexS);
}

# ************* compareStarts *****************                               
# Description: This sub routine compares two strings starting starting from 
#              the end of strings and returns 1 or 0 on match or mismatch
# Inputs: 2 strings and the index of second string
# Outputs: Returns either 0 or 1 for match or mismatch
# Subroutine call/Testing: my $statusStartL = 
#                           compareStarts($ab1Lstart, $ab2Lstart, $index_s);
# Date: 22 April 2015       
# Author: Saba
sub compareStarts
{
    my ($str1,$str2, $index) = @_;
    my $count = 0;
    my $match = '';
    for (my $i = 4; $i>=0; $i--) # Needs reverse loop to iterate from end
    {
	if ($str1 && $str2)
	{
	    if (substr($str1, $i, 1) eq 
		substr($str2, $index-1, 1) )
	    {
		$match .= substr($str2, $index-1, 1); # Match string
		$count++;
		$index--;
	    }
	    else
	    {
		next;
	    }
	}
	
    }
    my $revMatch = reverse($match); # Get real string 
    my $status = checkSubstring($str1, $revMatch);
    return $status;
}

# ************* compareEnds *****************                               
# Description: This sub routine compares two strings starting starting from 
#              the start of strings and returns 1 or 0 on match or mismatch
# Inputs: 2 strings 
# Outputs: Returns either 0 or 1 for match or mismatch
# Subroutine call/Testing:  my $statusEndL = compareEnds($ab1Lend, $ab2Lend);
# Date: 22 April 2015       
# Author: Saba
sub compareEnds
{
    my ($str1, $str2) = @_;
    my $match = '';
    for (my $i = 0; $i<=4; $i++)
    {
	if(!$str2)
	{
	    if (substr($str1, $i, 1) eq 
		substr($str2, $i, 1))
	    {
		$match .= substr($str2, $i, 1);
	    }
	    else
	    {
		next;
	    }
	}
	else
	{
	    next;
	}
    }
    my $status = checkSubstring($str1, $match);
    return $status;
}

# ************* checkSubstring *****************                               
# Description: Checks if first given string contains the second given string
# Inputs: 2 strings to be checked
# Outputs: Returns either 0 or 1 for match or mismatch
# Subroutine call/Testing:  my $status = checkSubstring($str1, $match);
# Date: 22 April 2015       
# Author: Saba
sub checkSubstring
{
    my ($str1, $str2) = @_;
    my $status;
    if (index($str1, $str2) != -1)
    {
	$status = 1;
    }
    else
    {
	$status= 0;
    }
    return $status;
}

# ************* isIdentical *****************                               
# Description: This sub routine compares the starts and ends of given strings.
#              It first calculates the start and ends of second string and then
#              compares these with the start and ends of first string (inputs) 
# Inputs: 4 strings 
# Outputs: Returns status (match or mismatch) for start and end sections
# Subroutine call/Testing:  ($statusStartL, $statusEndL) =
#                      isIdentical ($ab1Lmid, $ab2Light, $ab1Lstart, $ab1Lend);
# Date: 22 April 2015       
# Author: Saba
sub isIdenticalEnds
{
    my ($ab1Lmid, $ab2Light, $ab1Lstart, $ab1Lend) = @_; 
    # Obtain start and ends of ab2                                      
    my ($ab2Lstart, $ab2Lend, $index_s) = getEndsSeq($ab2Light, $ab1Lmid);
    # Compare the starts and ends and returns 0 or 1                    
    my $statusStart = compareStarts($ab1Lstart, $ab2Lstart, $index_s);
    my $statusEnd = compareEnds($ab1Lend, $ab2Lend);
    
    return $statusStart, $statusEnd; 
}

# ************* printClusters *****************                               
# Description: print the contents of hashes as redundant clusters
# Inputs: 2 hashes; 1) hash with cluster of redundant PDBs, 2) hash with unique
#                   PDBs having no redundant PDB
# Outputs: A text file with information of reduntant clusters
# Subroutine call/Testing: printClusters (\%redundantClusters, \%unique, $OUT);
# Date: 22 April 2015       
# Author: Saba
sub printClusters
{
    my ($redundantClustersRef, $uniqueRef, $OUT) = @_;
    # Dereferening
    my %redundantClusters = % {$redundantClustersRef};
    my %unique = % {$uniqueRef};

    my ($completeClusterRef, @completeCluster);
    
    # Print the redundadant cluster hash
    foreach my $key (sort { $a <=> $b } keys %redundantClusters)
    {
        #Since hash value is an array - So it needs to dereferenced as below:
        #The key and value are assigned to an array                            
        my $completeClusterRef = 
	    [$key, join(" ", @{$redundantClusters{$key}}) ];
	my @completeCluster = split ("\n", "@{$completeClusterRef}");
	
        # Each file is with.pdb extension, To avoid that in final output file, 
        # REGEX is used as below:                                            
        foreach my $line (@completeCluster)
        {
            $line =~ s/.pdb/,/g;
            print $OUT "$line\n";
        }
    }
    
    # print the hash with PDBs with no redundant PDB 
    foreach my $uniq (sort {$a<=>$b} keys %unique)
    {
        $uniq =~ s/.pdb//g;
        print $OUT "$uniq\n";
    }
}

