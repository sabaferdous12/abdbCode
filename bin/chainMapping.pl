use strict;
use warnings; 

# This script reads processed antibody complexes and maps 
# chain information with the original chain IDs from real 
# PDB file and writes a flat file for all the complexes in 
# the given directory.
# It uses header information from the given processed PDB 
# file 

# Usage: perl chainMapping.pl 
 
my $dir = "."; # Reads PDB diles 
my @files = readDirPDB($dir);
@files = sort @files;


foreach my $f (@files) # Reads each PDB file
{
    open (my $IN, '<', $f) or die "Can not open file !!!\n";

    my @filename = split (/\./, $f);
    
    my ($l, $b, $h, $d, $a);
    my (@Ls, @As);
    
    while (my $line = <$IN>) # Reading each line of the given file
    {
	# Checking for L chains
	if ($line =~ /REMARK 950 CHAIN L/)
	{
	    my @lineCols = split (" ", $line);
	  	    
	    $l = $lineCols[4];
	    $b = $lineCols[5];
	    push @Ls, $lineCols[5]; # For multiple L chains in a file
	}
	
	# Checking for H chains
	if ($line =~ /REMARK 950 CHAIN H/)
	{
	    my @lineCols = split (" ", $line);

	    $h = $lineCols[4];
	    $d = $lineCols[5]; # No H chaain dimers in a single file
	}

	if ($line =~ /REMARK 950 CHAIN A/)
	{
	    my @lineCols = split (" ", $line);
           
            $a = $lineCols[3];
            # $f = $lineCols[5];
            push @As, $lineCols[5]; # For multiple A chains in a file
	}

#    last;
}

# ******* Printing *******

    if ($a) # If bound with antigen
    {

	if ( ($l) and ($h)) # If it is L and H chain
	{ 
	    if (!$As[1]) # multiple antigen chains in a file
	    {
		print "$filename[0],$l:$b,$h:$d,$a:$As[0]\n";
	    }
	    else
	    {
		print "$filename[0],$l:$b,$h:$d,$a:$As[0],$a:$As[1]\n";
	    }

	} # L and H check end

	elsif ( ($l) and (!$h) ) # L only 
	{
	    if (!$As[1])
	    {
	    print "$filename[0],$l:$b,$a:$As[0]\n";
	    }

	    else
	    {
		print "$filename[0],$l:$b,$a:$As[0],$b:$As[1]\n";
	    }
   
	} # L only check ends 

	if ( ($h) and (!$l) ) # H only 
	{

            if (!$As[1])
            {
		print "$filename[0],$h:$d,$a:$As[0]\n";
	    }
            else
            {
                print "$filename[0],$h:$d,$a:$As[0],$a:$As[1]\n";
            }

	} # H only check ends




    } # bound antigen check ends 


    else # Not bound with antigen
    {
	if (($l) and ($h)) # L and H
	{
	print "$filename[0],$l:$b,$h:$d\n";
	} # L and H check ends

	if ( ($l) and (!$h)) # L only
	{
	    if (!$Ls[1]) # For multiple L chains in a file
	    {
		print "$filename[0],$l:$Ls[0]\n";
	    }
	    else
	    {
		print "$filename[0],$l:$Ls[0],$l:$Ls[1]\n";
	    }

	} # L only check ends

	if ( ($h) and (!$l) ) # H only
	{
	    print "$filename[0],$h:$d\n";
	} # H only check ends

    } # not bound antigen check ends

} # Each File check  


# ***** METHODS ********

sub readDirPDB
{
    my ($dir) = @_;
    my @list;
    opendir (DIR, "$dir") or die "Can't open directory: $!\n";
    while ( my $file = readdir(DIR)){
	next unless (-f "$dir/$file"); # Reads only files
	next unless ($file =~ m/\.pdb$/); # Reads files ending with pdb
	push(@list, $file);
    }
    closedir (DIR);
        
    return @list;
}
 
