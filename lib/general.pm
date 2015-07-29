package general; 

use Exporter qw (import);
our @EXPORT_OK = qw (getPdbCodeListFromXmlFile readDirPDB); 

sub getPdbCodeListFromXmlFile
{
    my ( $inputFile ) = @_;
    open ( my $IN, "<$inputFile" ) or die "Can not open $inputFile\n";
    my @pdbLines = 
grep /pdb/, <$IN>; # Read file in array / grep line with word pdb 
    my @pdbCodes;

    foreach my $line ( @pdbLines )
    {
my ( $start, $code) = split ( "\"", $line );
push ( @pdbCodes, $code );
    }

    close $IN;

    return @pdbCodes;

}


# Inputs: path of Directory
# Outputs: Array with all files with pdb extension only
# Subroutine call: read_dir($dir)
# Testing: my $dir_files = read_dir($dir)
# Date: 02 April 2014      
# Author: Saba 
sub readDirPDB{
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

# Reads all the files of given extension from directory 
# and returns them as a list in an array
sub readSpecificFileFromDir
{
    my ($dir, $ext) = @_;
    my @list;
    
    opendir (DIR, "$dir") or die "Can't open directory: $!\n";
    while ( my $file = readdir(DIR))
    {
	next unless (-f "$dir/$file"); ## Reads only files                                                         
	next unless ($file =~ m/\.$ext$/); ## Reads files ending with pdb                                           
	push(@list, $file);
    }
    closedir (DIR);
    return @list;
}
