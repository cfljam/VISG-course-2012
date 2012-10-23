#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

huge-merge.pl - Merge the results of multiple huge-sort generated files into a single sorted file. 

=head1 SYNOPSIS

huge-merge.pl output-directory 

=head1 DESCRIPTION

Combine the sorted bigram files generated by huge-sort.pl efficiently.

This program is used internally by huge-count.pl. 

=head1 USGAE

huge-merge.pl [OPTIONS] SOURCEDIR

=head1 INPUT

=head2 Required Arguments:

=head3 SOURCEDIR

Input to huge-merge.pl should be a single flat directory containing 
multiple plain text files generated by huge-sort.pl. The result file,
merge.* (* is a number, the final result file has the maximum number), 
is in the source directory.  

=head2 Optional Arguments:

=head4 --keep  

Switches ON the --keep option will keep all the intermediate merging 
files.

=head3 Other Options:

=head4 --help

Displays the help information.

=head4 --version

Displays the version information.

=head3 BUGS

There is a limitation in huge-merge.pl. When the size of the
corpus is very large (>16G)  and the some of the terms of the
bigrams is very long (>30 chars), the program could run out of
memory at huge-merge.pl step. This is because huge-merge use
two hashes to count the frequencies of the first and second
term of the bigrams. These two hashes could use up the memory
with the increase of the length of the terms and the increase
of the number of the terms. If just for normal text, terms
are within limited length and numbers, the software won't
use up the memory.

=head1 AUTHOR

Ying Liu, University of Minnesota, Twin Cities.
liux0395 at umn.edu

Ted Pedersen, University of Minnesota, Duluth.
tpederse at umn.edu

=head1 COPYRIGHT

Copyright (C) 2009-2011, Ying Liu and Ted Pedersen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=cut


###################################################################################

use Getopt::Long;

# first check if no commandline options have been provided... in which case
# print out the usage notes!
if ( $#ARGV == -1 )
{
    &minimalUsageNotes();
    exit;
}

# now get the options!
GetOptions( "keep", "version", "help" );

if ( defined $opt_keep)    { $opt_keep = 1; }
else                          { $opt_keep = 0; }


# if help has been requested, print out help!
if ( defined $opt_help )
{
    $opt_help = 1;
    &showHelp();
    exit;
}

# if version has been requested, show version!
if ( defined $opt_version )
{
    $opt_version = 1;
    &showVersion();
    exit;
}

# readin the file names under the folder
$dir = $ARGV[0]; 
if ( !($dir) )
{
    print STDERR "No input folder name supplied.\n";
    askHelp();
    exit;
}

opendir(BIN, $dir) or die "Can't open $dir: $!";
my @files;
while( defined ($file = readdir BIN) ) 
{
	if($file=~/sorted$/)
	{
		push (@files, "$dir/$file") if -T "$dir/$file";
	}
}
closedir(BIN);
#print "merge files: @files\n";

my $i = 0;
my $bigramTotal = 0;
my %f1_w1;
my %f1_w2;
my %f2_w1;
my %f2_w2;
my %w1; 
my %w2;

# without split the bigrams list, only one sorted file
if(@files==1)
{
	my $file1 = shift @files;

	open(FILE1, "<$file1") or die("Error: cannot open file '$file1'\n");		

	my $merge = "$dir" . "/merge." . "$i";
	open(MERGE, ">$merge") or die("Error: cannot open file '$merge'\n");		

	while(my $line= <FILE1>)
	{
		chomp($line);
		my @s = split ('<>', $line);
		my @fre = split (' ', $s[2]); 
		$bigramTotal += $fre[0];
	}

	print MERGE "$bigramTotal\n";
	seek FILE1, 0, 0;
	while(my $line= <FILE1>)
	{
		chomp($line);
		print MERGE "$line\n";
	}
	close MERGE;
	close FILE1;

	# remove the unsorted duplicated bigrams
    if ($opt_keep == 0)
    {
		#print "remove $file1 $file2\n";
		system ("rm $file1");
	}

}


# split the bigrams list, more than one sorted file
while(@files>1)
{
	$bigramTotal = 0;

	$i++;
	%f1_w1 = ( );
	%f1_w2 = ( );
	%f2_w1 = ( );
	%f2_w2 = ( );
	%w1 = ( );
	%w2 = ( );

	my $file1 = shift @files;
	my $file2 = shift @files;
	open(FILE1, "<$file1") or die("Error: cannot open file '$file1'\n");		
	open(FILE2, "<$file2") or die("Error: cannot open file '$file2'\n");		
	#print "merge files: $file1 $file2\n";

	my $temp = "$dir" . "/temp." . "$i";
	open(TEMP, ">$temp") or die("Error: cannot open file '$temp'\n");		

	my $merge = "$dir" . "/merge." . "$i";
	open(MERGE, ">$merge") or die("Error: cannot open file '$merge'\n");		
		
	my $flag = 0;


#print "merge: $i\n";


	while ( )
	{
		if (!eof(FILE1) and !eof(FILE2))
		{
			if ($flag == 1) # read file 1
			{
				$line1 = <FILE1>; chop ($line1);
				my @s = split ('<>', $line1);
				$b1[1] = "$s[0]<>$s[1]<>"; # @b1 holds the bigram string
				my @fre = split (' ', $s[2]); # bigram freq, w1 freq and w2 freq in @fre
				$b1[0] = $fre[0];  
				$f1_w1{$s[0]} = $fre[1];
				$f1_w2{$s[1]} = $fre[2];
				$bigramTotal += $fre[0];
			}
			elsif ($flag == 2) # read file 2
			{
				$line2 = <FILE2>; chop ($line2);
				my @s = split ('<>', $line2);
				$b2[1] = "$s[0]<>$s[1]<>";  # @b2 holds the bigram string
				my @fre = split (' ', $s[2]); 
				$b2[0] = $fre[0]; 
				$f2_w1{$s[0]} = $fre[1];
				$f2_w2{$s[1]} = $fre[2];
				$bigramTotal += $fre[0];
			}
			elsif ($flag == 0) # read file 1 and 2
			{
				$line1 = <FILE1>; chop ($line1);
				my @s1 = split ('<>', $line1);
				$b1[1] = "$s1[0]<>$s1[1]<>";
				my @fre1 = split (' ', $s1[2]); 
				$b1[0] = $fre1[0]; 
				$f1_w1{$s1[0]} = $fre1[1];
				$f1_w2{$s1[1]} = $fre1[2];
				$bigramTotal += $fre1[0];

				$line2 = <FILE2>; chop ($line2);
				my @s2 = split ('<>', $line2);
				$b2[1] = "$s2[0]<>$s2[1]<>";
				my @fre2 = split (' ', $s2[2]); 
				$b2[0] = $fre2[0]; 
				$f2_w1{$s2[0]} = $fre2[1];
				$f2_w2{$s2[1]} = $fre2[2];
				$bigramTotal += $fre2[0];
			}

			if ($b1[1] eq $b2[1]) # two string is the same, add their freqs
			{
				my $total = $b1[0] + $b2[0];
				print TEMP "$b1[1]$total\n";
				$flag = 0;
			}
			elsif ($b1[1] gt $b2[1]) # print string 2 
			{
				print TEMP "$b2[1]$b2[0]\n";		
				$flag = 2;	
				print TEMP "$b1[1]$b1[0]\n" if(eof(FILE1)and eof(FILE2));		
			}
			elsif ($b1[1] lt $b2[1]) # print string 1
			{
				print TEMP "$b1[1]$b1[0]\n";		
				$flag = 1;	
				print TEMP "$b2[1]$b2[0]\n" if (eof(FILE2)and eof(FILE1));		
			}
		}	
		elsif (eof(FILE1) and eof(FILE2))
		{
			last;
		}
		else
		{

			# $left_freq is already added in to $bigramTotal
			my $left_string = "";
			my $left_freq = 0;
			my $print_flag = 0;
			if ($flag == 1)
			{
				$left_string = $b2[1];
				$left_freq = $b2[0];

			}

			if ($flag == 2)
			{
				$left_string = $b1[1];
				$left_freq = $b1[0];

			}

			if(!eof(FILE1))
			{
				while($line1 = <FILE1>)
				{
					chop($line1);
					my @s = split ('<>', $line1);
					my @fre = split (' ', $s[2]); 
					$f1_w1{$s[0]} = $fre[1];
					$f1_w2{$s[1]} = $fre[2];

					my $bigram_string = "$s[0]<>$s[1]<>";
					if ($left_string ne "") #$flag = 1 or flag = 2
					{
						if ($print_flag==0)
						{
							if ($bigram_string eq $left_string)
							{
								my $frequency = $fre[0] + $left_freq;
								print TEMP "$bigram_string$frequency\n";
								$print_flag = 1;
							}
				
							if ($bigram_string lt $left_string)
							{
								print TEMP "$bigram_string$fre[0]\n";
							}	

							if ($bigram_string gt $left_string)
							{
								print TEMP "$left_string$left_freq\n"; 
								$print_flag = 1;
								print TEMP "$bigram_string$fre[0]\n";
							}	
						}
						else	
						{
							print TEMP "$bigram_string$fre[0]\n";
						}

						$bigramTotal += $fre[0];
					}
					else #$flag = 0
					{
						print TEMP "$bigram_string$fre[0]\n";
						$bigramTotal += $fre[0];
					}
				}

				if (($print_flag==0) and ($left_string ne ""))
				{
					print TEMP "$left_string$left_freq\n"; 
				}

			}
			elsif(!eof(FILE2))
			{
				while($line2= <FILE2>)
				{
					chop($line2);
					my @s = split ('<>', $line2);
					my @fre = split (' ', $s[2]); 
					$f2_w1{$s[0]} = $fre[1];
					$f2_w2{$s[1]} = $fre[2];

					my $bigram_string = "$s[0]<>$s[1]<>";

					if ($left_string ne "") #$flag = 1 or $flag = 2
					{
						if ($print_flag==0)
						{
							if ($bigram_string eq $left_string)
							{
								my $frequency = $fre[0] + $left_freq;
								print TEMP "$bigram_string$frequency\n";
								$print_flag = 1;
							}
				
							if ($bigram_string lt $left_string)
							{
								print TEMP "$bigram_string$fre[0]\n";
							}	

							if ($bigram_string gt $left_string)
							{
								print TEMP "$left_string$left_freq\n"; 
								$print_flag = 1;
								print TEMP "$bigram_string$fre[0]\n";
							}	
						}
						else	
						{
							print TEMP "$bigram_string$fre[0]\n";
						}


						$bigramTotal += $fre[0];

					}
					else #$flag = 0
					{
						print TEMP "$bigram_string$fre[0]\n";
						$bigramTotal += $fre[0];
					}
				}

				if (($print_flag==0) and ($left_string ne ""))
				{
					print TEMP "$left_string$left_freq\n"; 
				}

			}


		}

	} # end of while (), merge 2 files

	close FILE1;
	close FILE2;
	close TEMP;

	# merge single word w1 and w2 frequency of file1 and file2 
	%w1 = %f1_w1;
	foreach my $key2 (keys %f2_w1)			
	{
		if (exists $w1{$key2})
		{
			$w1{$key2} += $f2_w1{$key2};
		}
		else
		{
			$w1{$key2} = $f2_w1{$key2};
		}
	}

	%w2 = %f1_w2;
	foreach my $key2 (keys %f2_w2)			
	{
		if (exists $w2{$key2})
		{
			$w2{$key2} += $f2_w2{$key2};
		}
		else
		{
			$w2{$key2} = $f2_w2{$key2};
		}
	}


	# TEMP file for hold the merge file before get the 
	# correct freqeuncy of each word of the bigrams
	if (@files==0) # no files in the queue for merge
	{
		print MERGE "$bigramTotal\n";
	}
	
	open(TEMP, "<$temp") or die("Error: cannot open file '$temp'\n");		
	while (my $line = <TEMP>)
	{
		chop ($line);
		my @words = split('<>', $line);	
		print MERGE "$line $w1{$words[0]} $w2{$words[1]} \n"; 
	}
	close TEMP;
	close MERGE;

	system("rm $temp");
	push (@files, $merge);

	# remove the unsorted duplicated bigrams
    if ($opt_keep == 0)
    {
		#print "remove $file1 $file2\n";
		system ("rm $file1");
		system ("rm $file2");
	}

} # end of while (@files>1)





#-----------------------------------------------------------------------------
#                       User Defined Function Definitions
#-----------------------------------------------------------------------------

# function to output a minimal usage note when the user has not provided any
# commandline options
sub minimalUsageNotes
{
    print STDERR "Usage: huge-merge.pl [OPTIONS] SOURCEDIR\n";
    askHelp();
}

# function to output "ask for help" message when the user's goofed up!
sub askHelp
{
    print STDERR "Type huge-merge.pl --help for help.\n";
}

# function to output help messages for this program
sub showHelp
{
    print "\n";
    print "Usage: huge-merge.pl [OPTIONS] SOURCEDIR\n\n";
    print "huge-merge.pl takes the directory of sorted bigram files generated \n";
    print "by huge-sort.pl as input, and output one sorted bigram file. \n";
    print "The output file is in the source folder. The individual bigram files\n"; 
    print "that have been merged are deleted by default.\n\n";

    print "OPTIONS:\n\n";

    print "  --keep             Keep the unmerged files.\n";
    print "                     The unmerged files are deleted by default.\n\n";

    print "  --help             Prints this help message.\n";
    print "  --version          Prints this version message.\n";
}

# function to output the version number
sub showVersion
{
    print STDERR 'huge-merge.pl $Id: huge-merge.pl,v 1.26 2011/03/31 23:04:04 tpederse Exp $';
    print STDERR "\nCopyright (C) 2009-2011, Ying Liu\n";

}

