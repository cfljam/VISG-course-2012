#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

huge-split.pl - Split bigram files from huge-count.pl into pieces.

=head1 DESCRIPTION

See perldoc huge-split.pl 

=head1 USAGE

huge-split.pl [OPTIONS] SOURCE

=head1 INPUT

=head2 Required Arguments:

=head3 SOURCE

Input to huge-split.pl should be a file generated by huge-count.pl or 
count.pl with tokenlist option. The results files have the same name 
with the input source file and each split file has an extention 
sequence number.  

=head4 --split N 

This parameter should be set. huge-split will divide the output bigrmas
tokenlist generated by count.pl or huge-count.pl. Each part created with 
--split N will contain N lines. Value of N should be chosen such that 
huge-sort.pl can be efficiently run on any part containing N lines from 
the file contains all bigrams file. 

We suggest that N is equal to the number of KB of memory you have. If the 
computer has 8 GB RAM, which is 8,000,000 KB, N should be set to 8000000. 

=head3 Other Options :

=head4 --help

Displays this message.

=head4 --version

Displays the version information.

=head1 AUTHOR

Amruta Purandare, Ted Pedersen, Ying Liu.
University of Minnesota at Duluth.

=head1 COPYRIGHT

Copyright (c) 2004-2011

Ted Pedersen, University of Minnesota, Duluth.
tpederse@umn.edu

Ying Liu, University of Minnesota, Twin Cities.
liux0395@umn.edu

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to

The Free Software Foundation, Inc.,
59 Temple Place - Suite 330,
Boston, MA  02111-1307, USA.

=cut

###############################################################################

#                           ================================
#                            COMMAND LINE OPTIONS AND USAGE
#                           ================================

# command line options
use Getopt::Long;
GetOptions ("help","version","split=i");
# show help option
if(defined $opt_help)
{
        $opt_help=1;
        &showhelp();
        exit;
}

# show version information
if(defined $opt_version)
{
        $opt_version=1;
        &showversion();
        exit;
}


#############################################################################

#			========================
#			      CODE SECTION
#			========================

use Getopt::Long;

# first check if no commandline options have been provided... in which case
# print out the usage notes!
if ( $#ARGV == -1 )
{
    &showminimal();
    exit;
}

# now get the options!
GetOptions( "keep", "version", "help" );


my $file = $ARGV[0];
open SPLIT, "<$file" or die $!;

if (!defined $opt_split)
{
	print STDERR "Warning($0): You did not specify the split size. huge-split.pl\n";
	print STDERR "will not split the whole bigrams file into smaller pieces.\n";
	exit;
}

my $split_num = 0;
my $sub_i = 1;
my $sub_file = "$file" . ".$sub_i";
open(SUB, ">$sub_file") or die("Error: cannot open file '$sub_file' for output index.\n");
while (my $line = <SPLIT>)
{
	print SUB "$line";
   	$split_num++;

   	if ($split_num == $opt_split)
   	{
   		close SUB;
   		if (eof (SPLIT))
   		{
   			last; 
   		}
   		else
   		{
   			$sub_i++;
   			$split_num = 0;
   			my $sub_file = "$file" . ".$sub_i";
   			open(SUB, ">$sub_file") or die("Error: cannot open file '$sub_file' for output index.\n");
   		}
   }
}

close SPLIT;
#system("/bin/rm $file");


##############################################################################

#                      ==========================
#                          SUBROUTINE SECTION
#                      ==========================

#-----------------------------------------------------------------------------
#show minimal usage message
sub showminimal()
{
        print "Usage: huge-split.pl [OPTIONS] SOURCE";
        print "\nTYPE huge-split.pl --help for help\n";
}

#-----------------------------------------------------------------------------
#show help
sub showhelp()
{
	print "Usage:  huge-split.pl [OPTIONS] SOURCE\n";

	print "OPTIONS:\n\n";

    print "  --split N          Split the bigram file into smaller files. Every \n";
    print "                     smaller file contains N bigrams. N must be an integer. \n\n";

    print "  --help             Prints this help message.\n";
    print "  --version          Prints this version message.\n";
}

#------------------------------------------------------------------------------
#version information
sub showversion()
{
        print 'huge-split.pl $Id: huge-split.pl,v 1.12 2011/03/31 23:04:04 tpederse Exp $';
        print "\nCopyright (C) 2004-2011, Amruta Pruandare, Ted Pedersen & Ying Liu.\n";
}



#############################################################################

