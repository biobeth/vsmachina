#!/usr/bin/perl
#a script to stop me screwing up esl-sfetch 
#takes gff from STDIN, and pre-built esl-sfetch index
#-g names sequences "seqname_feature" for multi-genome, multi-gene gff files 
use strict; use warnings;
use Getopt::Long;

my ($help,$name);
&GetOptions(
	"h|help"	=>	\$help,
	"g|genome"	=>	\$name
	);

if($help){
	&help();
	exit(1);
}


if (-t STDIN){
	print "Error: No input from STDIN\n";
	&help();
        exit(1);
}


my $index=$ARGV[0];

while(<STDIN>){
	my $line = $_;
	my @columns = split('\t',$line);
	if($name){
		$columns[2] = $columns[0]."\_".$columns[2];
	}
	if($columns[6] eq "+"){
		system("esl-sfetch -n $columns[2] -c $columns[3]\.\.$columns[4] $index $columns[0]");
		}
	 elsif($columns[6] eq "-"){
		system("esl-sfetch -r -n $columns[2] -c $columns[3]\.\.$columns[4] $index $columns[0]");
		}
}

sub help{
	print "fetchGff.pl\nUsage [Gff from STDIN} | fetchGff.pl index_file\n\nRequires pre-built esl-sfetch index.\n\nOptions:\n-------\n-h --help\tDisplay this help\n\n";
}
