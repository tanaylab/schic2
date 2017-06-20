#!/usr/bin/env perl

use strict;
use warnings FATAL => qw( all );
use File::Basename;
use List::Util qw(max min);

if ($#ARGV == -1) {
	print "usage: $0 <fa file> <output dir> <step> <read length> <maximal reads per file>\n";
	exit 1;
}
my $fa = $ARGV[0];
my $odir = $ARGV[1];
my $step = $ARGV[2];
my $read_length = $ARGV[3];
my $max_reads = $ARGV[4];

open(IN, $fa) or die $fa;
my $chr;
my $coord = 0;
my %seqs;
while (my $line = <IN>)
{
	chomp($line);
	if ($line =~ /^>chr(.*)/) {
		$chr = $1;
		print STDERR "reading sequence of chr$chr\n";
		next;
	}
	$seqs{$chr} = "" if (!defined($seqs{$chr}));
	$seqs{$chr} .= $line;
}
close(IN);

my $index = 1;
my $ofile = $odir."/".$index.".fastq";

system("mkdir -p $odir");
open(OUT, ">", $ofile) || die $ofile;
print STDERR "generating file $ofile\n";

my $count = 1;
foreach my $chr (keys %seqs)
{
	my $seq = $seqs{$chr};
	my $end = length($seq);
	my $coord = $step;
	while ($coord < $end)
	{
		if ($coord + $read_length < $end)
		{	
			my $read = substr($seq, $coord-1, $read_length);
			#print "$chr\t$coord\t$read\n";

			print OUT "\@MAPABILITY_".$chr."_".$coord."\n";
			print OUT $read, "\n";
			print OUT "+\n";
			print OUT "<" x $read_length, "\n";
		}

		if ($max_reads > 0 and $count == $max_reads)
		{
			close(OUT);
			$index++;
			$count = 0;
			my $ofile = $odir."/".$index.".fastq";
			open(OUT, ">", $ofile) || die;
			print STDERR "generating file $ofile\n";
		}

		$count++;
		$coord += $step;
	}
}
close(OUT);
