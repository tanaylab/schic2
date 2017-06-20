use strict;

use IO::File;
use File::Basename;
use GDBM_File;

require map3c::ReDB;


#simplified version that take mapped pairs reads and create adj table
#
#In the full mode it report and adj for each mapped pair (fragment id, strand, offset)
#In the naive mode it generate an adj table for trackdb import : fend1,fend2,count, with duplicates elimination.
#
if($#ARGV < 4) {
	die "sampair2adj redb_dir re_seq sam_index_file output naive no_pcr_dup [min_q]\n";
}

my ($no_dup) = $ARGV[5] eq "T";

my($min_q) = 38;
if(defined($ARGV[6])) {
	$min_q = $ARGV[6];
}
open(OUT, ">$ARGV[3]") || die "cannot open output adj tab $ARGV[3]\n";
my($naive_track) = $ARGV[4] eq "naive";
if (!$naive_track) {
	print OUT "fid1\tfid2\tstrand1\tstrand2\toff1\toff2\n";
}

my($odir) = dirname($ARGV[3]);

my(%adj_count); #to store adjs in naive mode
#my(%unique_adj_ind); # storing yes/no if a specific fend-offset pair appeared (to eliminate PCR amplification)

#my $adj_count_fn = "/tmp/".($$)."_adj_count";
#my $unique_adj_ind_fn = "/tmp/".($$)."_unique_adj_ind";

#tie %adj_count, 'GDBM_File', $adj_count_fn, &GDBM_WRCREAT, 0640;
#tie %unique_adj_ind, 'GDBM_File', $unique_adj_ind_fn, &GDBM_WRCREAT, 0640;

my($redb) = ReDB::new("ReDB");
print STDERR "will init re db\n";
$redb->init_from_tabs($ARGV[0]."/$ARGV[1].frags");
print STDERR "done read tab\n";

my($search_scope) = 5000;

my($lowq_count) = 0;
my($trans_count) = 0;
my($self_count) = 0;
my($nomap_count) = 0;
my($nore_count) = 0;
my($dup_count) = 0;
my(@count_d);
my($cnt) = 0;
my($m1,$m2);
my($fend1, $fend2);

my($sam) = $ARGV[2];
open(SAM, $sam) || die("could not open sam index file");

my $skip = 0;
if (!$skip) {
while (<SAM>) {
	chomp;
	my($s1, $s2) = split("\t", $_);
	print STDERR "sam files: $s1 , $s2\n";
	my($sam1) = IO::File::new();
	$sam1->open($s1) || die "cannot open sam1 $s1\n";
	my($sam2);
	my($pair_embed) = 0;
	if($s2 eq ".") {
		$pair_embed = 1;
		$sam2 = $sam1;
	} else {
		$sam2 = IO::File::new();
		$sam2->open($s2) || die "cannot open sam2 $s2\n";
	}
	while(!$sam1->eof()) {
		$m1 = <$sam1>;
		if($sam1->eof() || $m1 =~/^\@/) {
			if(!$pair_embed) {
				$m2 = <$sam2>;
			}
			next;
		}

		$m2 = <$sam2>;

		my(@f1) = split("\t", $m1);
		my(@f2) = split("\t", $m2);
		if($f1[4] < $min_q || $f2[4] < $min_q) {
			$lowq_count++;
			next;
		}
		my($chr1) = $f1[2];
		my($chr2) = $f2[2];
		my($pos1) = $f1[3];
		my($pos2) = $f2[3];
		my($strand1) = $f1[1] & 16 ? -1 : 1;
		my($strand2) = $f2[1] & 16 ? -1 : 1;
		my($l1) = length($f1[9]);
		my($l2) = length($f2[9]);
		my($epos1) = $pos1 + $l1;
		my($epos2) = $pos2 + $l2;
		if($f1[5] =~/^\d+S/) {
			my($gap) = $f1[5]=~/^(\d+)S/;
			$pos1 += $gap;
			($gap) = $f1[5]=~/(\d+)S$/;
			$epos1 -= $gap;
		}	
		if($f2[5] =~/^\d+S/) {
			my($gap) = $f2[5]=~/^(\d+)S/;
			$pos2 += $gap;
			($gap) = $f2[5]=~/(\d+)S$/;
			$epos2 -= $gap;
		}

		my($fid1) = $redb->map_frag($chr1, $strand1 == 1 ? $pos1 + 20 : $epos1 - 20, $search_scope);
		my($fid2) = $redb->map_frag($chr2, $strand2 == 1 ? $pos2 + 20 : $epos2 - 20, $search_scope);

		$cnt++;

		if($fid1 eq "-1" || $fid2 eq "-1") {
			$nomap_count++;
			next;
		}

		if($fid1 eq $fid2) {
			$self_count++;
			next;
		}

		my($off1) = $strand1 eq 1 ? 
			$redb->frag_3_dist($fid1, $pos1) :
			$redb->frag_5_dist($fid1, $epos1);
		my($off2) = $strand2 eq 1 ? 
			$redb->frag_3_dist($fid2, $pos2) :
			$redb->frag_5_dist($fid2, $epos2);

		$fend1 = $fid1 * 2 + ($strand1 == 1 ? 1 :0);
		$fend2 = $fid2 * 2 + ($strand2 == 1 ? 1 :0);

		if ($fend2 < $fend1) {
			my($temp) = $fend1;
			$fend1 = $fend2;
			$fend2 = $temp;
			$temp = $off1;
			$off1 = $off2;
			$off2 = $temp;
		}
		
		if($naive_track) {
			if ($no_dup) {
				print OUT "$fend1\t$fend2\t$off1\t$off2\n";
			} else {
				$adj_count{"$fend1\t$fend2"}++;
			}
		} else {
			if ($fid1 < $fid2) {
				print OUT "$fid1\t$fid2\t$strand1\t$strand2\t$off1\t$off2\n";
			}
			else {
				print OUT "$fid2\t$fid1\t$strand2\t$strand1\t$off2\t$off1\n";
			}
		}
		if($chr1 ne $chr2) {
			$trans_count++;
		} elsif (abs($fend2 - $fend1) == 1) {
			$nore_count++;
		}
		else {
			my($d) = int(log(1+abs($pos1 -$pos2))/log(10));
			if($d < 2) {
				$d = 2;
			}
			$count_d[$d]++;
		}

		if($cnt % 5000000 == 0) {
			my($tot) = 1 + $trans_count + $self_count + $nore_count;
			for(my($d) = 2; $d < 9; $d++) {
				$tot += $count_d[$d];
			}
			printf STDERR ("CNT %dK", int($cnt)/1000);
			printf STDERR ("\tMAP %dK", int($lowq_count/1000));
			printf STDERR ("\tDUP %d", $dup_count);
		 
			printf STDERR ("\tRND %.3f", $nomap_count/$tot);
			printf STDERR ("\tTR  %.3f", $trans_count/$tot);
			printf STDERR ("\tSLF %.3f", $self_count/$tot);
			printf STDERR ("\tNoR %.3f", $nore_count/$tot);
			
			for(my($d) = 2; $d < 9; $d++) {
				printf STDERR ("\t%.3f", $count_d[$d]/$tot);
			}
			print STDERR "\n";
		}
	}
}}
print STDERR "Done reading...\n";

if($naive_track) {
	if ($no_dup) {
		close(OUT);
		system("mv $ARGV[3] $ARGV[3].tmp");
		open(OUT, ">$ARGV[3]") or die "failed to open $ARGV[3] for output";
		print OUT "fend1\tfend2\tcount\n";
		close(OUT);
		system("sort -k 1n,4n $ARGV[3].tmp | uniq | cut -f 1,2 | uniq -c | awk '{ print \$2, \$3, \$1;}' | tr ' ' '\t'  >> $ARGV[3]");
	}
	else {
		print OUT "fend1\tfend2\tcount\n";
		my($k);
		foreach $k (keys %adj_count) {
			print OUT "$k\t$adj_count{$k}\n";
		}
	}
}

my($tot) = 1 + $trans_count + $self_count + $nore_count;

for(my($d) = 2; $d < 9; $d++) {
	$tot += $count_d[$d];
}
printf STDERR ("CNT %dK", int($cnt)/1000);
printf STDERR ("\tMAP %dK", int($lowq_count/1000));
printf STDERR ("\tDUP %d", $dup_count);

printf STDERR ("\tRND %.3f", $nomap_count/$tot);
printf STDERR ("\tTR  %.3f", $trans_count/$tot);
printf STDERR ("\tSLF %.3f", $self_count/$tot);
printf STDERR ("\tNoR %.3f", $nore_count/$tot);

for(my($d) = 2; $d < 9; $d++) {
	printf STDERR ("\t%.3f", $count_d[$d]/$tot);
}
print STDERR "\n";
