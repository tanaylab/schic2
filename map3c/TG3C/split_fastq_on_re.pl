use strict;

sub gen_split {
	my($l1, $l2, $l3, $l4, $pr, $outf, $re) = @_;

	my($l_re) = length($re);
	my($i) = -1;
	my($base) = 0;
	my($pos) = index($l2, $re, $base);
	my($seg) = 0;

	my(@head) = split(" ", $l1);
	my($head_suf) = join(" ", @head[1..$#head]);

	while($pos > $base) {
		if($pos - $base + $l_re >= $::min_len) {
			print $outf "$head[0]:$pr"."SEG$seg:$base:".($pos + $l_re)." $head_suf\n";
			print $outf substr($l2, $base, $l_re + $pos - $base)."\n";
			print $outf "$l3\n";
			print $outf  substr($l4, $base, $l_re + $pos - $base)."\n";
			$::count->[$seg]++;
			$seg++; 
		}
		$base = $pos;	
		$pos = index($l2, $re, $base+1);
	}
	$pos = length($l2);
	if($pos - $base >= $::min_len) {
		print $outf "$head[0]:$pr"."SEG$seg:$base:$pos $head_suf\n";
		print $outf substr($l2, $base, $l_re + $pos - $base)."\n";
		print $outf "$l3\n";
		print $outf substr($l4, $base, $l_re + $pos - $base)."\n";
		$::count->[$seg]++;
		$seg++; 
	}
}
if($#ARGV < 3) {
	die "usage  split_fastq_re.pl in_dir in_fn_regexp out_fastq re_seq min_len pair_mode [first_end_code second_end_code]\n";
}

my(@pre_fns) = <$ARGV[0]/*.fastq>;
my($fn_regexp) = $ARGV[1];

print STDERR "looking at $ARGV[0] regexp is $fn_regexp\n";
my(@fns);
my($fn);
foreach $fn (@pre_fns) {
	if($fn =~/$fn_regexp/) {
		push(@fns, $fn);
	} 
}

my($out_fq) = $ARGV[2];
my($re_seq)= $ARGV[3];
$::min_len = $ARGV[4];

my($pair_mode) = 1;
if(defined($ARGV[5])) {
	$pair_mode = $ARGV[5];
}
my($first_end) = "R1";
my($second_end) = "R2";

if(defined($ARGV[6])) {
	$first_end = $ARGV[6];
	$second_end = $ARGV[7];
}


my($fn);

my($l1,$l2,$l3,$l4) = @_;

my($paired);

open(OUT,  ">$out_fq") || die "cannot write segmented fastq output at $out_fq\n";
foreach $fn (@fns) {
	my($n) = 0;
	$::stat_count = [];
	if($pair_mode && $fn !~/$first_end/) {
		next;
	}
	open(FASTQ, $fn) || die "Cannot open fastq at $fn\n";

	my($fn_pair) = $fn;
	$fn_pair =~ s/$first_end/$second_end/;
	
	$paired = 0;
	if(-e $fn_pair) {
		$paired = 1;
		print STDERR "also paired $fn_pair\n";
		open(FASTQ2, $fn_pair) || die "Cannot open fastq at $fn_pair\n";
	} 
	print STDERR "Segment processing $fn";

	while(<FASTQ>) {
		$l1 = $_; chop $l1;
		$l2 = <FASTQ>; chop $l2;
		$l3 = <FASTQ>; chop $l3;
		$l4 = <FASTQ>; chop $l4;

		gen_split($l1,$l2,$l3,$l4, "R1", \*OUT, $re_seq);
		if($paired) {
			$l1 = <FASTQ2>; chop $l1;
			$l2 = <FASTQ2>; chop $l2;
			$l3 = <FASTQ2>; chop $l3;
			$l4 = <FASTQ2>; chop $l4;
			gen_split($l1,$l2,$l3,$l4, "R2", \*OUT,$re_seq);
		}
		
		if($n % 100000 == 0) {
			print STDERR "cnt $n 0: $::count->[0], 1: $::count->[1], 2: $::count->[2], 3: $::count->[3]\n";
		}
		$n++;
	}
}
close OUT;

print STDERR "total:\t0: $::count->[0]\t1: $::count->[1]\t2: $::count->[2]\t3: $::count->[3]\n";
