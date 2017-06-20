use strict;

my(@fns) = split(",", $ARGV[0]);
my($out_fn) = $ARGV[1];
my($mode) = "full";
if(defined($ARGV[2])) {
	$mode = $ARGV[2];
}

open(OUT, ">$out_fn") || die "cannot write merged adj table at $out_fn\n";

if($#fns == 0 && $ARGV[0]!~/adj/) {
	@fns = <$ARGV[0]*adj>;
}

my($fn);
my($ni) = 0;
my(%all, %adj, %bin, %ndx_bin);
my(@ndxs);
foreach $fn (@fns) {
	my($ndx) = $fn=~/\/([^\/]+)[\/*]adj$/;
	push(@ndxs, $ndx);

	print STDERR "Processing $fn...\n";

	open(TAB, $fn) || die "cannot open file $fn\n";
	while(<TAB>) {
		chop;
		my(@f) = split("\t", $_);
		
		if ($f[0] !~ /^\d+/) {
			next;
		}

		if ($mode eq "raw_sum") {
			my $p = $f[0] > $f[1] ? "$f[1]\t$f[0]" : "$f[0]\t$f[1]";
			$bin{$p} += $f[2];
		}
		elsif ($mode eq "binary_per_adj") {
			my $p = $f[0] > $f[1] ? "$f[1]\t$f[0]" : "$f[0]\t$f[1]";
			if (exists($ndx_bin{$p}) and $ndx_bin{$p} == $ni) {
				next;
			}
			$ndx_bin{$p} = $ni;
			$bin{$p} += 1; 
		}
		else {
			$adj{"$f[0]\t$f[1]\t$ni"} = $f[2];
		}
		$all{"$f[0]\t$f[1]"} = 1;
	}
	$ni++;
}
my($a);
my($naked) = 0;
if($mode eq "to_imp" or $mode eq "binary_per_adj" or $mode eq "raw_sum") {
	$naked = 1;
	print OUT "fend1\tfend2\tcount\n";
} else {
	print OUT "fid1\tfid2\t".join("\t", @ndxs)."\tsum\tlsum\n";
}
foreach $a (keys %all) {
	my($sum) = 0;
	my($suml) = 0;
	if($all{$a} == 0) {
		next;
	}
	my($fr,$to) = split("\t", $a);
	my($rev_a) = "$to\t$fr";
	if(exists($all{$rev_a})) {
		$all{$rev_a} = 0;
	}
	if ($mode eq "binary_per_adj" or $mode eq "raw_sum") {
		my $p = $fr > $to ? $rev_a : $a;
		print OUT "$p\t" . $bin{$p} . "\n";
	} else {
		print OUT "$a";
		for(my($i) = 0; $i < $ni; $i++) {
			my($v) = (int($adj{"$rev_a\t$i"}+int($adj{"$a\t$i"})));
			if(!$naked) {
				print OUT "\t$v";
			}
			$sum += $v;
			$suml += log(1+$v);
		}
		if($naked) {
			print OUT "\t".(int(100*($ni*(exp($suml/$ni)-1)))/100)."\n";
		} else {
			print OUT "\t$sum\t".(int(100*($ni*(exp($suml/$ni)-1)))/100)."\n";
		}
	}
}
