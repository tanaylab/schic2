use strict;
use Getopt::Long;
sub pcnt {
	my @vec = @{$_[0]};
	
	my $min_base_qual = $_[1];
	my $i = 0;
	foreach (@vec) {
		next if $_ !~ /\d+/;
		$i++ if $_ >= $min_base_qual;
	}
	#print "i = $i; len = " . scalar(@vec) . "\n";
	return ($i / scalar(@vec));

}

sub pcnt_from_array {
	print STDERR "will extract percentile\n";
	my @matrix = @{$_[0]};
	my $min_base_qual = $_[1];
	my @out;
	
	foreach (@matrix) {
		my @vec = split(/ /, $_);
		
		#my $sum = map { $sum += $_ } @vec;
		#push(@avg_quals1, $sum / scalar(@vec));
		my $pcnt_gt_min = pcnt(\@vec, $min_base_qual);
		push (@out, $pcnt_gt_min);
		#print STDERR "$cnt reads calc...\t";
	}
	return(@out);
}

sub avg_read_qual {
	my @qual = @_;
	my $sum;
	foreach (@qual) {
		$sum += $_;
	}
	
	return($sum / scalar(@qual));
}

sub gen_split {
	my($l1, $l2, $l3, $l4, $pr, $outf, $re, $min_frag_len) = @_;
	
	my($l_re) = length($re);
	my($i) = -1;
	my($base) = 0;
	my($pos) = index($l2, $re, $base);
	my($seg) = 0;

	my(@head) = split(" ", $l1);
	my($head_suf) = join(" ", @head[1..$#head]);

	while($pos > $base) {
		if($pos - $base >= $min_frag_len) {
			print $outf "$head[0]:$pr"."SEG$seg:$base:$pos $head_suf\n";
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
	if($pos - $base > $min_frag_len) {
		print $outf "$head[0]:$pr"."SEG$seg:$base:$pos $head_suf\n";
		print $outf substr($l2, $base, $l_re + $pos - $base)."\n";
		print $outf "$l3\n";
		print $outf substr($l4, $base, $l_re + $pos - $base)."\n";
		$::count->[$seg]++;
		$seg++; 
	}
}

################code####################

die "usage: perl trim_fastqs_by_qual.pl --scope_dir dir --out_fn out_fn --regexp fn_regexp
--min_base_qual 30 --reads_to_sample 1000 --min_len 30 --min_read_qual 20
--re_seq GATC --min_frag_len 20\n" if scalar(@ARGV) < 1;
my ($scope_dir, $out_fn, $min_base_qual, $fn_regexp, $min_read_qual, $re_seq, $min_frag_len);
my $sample_n = 0;
my $min_len = 30;
my $read1_pad = 30;
my $pos_err_thres = 0.8;

GetOptions("scope_dir=s" => \$scope_dir,
			"regexp=s" => \$fn_regexp,
			"out_fn=s" => \$out_fn,
			"min_base_qual=i" => \$min_base_qual,
			"reads_to_sample=i" => \$sample_n,
			"min_len=i" => \$min_len,
			"min_read_qual:i" => \$min_read_qual,
			"re_seq=s" => \$re_seq,
			"min_frag_len=i" => \$min_frag_len)
	or die "not all args supplied\n";
$min_read_qual = 20 if $min_read_qual < 1;
	
	
my @fastqs = <$scope_dir/*fastq>;

print STDERR "looking at $scope_dir regexp is $fn_regexp\n";
my(@fns);
my($fn);
foreach $fn (@fastqs) {
	if($fn =~/$fn_regexp/) {
		push(@fns, $fn);
	} 
}

my($first_end) = "R1";
my($second_end) = "R2";

open(my $out_fh, ">", $out_fn) or die "cannot write to $out_fn";
foreach my $fn (@fns) {
	if ($fn !~ /$first_end/) {
		next;
	}
	my ($in_fh, $in_fh2);
	open ($in_fh, "<", $fn) or die "cannot open $fn\n";
	my $paired;
	my $fn2 = $fn;
	$fn2 =~ s/$first_end/$second_end/;
	if(-e $fn2) {
		$paired = 1;
		print STDERR "also paired $fn2\n";
		open ($in_fh2, "<", $fn2) or die "cannot open $fn2\n";
	} 
	
	print STDERR "proccessing $fn and $fn2\n";
	
	my $cnt = 0;
	my $cnt_q = 0; #reads that were used
	my @all_quals1;
	my @all_quals2;
	while (<$in_fh>) {
		
		if ($cnt % 123 != 0) {
			$cnt++;
			next;
		}	
	
		my $l11 = $_; chomp $l11;
		my $l12 = <$in_fh>; chomp $l12;
		my $l13 = <$in_fh>; chomp $l13;
		my $l14 = <$in_fh>; chomp $l14;
		my @qual1 = map {(ord) - 33 } split(//, $l14);
		
		#creates a matrix with nrows = length of read - each entry represents base pos; ncols = num of reads
		
		for(my $k = 0; $k <= $#qual1; $k++) {
			$all_quals1[$k] .= $qual1[$k] . " ";
		}
		
		if ($paired) {
			my $l21 = <$in_fh2>; chomp $l21;
			my $l22 = <$in_fh2>; chomp $l22;
			my $l23 = <$in_fh2>; chomp $l23;
			my $l24 = <$in_fh2>; chomp $l24;
			my @qual2 = map {(ord) - 33 } split(//, $l24);
		
			for(my $k = 0; $k <= $#qual2; $k++) {
				$all_quals2[$k] .= $qual2[$k] . " ";
			}
		}
		
		
		if ($cnt % 10000 == 0) {
			print STDERR "$cnt_q reads...\t";
		}
		
		$cnt++;
		$cnt_q++;
		if ($sample_n > 0 && $cnt_q >= $sample_n) {
			last;
		}
	}
	close($in_fh);
	close($in_fh2);

	################Calc base stats################
	#arrays with percentile of bases passing $min_base_qual in each pos
	my @pcnt1 = pcnt_from_array(\@all_quals1, $min_base_qual);
	my @pcnt2 = pcnt_from_array(\@all_quals2, $min_base_qual);
	my ($trim_base1, $trim_base2);
	my $trimmed1 = 0;
	my $trimmed2 = 0;
	#decide base to trim from
	for ($trim_base1 = -1 + $read1_pad + $min_len; $trim_base1 <= $#pcnt1; $trim_base1++ ) {
		if ($pcnt1[$trim_base1] < $pos_err_thres ) {
			print STDERR "will trim read1 from pos ";
			print STDERR 1+$trim_base1 . "; % < qual = $pcnt1[$trim_base1]\n";
			$trimmed1 = 1;
			last;
		}
	}
	print STDERR "no trimming done to read1\n" if $trimmed1 == 0;
	
	for ($trim_base2 = -1 + $min_len; $trim_base2 <= $#pcnt2; $trim_base2++ ) {
		if ($pcnt2[$trim_base2] < $pos_err_thres ) {
			print STDERR "will trim read2 from pos ";
			print STDERR 1+$trim_base2 . "; % < qual = $pcnt2[$trim_base2]\n";
			$trimmed2 = 1;
			last;
		}
	}
	print STDERR "no trimming done to read2\n" if $trimmed2 == 0;
	
	################Trim and split fastqs##########
	open ($in_fh, "<", $fn) or die "cannot open $fn\n";
	open ($in_fh2, "<", $fn2) or die "cannot open $fn2\n";
	my ($qual_filt, $read_cnt);
	
	while (<$in_fh>) {
		$read_cnt++;
		my (@qual1, @qual2);
		my ($avg_qual1, $avg_qual2);
		my $l11 = $_; chomp $l11;
		my $l12 = <$in_fh>; chomp $l12;
		my $l13 = <$in_fh>; chomp $l13;
		my $l14 = <$in_fh>; chomp $l14;
		
		if ($trimmed1) {
			$l12 = substr($l12, 0, $trim_base1 + 1);
			$l14 = substr($l14, 0, $trim_base1 + 1);
		}
		
		@qual1 = map {(ord) - 33 } split(//, $l14);   #filters out low qual reads - is there a faster way to do it???
		#for(my $i = 0; $i <= length($l14); $i++) {
		#	my $n = substr($l14, $i, 1);
		#	push @qual1, (ord($n) - 33);
		#}
		$avg_qual1 = avg_read_qual(@qual1);
		if ($avg_qual1 < $min_read_qual && !$paired) {
			$qual_filt++;
			next;
		}
		
		if ($paired) {
			my $l21 = <$in_fh2>; chomp $l21;
			my $l22 = <$in_fh2>; chomp $l22;
			my $l23 = <$in_fh2>; chomp $l23;
			my $l24 = <$in_fh2>; chomp $l24;
			
			if ($trimmed2) {
				$l22 = substr($l22, 0, $trim_base2 + 1);
				$l24 = substr($l24, 0, $trim_base2 + 1);
			}
			
			
			@qual2 = map {(ord) - 33 } split(//, $l24);
			#for(my $i = 0; $i <= length($l24); $i++) {
                	#        my $n = substr($l24, $i, 1);
        	        #        push @qual2, (ord($n) - 33);
	                #}

			$avg_qual2 = avg_read_qual(@qual2);
			
			if ($avg_qual1 < $min_read_qual || $avg_qual2 < $min_read_qual) {
				$qual_filt++;
				next;
			}
			
			gen_split($l11,$l12,$l13,$l14, "R1", $out_fh, $re_seq, $min_frag_len);
			gen_split($l21,$l22,$l23,$l24, "R2", $out_fh, $re_seq, $min_frag_len);
		} else {
			gen_split($l11,$l12,$l13,$l14, "R1", $out_fh, $re_seq, $min_frag_len);
		}
		
		if ($read_cnt % 100000 == 0) {
			printf STDERR "proccessed %uK reads... ", $read_cnt/1000;	
		}
		
	}
	
		
	print STDERR "\n$fn: filtered $qual_filt out of $read_cnt reads based on low qual\n";
	
	
}
