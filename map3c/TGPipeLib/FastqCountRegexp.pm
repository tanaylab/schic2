use strict;

package TGPipeLib::FastqSplicer;

1;


sub required() {
	return("re_code");
}

sub proc_line
{
	my($split_id, $params, $fastq, $fastq2) = @_;

	my($n) = 0;
	$fastq->first();
	my($re1) = $params->{regexp1};

	my(%count);
	if($fastq2 eq "NA") {
		while($fastq->valid()) {
			$fastq->next();
			if($re1 eq "" || $fastq->cur_seq()=~/$re1/) {
				my($res1) = ($re1 eq "") ? "" : $fastq->cur_seq()=~/$re1/;
				$count{$res1}++;
			} 
			$n++;
		}
	} else {
		$fastq2->first();
		my($re2) = $params->{regexp2};
		while($fastq->valid()) {
			$fastq->next();
			$fastq2->next();
			if(($re1 eq "" || $fastq->cur_seq()=~/$re1/)
			&& ($re2 eq "" || $fastq2->cur_seq()=~/$re2/)) {
				my($res1) = ($re1 eq "") ? "" : $fastq->cur_seq()=~/$re1/;
				my($res2) = ($re2 eq "") ? "" : $fastq2->cur_seq()=~/$re2/;
				$count{$res1.$res2}++;
			} 
			$n++;
		}
	}
	my($k);
	my($out_fn) = $params->{work_dir}."/stat.$split_id";
	if(open(RPT, ">$out_fn") ) {
		foreach $k (keys %count) {
			print RPT "$k\t$count{$k}\n";
		}
		close RPT;
		return($out_fn);
	} else {
		return("ERR:cannot write $split_id");
	}
}

sub reduce($) {
	my($params) = @_;

	my($wd) = $params->{work_dir};
	my(@out_fns) = <$wd/stat.*>;

	print STDERR "merge $#out_fns files\n";
	my(%all);
	my($i) = 0;
	my($fn);
	foreach $fn (@out_fns) {
		if(!open(CNT, $fn) ) {
			print STDERR "cannot open count report for split $i. $fn\n";
			next;
		}
		print STDERR "reduxing $i fn $fn\n";
		while(<CNT>) {
			chop;
			my($k, $v) = split("\t", $_);
			$all{$k} += $v;
		}
		$i++;
	}
	my($k);
	open(OUT, ">".($params->{out_fn})) || die "cannot write output!\n";
	foreach $k (keys %all) {
		print OUT "$k\t$all{$k}\n";
	}
}
