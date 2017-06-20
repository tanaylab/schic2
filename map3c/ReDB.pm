package ReDB;

use strict;

#use diagnostics;

sub new($) {
	my($clas) = @_;

	my($self) = {};
	
	bless($self,$clas);

	$self->{coord} = [];
	$self->{chrom} = [];
	$self->{frag_len} = [];
	$self->{binsize} = 10;
	$self->{coor_to_frag} = {};

	return($self);
}

sub init_from_tabs {
	my($self, $frag_tab) = @_;

	open(FRAG, $frag_tab) || die "cannot open frag tab $frag_tab\n";

	my($ctf) = $self->{coor_to_frag};
	my($flen) = $self->{frag_len};
	my($fchrom) = $self->{chrom};
	my($fcoord) = $self->{coord};

	my($bs) = $self->{binsize};

	print STDERR "loading the re db\n";
	while(<FRAG>) {
		chop;
		my($fid, $chrom, $start, $end) = split("\t", $_);

		$chrom=~s/chr//;

		if($end - $start < $bs * 1.5) {
			next;
		}

#		my($mid_bin) = int(($start+$bs)/$bs);	
#		$ctf->{"$chrom\t$mid_bin"} = "$fid\tM";

		my($start_bin) = int(($start+$bs)/$bs);	
		$ctf->{"$chrom\t$start_bin"} = "$fid\t5";

		my($end_bin) = int(($end-$bs)/$bs);	
		$ctf->{"$chrom\t$end_bin"} = "$fid\t3";

		$flen->[$fid] = $end-$start;
		$fchrom->[$fid] = $chrom;
		$fcoord->[$fid] = $start;
		if($fid % 100000 == 0) {
			print STDERR "..$fid";
		}
		#if($fid >1000000) {
		#	last;
		#}
	}
}

sub build_from_genome {
	my($self, $chrom_key, $re_seq, $f_offset, $r_offset, $max_frag_len, $out_dir) = @_;

	if (!-e "$out_dir") {
		mkdir("$out_dir");
	}
	
	open(FRAGS, ">$out_dir/$re_seq.frags") || die "cannot write frags\n";
	open(FENDS, ">$out_dir/$re_seq.fends") || die "cannot write fends\n";
	open(FRAGS_CONT, ">$out_dir/$re_seq.frags_cont") || die "cannot write frags\n";
	open(FRAGS_GC, ">$out_dir/$re_seq"."_gc") || die "cannot write frags\n";
	open(FRAGS_FLEN, ">$out_dir/$re_seq"."_flen") || die "cannot write frags\n";

	print FENDS "fend\tchr\tcoord\n";
	print FRAGS_GC "chrom\tstart\tend\tgc\n";
	print FRAGS_FLEN "chrom\tstart\tend\tflen\n";

	open(KEY, $chrom_key) || die "cannot open chrom key\n";

	my($fid) = 0;
	while(<KEY>) {
		chop;
		my($chrom, $fn) = split("\t", $_);
		print STDERR "will proc $chrom\n";
		open(SEQ, $fn) || die "cannot opem crhom $chrom seq at $fn\n";
		my($seq) = "";
		<SEQ>;
		while(<SEQ>) {
			chop;
			$seq .= uc($_);
		}
		print STDERR "done read ".length($seq)." bp";
		my($prev_pos) = 0;
		my($re_pos) = index($seq, $re_seq, 0);
		while($re_pos != -1) {
			if($re_pos != 0
			&& $re_pos - $prev_pos < $max_frag_len) {
				print FRAGS "$fid\t$chrom\t".($f_offset + $prev_pos)."\t".($r_offset+$re_pos)."\n";
				print FENDS (2*$fid)."\t$chrom\t".($f_offset + $prev_pos)."\n";
				print FENDS (2*$fid+1)."\t$chrom\t".($r_offset + $re_pos - 1)."\n";
				my($na) = 0;
				my($nc) = 0;
				my($ng) = 0;
				my($nt) = 0;
				for(my($pos_i) = $prev_pos; 
				    $pos_i < $re_pos; 
				    $pos_i++) {
					my($c) = uc(substr($seq, $pos_i, 1));
					if($c eq "A") {
						$na++;
					}elsif($c eq "C") {
						$nc++;
					} elsif($c eq "G") {
						$ng++;
					} elsif($c eq "T") {
						$nt++;
					}
				}
				my($tot) = $na+$nc+$ng+$nt;
				if($tot == 0) {
					$tot = 1;
				}
				print FRAGS_CONT "$fid\t$chrom\t".($f_offset + $prev_pos)."\t".($r_offset+$re_pos)."\t$na\t$nc\t$ng\t$nt\n";
				print FRAGS_GC "chr$chrom\t".($f_offset + $prev_pos)."\t".($r_offset+$re_pos)."\t".(int(1000*($nc+$ng)/($tot))/1000)."\n";
				print FRAGS_FLEN "chr$chrom\t".($f_offset + $prev_pos)."\t".($r_offset+$re_pos)."\t".($re_pos-$prev_pos)."\n";
			}
			$prev_pos = $re_pos;
			$re_pos = index($seq, $re_seq, $re_pos+1);
			$fid++;
		}
		print STDERR "..frag $fid\n";
	}
}

sub get_fid_coord($$) {
	my($self, $fid) = @_;
	return($self->{chrom}->[$fid],$self->{coord}->[$fid]);
}

sub get_fend_coord($$) {
	my($self, $fend) = @_;
	my($fid) = int($fend)/2;
	my($strand) = $fend % 2;
	return($self->{chrom}->[$fid],
		$self->{coord}->[$fid] + 
			($strand ? $self->{frag_len}->[$fid] : 0));
}


sub map_frag($$$$) {
	my($self, $chrom, $coord, $scope) = @_;
	$chrom =~s/chr//;
	my($binsize) = $self->{binsize};
	my($bin) = int($coord/$binsize);
	my($fid) = -1;
	my($off) = 0;
	my($max_dev) = int($scope/$binsize);

	my($hit_bin) = -1;

	for(my($dev) = 0; $dev < $max_dev; $dev++) {
		my($dbin) = $bin + $dev;
		if(exists($self->{coor_to_frag}->{"$chrom\t$dbin"})) {
			$hit_bin = $dbin;
			last;
		}
		$dbin = $bin - $dev;
		if(exists($self->{coor_to_frag}->{"$chrom\t$dbin"})) {
			$hit_bin = $dbin;
			last;
		}
	}
	my($side) = "U";
	if($hit_bin != -1) {
		($fid,$side) = split("\t", 
			$self->{coor_to_frag}->{"$chrom\t$hit_bin"});
		$off = $side eq "5" ?
			($coord - $self->{coord}->[$fid]) :
			($self->{coord}->[$fid]+ $self->{frag_len}->[$fid] - $coord);
	}
	return(($fid, $off, $side));
}

sub frag_5_dist($$$) {
	my($self, $fid, $coord) = @_;

	my($c5) = $self->{coord}->[$fid];
	return($coord - $c5);
}
sub frag_3_dist($$$) {
	my($self, $fid, $coord) = @_;

	my($c3) = $self->{coord}->[$fid] + $self->{frag_len}->[$fid];
	return($c3 - $coord);
}

sub calc_multi_re_frag_lengths {
	my($self, $chrom_key, $re_seqs, $binsize, $max_len, $out_fn) = @_;

	my(@re_list) = split(/,/, $re_seqs);
	my(%len_counts);
	my($len,$bin);

	open(KEY, $chrom_key) || die "cannot open chrom key\n";

	while(<KEY>) {
		chop;
		my($chrom, $fn) = split("\t", $_);
		print STDERR "will proc $chrom\n";
		open(SEQ, $fn) || die "cannot opem crhom $chrom seq at $fn\n";
		my($seq) = "";
		<SEQ>;
		while(<SEQ>) {
			chop;
			$seq .= uc($_);
		}
		print STDERR "done read ".length($seq)." bp\n";

		my($prev_pos) = 0;
		my($re_pos) = get_next_re_pos(\$seq, 0, \@re_list);

		while($re_pos != -1) {
			if($re_pos != 0) {
				$len = $re_pos - $prev_pos + 1;
				$len = $len > $max_len ? $max_len : $len;
				$bin = int($len/$binsize) * $binsize;
				$len_counts{$bin}++;
			}
			$prev_pos = $re_pos;
			$re_pos = get_next_re_pos(\$seq, $re_pos+1, \@re_list);
		}
	}

	open(OUT, ">", $out_fn);
	print OUT "segment_length\tcount\n";
	foreach $bin (sort { $a <=> $b} keys %len_counts) {
		print OUT "$bin\t".$len_counts{$bin}."\n";
	}
	close(OUT);
}

sub get_next_re_pos($$@)
{
	my ($seq_ref,$pos,$re_ref) = @_;
	my $npos = 10**10;
	my $tpos = -1;
	foreach (@$re_ref) {
		$tpos = index($$seq_ref, $_, $pos);
		$npos = $npos < $tpos ? $npos : $tpos;
	}
	return ($npos);
}

1;
