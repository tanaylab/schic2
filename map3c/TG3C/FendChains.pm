#support for analysis of bis amplicons and extraction of methylation state
package TG3C::FendChains;

use strict;
use Data::Dumper;
 
1;
#use diagnostics;

sub new($) {
	my($clas) = @_;

	my($self) = {};
	
	bless($self,$clas);

#map fend chains, with read counts
	$self->{fcs} = {};

#overall coverage per fid
	$self->{cov} = {};

#fids that are adj to each fid, with mol count, read count
	$self->{fid_adj} = {};

	$self->{pairs_r} = {};
	$self->{pairs_m} = {};
	$self->{pairs_direct} = {};
	$self->{pairs_indirect} = {};
	$self->{pairs_strandmis} = {};
	$self->{pairs_max_r} = {};
	$self->{pairs_all_r} = {}; #count also filtered reads as total_reads
	$self->{foc_bait_fend};
	return($self);
}

sub read_fendchain {
	my($self, $fc_tab, $read_count, $foc_fend, $key2fid_fn, $id2key_fn) = @_;

	my $filter_bait = ($#_ > 2 and $foc_fend != -1) ? 1 : -1 ;

	my($key2fid);
	my($id2key);

	if ($key2fid_fn) {
		$key2fid = tab2hash($key2fid_fn, "key");
		die "Missing id2key file." if (! -e $id2key_fn);
		$id2key = tab2hash($id2key_fn, "id");
	}

	open(FC, $fc_tab) || die "cannot open fendchain tab $fc_tab\n";
	my($fcs) = $self->{fcs};
	my($count) = 0;

	while(<FC>) {
		chop;
		# if filtering by bait fend - skip line for irrelevant fendchains
		if($filter_bait == 1){
			my $valid = find_bait($_, $foc_fend);
			if($valid != -1){
				$self->{foc_bait_fend} = $valid;
			} else { 
				next;
			}
		} 

		my(@f) = split("\t", $_);

		my $id = shift(@f);

		if (defined($id2key)) {
			if (!defined($id2key->{$id})) {
				die "\nKey for id $id not found in $id2key_fn\n";
			}
			
			my($hkey) = $id2key->{$id};
			
            # in this mode, skipping fendchains without a matching fid
			next if (!defined($key2fid->{ $hkey->{key} }));

			my $hrep_fid = $key2fid->{ $hkey->{key} }; 
			my $lfc = shift(@f);
			my $maxpos = $hrep_fid->{frag_len} > 50 ? 50 : $hrep_fid->{frag_len};
			unshift(@f, sprintf("%d:1-:%d:0:%d:%d:42", $hrep_fid->{frag}, $maxpos, $maxpos, $hrep_fid->{frag_len} - $maxpos));
			unshift(@f, $lfc + 1);
		}

		unshift(@f, $id);

		my(@f_rc) = fcs_rc(@f);

		if($f[1] == 1) {
			next;
		}
		if($count % 100000 == 0) {
			print STDERR "count ".int($count/100000)."M..";
		}
		$count++;	
		# TODO:: remove this mols?
		my($mols) = 1;
		my($reads) = 1;
		my($maxr) = 1;
		
		for(my($i) =2; $i <= $#f; $i++) {
			$f[$i] =~s/:\d+$//;
			$f_rc[$i] =~ s/:\d+$//;
		}
		
		my($key) = join("\t", @f[2..$#f]);
		if($filter_bait == 1){
			$key =  filter_nondig($key, 2);
			next if(length($key)==0);
		}
		
		my($rc_key) = join("\t", @f_rc[2..$#f]);
	
		#Avoid counting reverse comp as mols but count them as reads 
		if (exists($fcs->{$rc_key})) {
			$fcs->{$key} = $fcs->{$rc_key};
			$maxr = $fcs->{$key};
			delete $fcs->{$rc_key};
		}

		if($read_count) {
			$reads = pop @f;
			$maxr = $reads;
		} else {
			if(exists($fcs->{$key})) {
				# TODO:: remove this mols?
				$mols = 0;
				$fcs->{$key}++;
				$maxr = $fcs->{$key};
			} else {
				$fcs->{$key} = 1;
			}
		}

		#test unique

		if($filter_bait==1){
			my @freads = split("\t", $key);
			@f = ($f[0], $#freads+1, @freads);
		}
		for(my($i) = 2; $i < $#f; $i++) {
			my(@fe1) = split(":", $f[$i]);
			my(@fe2) = split(":", $f[$i+1]);
			if($fe1[0] eq "-1" || $fe2[0] eq "-1") {
				next;
			}
			# transform the fids of each fragment into fends
			my($fe_id1) = $fe1[0]*2 + (($fe1[1]=~/\+/) ? 1 : 0);
			my($fe_id2) = $fe2[0]*2 + (($fe2[1]=~/\+/) ? 0 : 1);

			my($fpr) = "$fe_id1\t$fe_id2";

			if(int($self->{pairs_max_r}->{$fpr}) < $maxr) {
				$self->{pairs_max_r}->{$fpr} = $maxr;
			}
		}
	}
	print STDERR "fc2ad: Done\n";
	print STDERR "fc2ad: fcs $count (at least 2 segs)\n";

}

sub fcs_rc {	
	my @f = @_;
	@f[2..$#f] = reverse @f[2..$#f];
	for(my($i) =2; $i <= $#f; $i++) {
			my @flds = split(/:/, $f[$i]);
			$flds[1] =~ tr/12\+\-/21\-\+/;
			$flds[2] = $flds[2];
			my $tmp = $flds[3];
			$flds[3] = $flds[5];
			$flds[5] = $tmp;
			$flds[4] = $flds[2]+$flds[3];
			$f[$i] = join(':', @flds);
	}

	return(@f);
}
sub gen_filt_adj($$$) {
	my($self, $switch_ratio, $filt_pair_pad, $canonic) = @_;

	if(!defined($filt_pair_pad)) {
		$filt_pair_pad = -1;
	}

	my($adj) = $self->{fid_adj};
	my($fcs) = $self->{fcs};

	print STDERR "into adj generation switch filt $switch_ratio\n";
	my($mcount) = 0;
	my($filt_cnt) = 0;
	my($filt_pair) = 0;
	my($fc);
	my(%adj_keys);
	foreach $fc (keys %$fcs) {
		chop;
		my(@f) = split("\t", $fc);

		if($mcount % 100000 == 0) {
			print STDERR "mcnt ".int($mcount/1000)."K..";
			print STDERR "filtpair ".int($filt_pair/1000)."K..";
			print STDERR "filtswitch ".(int($filt_cnt/1000))."K\n";
		}
		my($mols) = 1;
		my($reads) = $fcs->{$fc};

		#test unique
		for(my($i) = 0; $i < $#f; $i++) {
			my(@fe1) = split(":", $f[$i]);
			my(@fe2) = split(":", $f[$i+1]);
			if($fe1[0] eq "-1" || $fe2[0] eq "-1") {
				next;
			}

			my($fe_id1) = $fe1[0]*2 + (($fe1[1]=~/\+/) ? 1 : 0);
			my($fe_id2) = $fe2[0]*2 + (($fe2[1]=~/\+/) ? 0 : 1);
			my($fpr) = "$fe_id1\t$fe_id2";
			#if canonic is true; output sorted fend pairs list (NOT FOR 4C)
			if ($fe_id2 > $fe_id1 && $canonic) {
				 my (@tmp) = @fe1;
				 @fe1 = @fe2;
				 @fe2 = @tmp;
				 $fpr = "$fe_id2\t$fe_id1";
			 }
			$self->{pairs_all_r}->{$fpr} += $reads; #Count ALL reads for fendpair	
			if(int($self->{pairs_max_r}->{$fpr}) * $switch_ratio > $reads) {
				$filt_cnt++;
				next;
			}
			if($filt_pair_pad != -1) {
				my($x) = $fe1[3];
				my($y) = $fe2[5];
				my($hit) = 0;
				for(my($o1) = -$filt_pair_pad; $o1 <= $filt_pair_pad; $o1++) {
					for(my($o2) = -$filt_pair_pad; $o2 <= $filt_pair_pad; $o2++) {
						my($pairkey) = "$fpr\t".($x+$o1)."\t".($y+$o2);
						if(exists($adj_keys{$pairkey})) {
							$hit = 1;
							last;
						}
					}
					if($hit) {
						last;
					}
				}
				if($hit) {
					$filt_pair++;
					next;
				}
				for(my($o1) = -$filt_pair_pad; $o1 <= $filt_pair_pad; $o1++) {
					for(my($o2) = -$filt_pair_pad; $o2 <= $filt_pair_pad; $o2++) {
						my($pairkey) = "$fpr\t".($x+$o1)."\t".($y+$o2);
						$adj_keys{$pairkey} = 1;
					}
				}
			}

			if(!exists($self->{pairs_m}->{$fpr})) {
				if(!exists($adj->{$fe_id1})) {
					$adj->{$fe_id1} = [];
				}	
				push(@{$adj->{$fe_id1}}, $fe_id2);
			}
			$self->{pairs_r}->{$fpr}+= $reads;
			$self->{pairs_m}->{$fpr}+= $mols;
			$mcount++;

			my($end1) = substr($fe1[1],0,1);
			my($end2) = substr($fe2[1],0,1);
			my($strand1) = substr($fe1[1],1,1);
			my($strand2) = substr($fe2[1],1,1);
			if(!($end1 eq "1" && $end2 eq "2")) {
				if($strand1 ne $strand2) {
					$self->{pairs_strandmis}->{$fpr}+= $mols;
				} else {
					$self->{pairs_direct}->{$fpr}+= $mols;
					if($i+2 <= $#f) {
						my(@fe3) = split(":", $f[$i+2]);
						my($fe_id3) = $fe3[0]*2 + (($fe3[1]=~/\+/) ? 0 : 1);
						if($fe3[0] ne "-1" 
						&& $fe1[0] ne $fe2[0]
						&& $fe2[0] ne $fe3[0]) {
							$self->{pairs_indirect}->{"$fe_id1\t$fe_id3"}+= $mols;
						}
					}
				}
			}
		}
	}
	print STDERR "fc2ad: Done\n";
	print STDERR "fc2ad: unique $mcount\n";
	print STDERR "fc2ad: filt $filt_cnt\n";
	print STDERR "fc2ad: filt pair $filt_pair\n";
}

sub write_pair_cov {
	my($self, $pair_fn, $print_bait_id) = @_;
	# tot_reads includes filtered
	my $header = "fend1\tfend2\tmol\tdirect\tindirect\tstrandmis\treads\tpair_max_r\ttot_reads";
		
	if($print_bait_id == -1) {
		open(OUT, ">$pair_fn") || die "cannot write pair fn $pair_fn\n";
		print OUT "$header\n";
	} else {
		open(OUT, ">>$pair_fn") || die "cannot write pair fn $pair_fn\n";
		if(-s $pair_fn ==0){ 
			print OUT "$header\tbait_fend\n";
		}
	}

	my($k);
	foreach $k (keys %{$self->{pairs_m}}) {
		next if ( $print_bait_id && (split("\t", $k))[0] == $self->{foc_bait_fend});
		print OUT "$k\t".($self->{pairs_m}->{$k}).
			"\t".int($self->{pairs_direct}->{$k}).
			"\t".int($self->{pairs_indirect}->{$k}).
			"\t".int($self->{pairs_strandmis}->{$k}).
			"\t".($self->{pairs_r}->{$k}).
			"\t".($self->{pairs_max_r}->{$k}).
			"\t".($self->{pairs_all_r}->{$k});
			if($print_bait_id != -1) {
			         print OUT "\t$self->{foc_bait_fend}\n";
			} else {
				print OUT "\n";
			}
	}
}

sub write_pairs_for_track {
	my($self, $pair_fn, $use_reads) = @_;

	if(!defined($use_reads)) {
		$use_reads = 0;
	}

	open(OUT, ">$pair_fn") || die "cannot write pair fn $pair_fn\n";
	print OUT "fend1\tfend2\tcount\n";
	my($k);
	my(%block);
	foreach $k (keys %{$self->{pairs_m}}) {
		if(exists($block{$k})) {
			next;
		}
		
		my($fe1, $fe2) = split("\t", $k);
		my($m);
		if($use_reads) {
			$m = $self->{pairs_r}->{$k};
		} else {
			$m = $self->{pairs_m}->{$k};
		}
		if(exists($self->{pairs_m}->{"$fe2\t$fe1"})) {
			if($use_reads) {
				$m += $self->{pairs_r}->{"$fe2\t$fe1"};
			} else {
				$m += $self->{pairs_m}->{"$fe2\t$fe1"};
			}
			$block{"$fe2\t$fe1"} = 1;
		}
		print OUT "$k\t$m\n";
	}
}

sub tab2hash {
	my($tab, $key_column) = @_;
	open(IN, $tab) or die "Failed to open $tab for reading.\n";

	my $header = <IN>;
	chop $header;
	my @h = split(/\t/, $header);
	my $kc = -1;
	for (my $i = 0; $i <= $#h; $i++) {
		if ($h[$i] eq $key_column) {
			if ($kc >= 0) {
				die "Multiple appearance of key column $key_column in $tab ($kc and $i).\n";
			}
			$kc = $i;
		}
	}
	my $o;
	my $k;

	while(<IN>) {
		chop;
		my(@f) = split("\t", $_);
		$o->{ $f[$kc] } = {};
		for (my $i = 0; $i <= $#f; $i++) {
			next if $i == $kc;
			$o->{ $f[$kc] }->{ $h[$i] } = $f[$i];
		}
	}
	close(IN);
	return $o;
}

sub find_bait 
# fendchains are actually in the fid space
{
	my($fcline, $bait_fend) = @_;
	my @f = split("\t", $fcline);
	my @fends;
	my @strands;

	for(my $i = 2; $i<= $#f; $i++){
		push(@fends,   (split(":", $f[$i]))[0]);
		push(@strands, (split(":", $f[$i]))[1]);
	}
	# transform to fends
	for( my $i=0; $i<=$#fends; $i++){
			$fends[$i] = $fends[$i]*2 + (($strands[$i]=~/\+/) ? 1 : 0);
	}
	# filter out pairwise contacts and fendchains not associated with bait
	return(($bait_fend eq $fends[0] & $#fends > 1 ? $bait_fend : -1));
}

sub filter_nondig
{
	# fcline is in fids
	my($fcline, $f_delta) = @_;
	my @f = split("\t", $fcline);
	my @fids;
	my @v_fids;
	my @deltas;
	for(my $i = 0; $i<= $#f; $i++){
		push(@fids,(split(":", $f[$i]))[0]);
	}
	for(my $i=0; $i < $#f; $i++){
		my $delta = abs($fids[$i]-$fids[$i+1]);
		push(@deltas, $delta);
		if($delta > $f_delta){
			push(@v_fids, $i);
			if($i == $#f-1){ push(@v_fids, $i+1)}
		} 
	}
	my @v_key = @f[@v_fids];
	return(join("\t",@v_key));
}

