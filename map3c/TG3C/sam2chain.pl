use strict;

#we read the SAM

sub resolve_hits {
	my($id, $hits1, $hits2) = @_;

#hits are tab delim string with the fields  $segi, chrom, from, strand, mapq
#we generate a chain string consisting of tab delimied fields: Rx-y:chr:fr:Q:R!2
#Rx-y is the location in the read, with and integrated R1 R2 space.
#Q is the mapq value. R12 is 1,2 or B (accordig to read coverage)


#if one hit - we do it quickly
	if($#$hits1 + $#$hits2 == -1) {
		my($end, $segi, $rfr, $rto, $chr, $pos, $strand, $qual) = split("\t", $#$hits1 == 0 ? $hits1->[0] : $hits2->[0]);
		return("1\tR$rfr-$rto:$chr:$pos:$strand:$qual:$end");
	}
	
#otherwise we go over the hits

#go over read 1 hits, store and index
	my($err) = undef;

	my(@r_rfr, @r_rto, @r_chr, @r_pos, @r_strand, @r_qual);

	for(my($i1) = 0; $i1 <= $#$hits1; $i1++) {
		my($hit) = $hits1->[$i1];
		my($end, $segi, $rfr, $rto, $chr, $pos, $strand, $qual) = split("\t", $hit);
		push(@r_rfr, $rfr);
		push(@r_rto, $rto);
		push(@r_chr, $chr);
		push(@r_pos, $pos);
		push(@r_strand, $strand);
		push(@r_qual, $qual);
		
		print STDERR "push R $rfr $rto at $chr:$pos:$strand\n";
	}
	print STDERR "pushed $#r_rfr elements to r1\n";

	my($max_r1_hit) = $#r_rfr;
	my($min_r2_hit) = $#r_rfr + 1;	#may be updated

	for(my($i2) = $#$hits2; $i2 >= 0; $i2--) {
		my($hit) = $hits2->[$i2];
		my($end, $segi, $rfr2, $rto2, $chr2, $pos2, $strand2, $qual2) = split("\t", $hit);
		my($maxpos2) = $pos2 + $rto2 - $rfr2;
		print STDERR "search for intersect R2 $rfr2 $rto2, $chr2:$pos2:$strand2\n";
		my($hit_at_1) = 0;
		for(my($i1) = 0; $i1 <= $max_r1_hit; $i1++) {
			if($chr2 ne $r_chr[$i1]) {
				print STDERR "diff chrom $chr2 $r_chr[$i1]\n";
				next;
			}
			my($pos1) = $r_pos[$i1];
			my($maxpos1) = $pos1 + ($r_rto[$i1] - $r_rfr[$i1]);
			print STDERR "test intersect i1 $pos1 $maxpos1 i2 $pos2 $maxpos2\n";
			
			if(($pos2 >= $pos1 && $pos2 <= $maxpos1)
			|| ($pos1 >= $pos2 && $pos1 <= $maxpos2)) {
				$hit_at_1 = 1;
				if($i2 != $#$hits2 && $min_r2_hit > $i1) {
					$err = "ERR\tr2_skip\t$chr2\t$pos2";
				} elsif($r_strand[$i1] eq $strand2) {
					$err = "ERR\tstrand_endmis\t$chr2\t$pos2";
				} else {
					$min_r2_hit = $i1;
					if($strand2 eq "-") {
						if($maxpos2 > $maxpos1) {
							print STDERR "Read2 - extend forward i1 $i1 delta ".($maxpos2-$maxpos1)."\n";
							$r_rto[$i1] += ($maxpos2-$maxpos1);
						}
						if($pos2 < $pos1) {
							print STDERR "Read2 - chewback? pos1 $pos1 pos2 $pos2\n";
							$r_rfr[$i1] -= ($pos1-$pos2);
							$r_pos[$i1] = $pos2; #OK??
						}
					} else {
						if($pos2 < $pos1) {
							print STDERR "Read2 + extend forward i1 $i1 delta ".($pos1-$pos2)."\n";
							$r_rto[$i1] += ($pos1-$pos2);
							$r_pos[$i1] = $pos2;
						}
						if($maxpos2 > $maxpos1) {
							print STDERR "Read2 + chewback? maxpos2 $maxpos2 maxpos1 $maxpos1\n";
							$r_rfr[$i1] -= ($maxpos2-$maxpos1);
						}
					}
					if($r_qual[$i1] < $qual2) {
						$r_qual[$i1] = $qual2;
					}
				}
				last;
			}
		}
		
		if(defined($err)) {
			last;
		}
		if(!$hit_at_1) {
			push(@r_chr, $chr2);
			push(@r_rfr, 1000-$rfr2);
			push(@r_rto, 1000-$rto2);
			push(@r_pos, $pos2);
			push(@r_strand, $strand2 eq "+" ? "-" : "+");
			push(@r_qual, $qual2);
		}
	}
	if(defined($err)) {
		print STDERR "$err\n";
		return($err);
	}

#first find the read2 to read1 mapping if exists
	my($read_offset) = undef;

	my($chain) = ($#r_chr+1);
#(sort by location?
#finally we go over the read and report all hits
	for(my($ri) = 0; $ri <= $#r_chr; $ri++) {
		my($end) = 2;
		if($ri <= $max_r1_hit) {
			if($ri >= $min_r2_hit) {
				$end = "B";
			} else {
				$end = "1";
			}
		}

		$chain .= "\tR$r_rfr[$ri]-$r_rto[$ri]:$r_chr[$ri]:$r_pos[$ri]:$r_strand[$ri]:$r_qual[$ri]:$end";
	}
	print STDERR "CHAIN\t$chain\n";

	return($chain);
}

if($#ARGV != 3 ){
	die "usage sam2chain sam_fn out_fn min_qual min_len\n";
}

# Notice that the sam file is expected to be sorted by alignment id (first column)
open(SAM, $ARGV[0]) || die "cannot open SAM input $ARGV[0])";
my($min_qual) = $ARGV[2];
my($min_len) = $ARGV[3];
my($cur_tag);
my($prev_id);
my($hits1) = [];
my($hits2) = [];
open(OUT, ">$ARGV[1]") || die "cannot write output chain\n";
my($id_i) = 0;
while(<SAM>) {
	chop;
	my(@f) = split("\t", $_);
	my($id, $end, $segi, $rfr, $rto) = $f[0]=~/(.+):R(.)SEG(.):(\d+):(\d+)/;
	if($id ne $prev_id) {
		if($#$hits1 + $#$hits2 >= -1) {
			$id_i++;
			my($chain) = resolve_hits($id_i, $hits1, $hits2);
			print OUT "$prev_id\t$chain\n";
		}
		$hits1 = [];
		$hits2 = [];
		$prev_id = $id;
	}
	if($f[1] != 4 && $f[4] > $min_qual && $rto-$rfr >= $min_len) {
		my($strand) = "+";
		if($f[1] & 16) {
			$strand = "-";
		}
		if($end eq "1") {
			push(@$hits1, "$end\t$segi\t$rfr\t$rto\t$f[2]\t$f[3]\t$strand\t$f[4]");
		} else {
			push(@$hits2, "$end\t$segi\t$rfr\t$rto\t$f[2]\t$f[3]\t$strand\t$f[4]");
		}
	}
}
if($#$hits1 + $#$hits2 >= -1) {
	$id_i++;
	my($chain) = resolve_hits($id_i, $hits1, $hits2);
	print OUT "$prev_id\t$chain\n";
}


#f0 seqid M00393:117:000000000-A6L8E:1:1101:16993:1442:R1SEG0:0:248	
#f1 code  4: unmapped, 16 reverse
#f2 chr1	
#f3 coord 144009694	
#f4 MAPQ 44	
#f5 CIGAR 248M 	
#f6 MRN/RNEXT *	
#f7 PNEXT 0	
#f8 TLEN 0	
#f9 SEQ 
#f10 QUAL 	
#f11- TAGS (XX:y:VAL)

