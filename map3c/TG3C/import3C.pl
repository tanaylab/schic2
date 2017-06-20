use strict;
use File::Path("make_path");
require map3c::ConfOpt;

#valid configuration option names

my($opt) = ConfOpt::new("ConfOpt");

if(-e "3c.conf") {
	$opt->read_conf("3c.conf");
} elsif(-e ".3c.conf") {
	$opt->read_conf(".3c.conf");
}

$opt->parse_argv_conf(\@ARGV);

my($done) = "";
my($ret) = 0;

#========================
#Construct working directory

my($dir_lscripts) = $opt->get_opt("TG3C.lscripts");
my($wd) = $opt->get_opt("TG3C.workdir");
my($exp_nm) = $opt->get_opt("TG3C.3C_exp_nm");
my($foc_ndx) = $opt->get_opt("TG3C.foc_ndx", "");

my ($overwrite) = $opt->get_opt("TG3C.overwrite_existing", 0);

if(-e "$wd/done_ok") {
	unlink("$wd/done_ok");
}

if(!-e "$wd/$exp_nm.$foc_ndx") {
	make_path("$wd/$exp_nm.$foc_ndx");
}
$wd = "$wd/$exp_nm.$foc_ndx";

if(!-e "$wd/logs") {
	make_path("$wd/logs");
}

if(!-e "$wd/conf") {
	make_path("$wd/conf");
}

if(!-e "$wd/track_vars") {
	make_path("$wd/track_vars");
}

#=========================
# Stripping reads is recommended only in applications that work with 
# libraries of integrated constructs in which some non-genomic DNA must
# be
my($did_strip) = 0;
if($opt->get_opt("TG3C.do_strip", 0)) {
	print STDERR "\n3CPipe: will do stripping\n";
	if ($foc_ndx eq "") { die "\n======================\nERROR: stripping failed, TG3C.foc_ndx is not defined in conf file\n"; }
	my($fastq_dir) = $opt->get_opt("TG3C.3C_dir_$foc_ndx");
	my($fastq_fn_regexp) = $opt->get_opt("TG3C.3C_fn_regexp_$foc_ndx");
	my($seq_pref) = $opt->get_opt("TG3C.strip_seq_pref");
	my($mode) = "opt";
	if($seq_pref eq "NAMED") {
		$mode="named";
		my($nm) = $opt->get_opt("TG3C.bait_$foc_ndx");
		$seq_pref = $opt->get_opt("TG3C.strip_seq_pref_$nm");
	}
	my($is_paired) = $opt->get_opt("TG3C.is_paired");
	my($keep_pref) = $opt->get_opt("TG3C.keep_pref");
	my($keep_pref2) = $opt->get_opt("TG3C.keep_pref2",1);
	
	my @seq_pref = split(/,/, $seq_pref);
	my $seq_pref_regex = "(^" .  join("|^", @seq_pref) . ")(.+)";
	
	open(CONF, ">$wd/strip.conf");
	print CONF "inp_base=$fastq_dir\n";
	print CONF "inp_re=$fastq_fn_regexp\n";
	print CONF "out_fn1=$wd/strip_R1.fastq\n";
	print CONF "inp_r1_re=R1\n";
	print CONF "re1=\$R1=~/$seq_pref_regex/i\n";
	if($keep_pref) {
		print CONF "O1=\$O1=\$S[0].\$S[1]\n";
	} else {
		print CONF "O1=\$O1=\$S[1]\n";
	}
	if($is_paired) {
		print CONF "out_fn2=$wd/strip_R2.fastq\n";
		print CONF "inp_r2_re=R2\n";
		print CONF "r2=1\n";
		my($seq2_pref) = $opt->get_opt("TG3C.strip_seq_pref2", "NA");
		if($mode eq "named") {
			my($nm) = $opt->get_opt("TG3C.bait2_$foc_ndx");
			$seq2_pref = $opt->get_opt("TG3C.strip_seq_pref2_$nm");
		}
		if($seq2_pref ne "NA") {
			print CONF "re2=\$R2=~/^($seq2_pref)(.+)\$/\n";
			if($keep_pref2) {
				my($bb) = $b+1;
				print CONF "O2=\$O2=\$S[$b].\$S[$bb]\n";
			} else {
				print CONF "O2=\$O2=\$S[2]\n";
			}
		} else {
			print CONF "O2=\$O2=\$R2\n";
		}
	}
	close CONF;
	system("perl $dir_lscripts/TGPipeLib/fastq_splicer.pl \@$wd/strip.conf");
	 
	if(!-e "$wd/strip_R1.fastq") {
		die "\n======================\nERROR: stripping failed, output does not exist\n";
	}
	if(-z "$wd/strip_R1.fastq") {
		die "\n======================\nERROR: stripping failed - output empty\n";
	}
	$did_strip = 1;
}

#seg on RE
if($opt->get_opt("TG3C.do_seg", 1)) {
	print STDERR "\n3CPipe: will do segmentation\n";
	my($fastq_dir) = $opt->get_opt("TG3C.3C_dir_$foc_ndx");
	my($fastq_fn_regexp) = $opt->get_opt("TG3C.3C_fn_regexp_$foc_ndx", "");
	my($first_end_code) = $opt->get_opt("TG3C.fastq_first_code", "R1");
	my($second_end_code) = $opt->get_opt("TG3C.fastq_second_code", "R2");
	if($fastq_fn_regexp eq "") {
		$fastq_fn_regexp  = ".";
	}
	if($did_strip) {
		$fastq_dir = "$wd/";
		$fastq_fn_regexp = "strip";
	}
	print STDERR "seg $fastq_dir regexp $fastq_fn_regexp\n";
	my($RE_seq) = $opt->get_opt("TG3C.RE_seq");
	my($min_len) = $opt->get_opt("TG3C.segment_min_len");
	my $use_trim_split = $opt->get_opt("TG3C.trim_and_filter_low_qual_reads", 0);
	my $min_avg_read_qual = $opt->get_opt("TG3C.min_avg_read_qual", 20);
	my $min_base_qual = $opt->get_opt("TG3C.min_base_qual", 30);
	
	my($cmd) = "perl $dir_lscripts/TG3C/split_fastq_on_re.pl $fastq_dir $fastq_fn_regexp $wd/split.fastq $RE_seq $min_len 1 $first_end_code $second_end_code 2>$wd/logs/seg.log";
	
	if ($use_trim_split) {
		$cmd = "perl $dir_lscripts/TG3C/trim_fastqs_by_qual.pl --scope_dir $fastq_dir --regexp $fastq_fn_regexp --out_fn $wd/split.fastq --re_seq $RE_seq --min_frag_len $min_len --min_base_qual $min_base_qual --reads_to_sample 10000 --min_len 30 --min_read_qual $min_avg_read_qual 2>$wd/logs/seg.log";
	}
	
	print STDERR "\n$cmd\n";
	$ret = system($cmd);
	if($ret || !-e "$wd/split.fastq"
	|| -z "$wd/split.fastq") {
		die "\n======================\nERROR: segmentation on RE seq failed\n";
	}
}

#map
if($opt->get_opt("TG3C.do_map", 1)) {
	print STDERR "3CPipe: bowtie mapping\n";

	my($fq_lst);
	my($fq_lst_r2);
	my($map_inp) = $opt->get_opt("TG3C.map_inp");
	if($map_inp eq "split") {
	 	my(@fns) = <$wd/split.fastq>;
		$fq_lst = join(",", @fns);
	} else {
		my (@fns);
		if ($map_inp eq "fq_dir") {
			my($fastq_dir) = $opt->get_opt("TG3C.3C_dir_$foc_ndx");
			my($fastq_fn_regexp) = $opt->get_opt("TG3C.3C_fn_regexp_$foc_ndx", "");
			
			print STDERR "looking for: $fastq_dir/*${fastq_fn_regexp}*\n";
			@fns = <$fastq_dir/*${fastq_fn_regexp}*>
		}
		else {
			@fns = <$map_inp/*>;
		}
		my($fq_re1) = $opt->get_opt("TG3C.map_inp_regexp1");
		my($fq_re2) = $opt->get_opt("TG3C.map_inp_regexp2");
		my($fn);
		for $fn (@fns) {
			if($fn =~/$fq_re1/) {
				$fq_lst .= ",$fn";
			}
			if($fn =~/$fq_re2/) {
				$fq_lst_r2.= ",$fn";
			}
		}
	}
	my($bt2_bin) = $opt->get_opt("TG3C.bowtie2_bin");
	my($bt2_ndx) = $opt->get_opt("TG3C.bowtie2_ndx");
	my($nthreads) = $opt->get_opt("TG3C.bowtie2_threads");
	my $align_mode = $opt->get_opt("TG3C.bowtie2_align_mode", "local");
	
	my($sam_fn1) = "$wd/read1.sam";
	print STDERR "invoking bowtie2\n";
	my($cmd) = "$bt2_bin -p $nthreads --reorder --$align_mode -x $bt2_ndx -U $fq_lst -S $sam_fn1 2>$wd/logs/bt.log";

	if (! -e $sam_fn1 or $overwrite) {
		print STDERR "\n$cmd\n";
		$ret = system($cmd);
		if ($ret) { die "\nERROR:\n===================\nbowtie mapping failed\n"; }
	}
	if(defined($fq_lst_r2)) {
		my($sam_fn2) = "$wd/read2.sam";
		if (! -e $sam_fn2 or $overwrite) {
			print STDERR "invoking bowtie2 on read2\n";
			$cmd = "$bt2_bin -p $nthreads --reorder --$align_mode -x $bt2_ndx -U $fq_lst_r2 -S $sam_fn2 2>$wd/logs/bt.log";
			print STDERR "$cmd\n";
			$ret = system($cmd);
			if ($ret) { die "\nERROR:\n===================\nbowtie mapping failed\n"; }
		}
	}
}

# direct sam to adj (select either this or do_fendchain + do_adj)

if($opt->get_opt("TG3C.do_sam2adj", 1)) {
	print STDERR "3CPipe: direct sam to adj\n";

	open(SAMIDX, ">$wd/conf/sam2adj.idx") or die ("failed to open $wd/conf/sam2adj.idx");
	print SAMIDX "$wd/read1.sam\t$wd/read2.sam\n";
	close SAMIDX;

	my($dir_redb) = $opt->get_opt("TG3C.redb");
	my($re) = $opt->get_opt("TG3C.RE_seq");
	my($min_qual) = $opt->get_opt("TG3C.min_qual");

	#open (QYL, ">>$ENV{PIPELINE_HOME}/q_cmds");
	my($cmd) = "perl $dir_lscripts/TG3C/sampair2adj.pl $dir_redb $re $wd/conf/sam2adj.idx $wd/adj naive T $min_qual 2>$wd/logs/sampair2adj.log";
	#print QYL "$cmd\n";
	$ret = system($cmd);
	if ($ret) { die "\nERROR:\n===================\nsampair2adj.pl failed\n"; }

}

#chain
if($opt->get_opt("TG3C.do_fendchain", 0)) {
	print STDERR "3CPipe: sam to fendchain\n";
	my($dir_redb) = $opt->get_opt("TG3C.redb");
	my($re) = $opt->get_opt("TG3C.RE_seq");
	my($min_qual) = $opt->get_opt("TG3C.min_qual");
	my($min_len) = $opt->get_opt("TG3C.segment_min_len");
	my($calc_chain_mul) = $opt->get_opt("TG3C.calc_chains_multiplicity", 0);

	print STDERR "chaining\n";
	$ret = system("perl $dir_lscripts/TG3C/sam2chain.pl $wd/read1.sam $wd/chain $min_qual $min_len 1>$wd/chain 2>$wd/logs/sam2chain.log");
	if ($ret) { die "\nERROR:\n===================\nsam2chain.pl failed\n"; }
	print STDERR "chain to fends";
	my ($chain_mul_fn) = $calc_chain_mul ? "$wd/track_vars/chain_mul.Rtable" : "";
	$ret = system("perl $dir_lscripts/TG3C/chain2fendchain.pl $dir_redb $re $wd/chain $wd/fendchain $chain_mul_fn 2>$wd/logs/fendchain.log");
	if ($ret) { die "\nERROR:\n===================\nchain2fendchain.pl failed\n"; }
}

#adj
if($opt->get_opt("TG3C.do_adj", 0)) {
	print STDERR "\n3CPipe: fendchain to adj\n";
	my($switch_ratio) = $opt->get_opt("TG3C.switch_ratio");
	my($out_mode) = $opt->get_opt("TG3C.adj_out_mode", "no");
	print STDERR "out mode is $out_mode\n";
	my($filter_near_sonic) = $opt->get_opt("TG3C.filter_near_sonic");
	my $canonic = $opt->get_opt("TG3C.output_sorted_adj", 0);
	my($key_out_pref) = $opt->get_opt("TG3C.key_out_pref", "read_key");

	my $adjust_fc_by_key2fid = $opt->get_opt("TG3C.adjust_fc_by_key2fid", 0);
	my $adjust_fc_params = $adjust_fc_by_key2fid ? " $wd/${key_out_pref}_key2fid.txt $wd/${key_out_pref}.id2key" : "";

	my $cmd = "perl $dir_lscripts/TG3C/fendchain2adj.pl $wd/fendchain $wd/adj $switch_ratio $filter_near_sonic $canonic $out_mode $adjust_fc_params 2>$wd/logs/fendchain2adj.log";
	print STDERR "\n$cmd\n";
	$ret = system($cmd);
	if ($ret) { die "\nERROR:\n===================\nfendchain2adj.pl failed\n"; }
	if ($out_mode =~ /full/) {
		$out_mode = "full";
		print STDERR "Will also create adj.full file\n";
		$ret = system("perl $dir_lscripts/TG3C/fendchain2adj.pl $wd/fendchain $wd/adj.full $switch_ratio $filter_near_sonic $canonic $out_mode $adjust_fc_params 2>$wd/logs/fendchain2adj.full.log");
		if ($ret) { die "\nERROR:\n===================\nfull fendchain2adj.pl failed\n"; }
	}
}
if($opt->get_opt("TG3C.do_adj_coord", 0)) {
	print STDERR "\n3CPipe: adj 2 coord\n";
	my($dir_redb) = $opt->get_opt("TG3C.redb");
	my($re) = $opt->get_opt("TG3C.RE_seq");
	if (-e "$wd/adj.full") {
		$ret = system("perl $dir_lscripts/TG3C/adj2coords.pl $dir_redb $re $wd/adj,$wd/adj.full 2>$wd/logs/adj2coord.log");
	} else {
		$ret = system("perl $dir_lscripts/TG3C/adj2coords.pl $dir_redb $re $wd/adj 2>$wd/logs/adj2coord.log");
	}
	if ($ret) { die "\nERROR:\n===================\nadj2coords.pl failed\n"; }
}

if(-e "$wd/adj") {
	open(OK, ">$wd/done_ok");
	print OK "ok";
	close OK;
}
