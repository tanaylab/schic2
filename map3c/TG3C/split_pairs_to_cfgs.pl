
use strict;

if ($#ARGV < 2) {
	die "usage: split_pairs_to_confs.pl <fastq dir> <ofn_pref> <n per conf> [ <pref name> suff_line1 suff_line2 ... ]\n\tReads fastq_regexp - idx_name pairs from STDIN, create config files for gtrack.2d.create_from_3Cseq with n_per_conf fastq pairs per config file\n";
}

my $fastq_dir = shift @ARGV;
my $out_pref = shift @ARGV;
my $n_per_conf = shift @ARGV;
my $pref_name = $#ARGV >= 0 ? shift @ARGV : "";

my $line;
my @idxs;
my $p2;
my ($suff) = "\n\n#include ". $ENV{PIPELINE_HOME} . "/map3c/config/scell_shared.conf";

if ($pref_name ne "") {
	$suff = "$suff\nTG3C.3C_exp_nm=$pref_name\n";
}
if ($#ARGV >= 0) {
	$suff .= join("\n", @ARGV);
}

my $n = 0;

while($line = <STDIN>) {
	chop $line;
	my ($fre, $name) = split('\t', $line);
	
	if ($n % $n_per_conf == 0) {
		if ($#idxs >= 0) {
			print OUT "TG3C.3C_indices=". join(',', @idxs) ."\n\n$p2\n$suff\n";
			close(OUT);
		}

		open(OUT, sprintf(">%s_%s.conf", $out_pref, $n / $n_per_conf + 1)) or die "failed to open output file";
		@idxs = ();
		$p2 = "";
	}
	push(@idxs, $name);
	$p2 .= sprintf("TG3C.3C_dir_%s=%s\nTG3C.3C_fn_regexp_%s=%s\n", $name, $fastq_dir, $name, $fre);
	$n++;
}
if ($#idxs >= 0) {
	print OUT "TG3C.3C_indices=". join(',', @idxs) ."\n\n$p2\n$suff\n";
	close(OUT);
}

