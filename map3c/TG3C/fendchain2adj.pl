use strict;

require map3c::TG3C::FendChains;

if($#ARGV < 3) {
	die "fendchain2adj fendchain adj_out switch_ratio filt_pair_pad sort_canonic(not 4C) [full]\n";
}

my($fc) = TG3C::FendChains::new("TG3C::FendChains");

my($switch_ratio) = $ARGV[2];
my($filt_pad) = $ARGV[3];
my ($sort_canonic) = $ARGV[4];
my($mode) = $ARGV[5];
my($key2fid) = $#ARGV == 7 ? $ARGV[6] : 0;
my($id2key) = $#ARGV == 7 ? $ARGV[7] : 0;


print STDERR "args sw $switch_ratio filt $filt_pad mode $mode key2fid $key2fid id2key $id2key\n";

if(!defined($mode)) {
	$mode = "";
}

print STDERR "will process $ARGV[0]\n";

$fc->read_fendchain($ARGV[0], 0, -1, $key2fid, $id2key);
$fc->gen_filt_adj($switch_ratio, $filt_pad, $sort_canonic);

if($mode eq "full") {
	$fc->write_pair_cov($ARGV[1]);
} else {
	if($mode eq "reads") {
		print STDERR "work in reads mode\n";
	}
	$fc->write_pairs_for_track($ARGV[1], $mode eq "reads" ? 1 : 0);
}
