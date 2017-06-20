use strict;
require map3c::ReDB;
require map3c::TG3C::FendChains;

if($#ARGV < 2) {
	die "usage: adj2coord.pl redb_dir re adj_file1,adj_file2...\n";
}

my($re) = $ARGV[1];
my($in_files) = $ARGV[2];
my @in_files = split(/,/, $in_files);
#my($out_file) = $ARGV[3];

my $redb = ReDB::new("ReDB");
$redb->init_from_tabs($ARGV[0]."/$re.frags");
foreach my $in_file (@in_files) {
	#read fend tb:
	open(ADJ, $in_file) or die "cannot open $in_file\n";
	open(OUT, ">$in_file.coord")  or die "Cannot write $in_file.coord\n";



	my($col1) = 0;
	my($col2) = 1;

	if($#ARGV > 3) {
		$col1 = $ARGV[4];
		$col2 = $ARGV[5];
	}
	my($h);
	$h = <ADJ>;
	chop $h;
	my(@hd) = split("\t", $h);
	$hd[$col1] = "fend1\tchr1\tstart1";
	$hd[$col2] = "fend2\tchr2\tstart2";
	print OUT join("\t", @hd)."\n";

	while (<ADJ>) {
		chomp $_;
		my @f = split("\t", $_);

		my($fend1) = $f[$col1];
		my($fend2) = $f[$col2];

		my($chr1, $coord1) = $redb->get_fend_coord($fend1);
		my($chr2, $coord2) = $redb->get_fend_coord($fend2);
		if(!defined($chr1)) {
			print STDERR "undef chr at fend $fend1\n";
		}
		$f[$col1] = "$fend1\t$chr1\t$coord1";
		$f[$col2] = "$fend2\t$chr2\t$coord2";
		print OUT join("\t", @f)."\n";
	}


}

