require map3c::ReDB;

die "usage: re_seq chrom_key out_dir\n" if $#ARGV < 2;

my($redb) = ReDB::new("ReDB");

$redb->build_from_genome($ARGV[1], $ARGV[0], 2, 2, 10000000, $ARGV[2]);

