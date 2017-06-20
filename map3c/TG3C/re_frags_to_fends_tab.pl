use strict;

my($re_dir) = $ARGV[0];
my($re_seq) = $ARGV[1];

open(GC, "$re_dir/$re_seq"."_gc") || die "Cannot open GC content table at $re_dir seq $re_seq\n";
open(FLEN, "$re_dir/$re_seq"."_flen") || die "Cannot open flen table at $re_dir seq $re_seq\n";
open(MAP, "$re_dir/$re_seq"."_map") || die "Cannot open mapab table at $re_dir seq $re_seq\n";
open(FE_GC, ">$re_dir/fe$re_seq"."_gc") || die "Cannot open GC content table at $re_dir seq $re_seq\n";
open(FE_FLEN, ">$re_dir/fe$re_seq"."_flen") || die "Cannot open flen table at $re_dir seq $re_seq\n";
open(FE_MAP, ">$re_dir/fe$re_seq"."_map") || die "Cannot open mapab table at $re_dir seq $re_seq\n";

my($l);
$l = <GC>;
print FE_GC $l;
while(<GC>) {
	chop;
	my(@f) = split("\t", $_);
	
	print FE_GC "$f[0]\t$f[1]\t".($f[1]+1)."\t$f[3]\n";
	print FE_GC "$f[0]\t".($f[2]-1)."\t$f[2]\t$f[3]\n";
}
$l = <FLEN>;
print FE_LEN $l;
while(<FLEN>) {
	chop;
	my(@f) = split("\t", $_);
	
	print FE_FLEN "$f[0]\t$f[1]\t".($f[1]+1)."\t$f[3]\n";
	print FE_FLEN "$f[0]\t".($f[2]-1)."\t$f[2]\t$f[3]\n";
}
$l = <MAP>;
print FE_MAP $l;
while(<MAP>) {
	chop;
	my(@f) = split("\t", $_);
	
	print FE_MAP "$f[0]\t$f[1]\t".($f[1]+1)."\t$f[3]\n";
	print FE_MAP "$f[0]\t".($f[2]-1)."\t$f[2]\t$f[3]\n";
}
