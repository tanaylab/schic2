use strict;
require map3c::ReDB;
require map3c::TG3C::FendChains;

if($#ARGV < 2) {
	die "usage: fendchain_report.pl redb_dir re fendchain cis_cutoff ofn\n";
}

my($re) = $ARGV[1];
my($fc_fn) = $ARGV[2];
my($cis_cutoff) = $ARGV[3];
my($ofn) = $ARGV[4];

my $redb = ReDB::new("ReDB");
$redb->init_from_tabs($ARGV[0]."/$re.frags");

my($fc) = TG3C::FendChains::new("TG3C::FendChains");

$fc->read_fendchain($fc_fn, 0, -1, 0, 0);

my($fcs) = $fc->{fcs};

print STDERR "into fendchain report: #trans/far_cis/near_cis adjs per fendchain\n";

my $counts = {};
my($n_close, $n_far, $n_trans);
foreach my $fc (keys %$fcs) {		
	chop;
	my(@f) = split("\t", $fc);
	
	$n_close = $n_far = $n_trans = 0;
	
	for(my($i) = 0; $i < $#f; $i++) {
		my(@fe1) = split(":", $f[$i]);
		my(@fe2) = split(":", $f[$i+1]);
		if($fe1[0] eq "-1" || $fe2[0] eq "-1") {
			next;
		}
		
		my($fe_id1) = $fe1[0]*2 + (($fe1[1]=~/\+/) ? 1 : 0);
		my($fe_id2) = $fe2[0]*2 + (($fe2[1]=~/\+/) ? 0 : 1);

		my($chr1, $coord1) = $redb->get_fend_coord($fe_id1);
		my($chr2, $coord2) = $redb->get_fend_coord($fe_id2);		

		if ($chr1 ne $chr2) {
			$n_trans++;
		}
		elsif (abs($coord1 - $coord2) <= $cis_cutoff) {
			$n_close++;
		}
		else {
			$n_far++;
		}
	}
	if(!exists($counts->{$n_trans})) {
		$counts->{$n_trans} = {};
	}
	if(!exists($counts->{$n_trans}->{$n_close})) {
		$counts->{$n_trans}->{$n_close} = {};
	}
	if(!exists($counts->{$n_trans}->{$n_close}->{$n_far})) {
		$counts->{$n_trans}->{$n_close}->{$n_far} = 0;
	}
	$counts->{$n_trans}->{$n_close}->{$n_far}++;
}

open(OUT, ">$ofn") or die "failed to open $ofn for writing\n";
print OUT "n_trans\tn_close\tn_far\tcount\n";

print STDERR "Writing output...\n";

foreach my $trans (keys %{$counts}) {		
	foreach my $close_cis (keys %{ $counts->{$trans} } ) {		
		foreach my $far_cis (keys %{ $counts->{$trans}->{$close_cis} } ) { 
			print OUT "$trans\t$close_cis\t$far_cis\t". $counts->{$trans}->{$close_cis}->{$far_cis} . "\n";
		}
	}
}
close(OUT)




