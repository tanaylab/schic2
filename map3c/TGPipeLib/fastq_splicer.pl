use strict;

require map3c::TGPipeLib::FastqSplicer;
require map3c::ConfOpt;
require map3c::TGPipeLib::VFastq;

my($conf) = ConfOpt::new("ConfOpt");
$conf->parse_argv_conf(\@ARGV);

my($fq) = TGPipeLib::VFastq::new("TGPipeLib::VFastq");
my($fq2) = undef;
my($fq3) = undef;

$fq->set_scope($conf->get_opt("inp_base"), 
		$conf->get_opt("inp_re"), 
		$conf->get_opt("inp_r1_re","R1"));

if($fq->scope_size() == 0) {
	if($conf->get_opt("inp_base") !~/\/$/) {
		print STDERR "ERROR: note that directoy names must end with a \/ \n";
	}
	print STDERR "ERROR: no fastq input files\n";
}

if($conf->get_opt("r2", "NA") ne "NA"
|| $conf->get_opt("r3", "NA") ne "NA") {
	$fq2 = TGPipeLib::VFastq::new("TGPipeLib::VFastq");
	$fq2->set_scope($conf->get_opt("inp_base"), 
			$conf->get_opt("inp_re"), 
			$conf->get_opt("inp_r2_re","R2"));
}
if($conf->get_opt("r3", "NA") ne "NA") {
	$fq3 = TGPipeLib::VFastq::new("TGPipeLib::VFastq");
	$fq3->set_scope($conf->get_opt("inp_base"), 
			$conf->get_opt("inp_re"), 
			$conf->get_opt("inp_r3_re","R3"));
}

#$conf->dump_all();
print STDERR "will apply splicer\n";
TGPipeLib::FastqSplicer::proc_line(-1, $conf->{conf}, $fq, $fq2, $fq3);
