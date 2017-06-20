use strict;

require map3c::ConfOpt;
require map3c::TGPipeLib::DistribLineProc;

if($#ARGV < 4) {
	print STDERR "usage DLP.pl proc_lib file_lib N work_dir scope_base [\@conf_file|-a b|]\n";
	print STDERR "example: perl DLP.pl TGPipeLib::FastqRegexp TGPipeLib::VFastq 10 work/ results/mydata/ \@safta.conf";
	print STDERR "VFastq can be abbreviated fq\n";
	print STDERR "TGPipeLib prefix can be ommited\n";
	die "";
}


my($proc_lib) = $ARGV[0];
if($proc_lib =~/^::/) {
	$proc_lib = "TGPipeLib".$proc_lib;
}
my($file_lib) = $ARGV[1];
if($file_lib =~/^::/) {
	$file_lib = "TGPipeLib".$proc_lib;
} elsif($file_lib eq "fq") {
	$file_lib = "TGPipeLib::VFastq";
}

eval "use $file_lib";
if($@ != "") {
	die "$@\nparse error in the code of $file_lib\n";
}

my($N) = $ARGV[2];
my($work_dir) = $ARGV[3];
if($work_dir !~/^\//) {
	$work_dir = $ENV{"PWD"}."/$work_dir";
}
my($scope_base) = $ARGV[4];
if($scope_base !~/^\//) {
	$scope_base= $ENV{"PWD"}."/$scope_base";
}

my($conf) = ConfOpt::new("ConfOpt");
@ARGV = @ARGV[5..$#ARGV];
$conf->parse_argv_conf(\@ARGV);

my($code_dir) = $conf->get_opt("code_dir", $ENV{"PWD"}."/lscripts");
if(!-e "$code_dir/TGPipeLib/DistribLineProcJob.pl") {
	die "code_dir was set to $code_dir, but DLP Job script (TGPipeLib/DistribLineProcJob.pl) is not present there\n";
}

my($lines1) = &{\&{$file_lib."::new"}}($file_lib);
$lines1->set_scope($scope_base, $conf->get_opt("scope1"), $conf->get_opt("scope1_re2", "NA"));
my($lines2);
my($scope2) = $conf->get_opt("scope2", "NA");
if($scope2 ne "NA") {
	$lines2 = &{\&{$file_lib."::new"}}($file_lib);
	$lines2->set_scope($scope_base, $scope2, $conf->get_opt("scope2_re2", "NA"));
}

my($dlp) = TGPipeLib::DistribLineProc::new(
	    "TGPipeLib::DistribLineProc",
	    $N, 
	    $code_dir, 
	    $work_dir, 
	    $proc_lib,
	    $file_lib,
	    $lines1, $lines2, 
	    $conf->{conf});

$dlp->run_on_qsub();
