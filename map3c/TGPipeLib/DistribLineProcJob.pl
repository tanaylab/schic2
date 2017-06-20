use strict;

require map3c::ConfOpt;

if($#ARGV != 6) {
	die "usage main work_dir proc_lib file_lib split1 split2|NA conf_file split_i\n";
}

my($work_dir) = $ARGV[0];
my($proc_lib) = $ARGV[1];
my($file_lib) = $ARGV[2];
my($split1) = $ARGV[3];
my($split2) = $ARGV[4];
my($conf_fn) = $ARGV[5];

my($conf) = ConfOpt::new("ConfOpt");
$conf->read_conf($conf_fn);

my($split_i) = $ARGV[6];

eval "use $proc_lib";
eval "use $file_lib";
#gene a new VTextRecord
#use the lib;
#
my($r1) = &{\&{$file_lib."::new"}}($file_lib);
$r1->set_subscope($split1, $split_i);
my($r2) = "NA";
if($split2 ne "NA") {
	$r2 = &{\&{$file_lib."::new"}}($file_lib);
	$r2->set_subscope($split2, $split2);
}


&{\&{$proc_lib."::proc_line"}}($split_i, $conf->{conf}, $r1, $r2);

if(!-e "$work_dir/status") {
	mkdir("$work_dir/status");
}
open(TAG, ">$work_dir/status/$split_i") || die "cannot write status file!\n";
print TAG "OK";
close TAG;
