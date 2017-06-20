package TGPipeLib::DistribLineProc;
use strict;

1;

sub new {
	my($clas, 
	    $N, 
	    $code_dir, 
	    $base_dir, 
	    $proc_lib,
	    $file_lib,
	    $lines1, $lines2, 
	    $params) = @_;

	my($self) = {};
	
	bless($self,$clas);

	my($runi) = 0;
	while(-e "$base_dir/work$runi") {
		$runi++;
	}
	my($work_dir) .= "$base_dir/work$runi";
	$self->{work_dir} = $work_dir;
	$self->{code_dir} = $code_dir;
	mkdir($work_dir);
	$params->{work_dir} = $work_dir;

	$self->{split1_tab_fn} = "$work_dir/split_tab1";
	$self->{split2_tab_fn} = "$work_dir/split_tab2";

	$self->{N} = $N;
	$self->{proc_lib} = $proc_lib;
	$self->{file_lib} = $file_lib;

	$self->{lines1} = $lines1;
	$lines1->gen_split_sscope($N, $self->{split1_tab_fn});
	if(defined($lines2)) {
		$self->{lines2} = $lines2;
		$lines2->gen_split_sscope($N, $self->{split2_tab_fn});
	} else {
		$self->{split2_tab_fn} = "NA";
	}

	$self->{params} = $params;
	$self->{max_restarts} = 4;
	if(exists($params->{max_restarts})) {
		$self->{max_restarts} = $params->{max_restarts};
	}
	$self->{max_time_gap} = 60;
	if(exists($params->{max_time_gap})) {
		$self->{max_time_gap} = $params->{max_time_gap};
	}
	$self->{max_slots} = 1000;
	if(exists($params->{max_slots})) {
		$self->{max_slots} = $params->{max_slots};
	}
	my($queue_list) = "all.q";

	$self->test_params();
	return($self);
}

sub run_on_qsub($) {
	my($self) = @_;

	eval "use ".($self->{proc_lib});

	$self->{time} = 0;
	$self->prepare_qscripts();

	$self->{pending} = [];
	$self->{pending_restart} = []; 
	$self->{running} = [];
	$self->{done} = [];

	$self->{job_id} = [];
	$self->{exit_stat} = [];	
	$self->{job_restarts} = [];
	$self->{id_to_job_i} = {};
	$self->{dead_on_queue} = {};

	for(my($i) = 0; $i < $self->{N}; $i++) {
		$self->{job_restarts}->[$i] = 0;	
		$self->{job_id}->[$i] = -1;
		push(@{$self->{pending}}, $i);
	}

	while(($#{$self->{running}} + $#{$self->{pending}}) >= 0) {
		$self->track_jobs();
		sleep(3);
		$self->{time}++;
	}
	$self->wrap_up();
	$self->clean();
}

sub prepare_qscripts($) {
	my($self) = @_;

	my($work_dir) = $self->{work_dir};
	my($code_dir) = $self->{code_dir};
	my($proc_lib) = $self->{proc_lib};
	my($file_lib) = $self->{file_lib};

	open(CONF, ">$work_dir/proc.conf") || die "Cannot open proc.conf file for writing\n";

	my($k);
	foreach $k (keys %{$self->{params}}) {
		print CONF "$k=".($self->{params}->{$k})."\n";
	}
	close CONF;

#main.pl - include relevantm modules and the specific proc, reduce functions
#
	my($proc_nm) = $proc_lib;
	my($split_tab2) = $self->{split2_tab_fn};
	$proc_nm =~s/::/./g;
	mkdir("$work_dir/qscripts");
	mkdir("$work_dir/qlogs");
	for(my($i) = 0; $i < $self->{N}; $i++) {

		open(SCR, ">$work_dir/qscripts/do$i.sh");
		print SCR "#!/bin/sh\n";
		print SCR "#\$ -N DLP$i.$proc_nm\n";
		print SCR "#\$ -e $work_dir/qlogs/err$i\n";
		print SCR "#\$ -o $work_dir/qlogs/out$i\n";
		print SCR "#\$ -S /bin/sh\n";
		print SCR "cd $work_dir\n";
		print SCR "export PERL5LIB=$code_dir\n";
		print SCR "perl $code_dir/TGPipeLib/DistribLineProcJob.pl $work_dir $proc_lib $file_lib $work_dir/split_tab1 $split_tab2 $work_dir/proc.conf $i\n";
	}
}

sub track_jobs($) {
	my($self) = @_;

	my($max_slots) = $self->{max_slots};
	my($max_restarts) = $self->{max_restarts};
	#extract status from qstat
	$self->query_q();	

	my($running) = $self->{running};
	for(my($i) = 0; $i <= $#$running; $i++) {
		my($job_i) = $running->[$i];
		my($stat) = $self->update_job_status($job_i);
		if($stat eq "done") {
			push(@{$self->{done}}, $job_i);
			$self->{exit_stat}->[$job_i] = "ok";
			print STDERR "job $job_i finished OK\n";
			splice(@$running, $i, 1);
			$i--;
		} elsif($stat eq "dead") {
			if($self->{job_restarts}->[$job_i] < $max_restarts) {
				print STDERR "job $job_i dead - try again\n";
				$self->{job_restarts}->[$job_i]++;
				push(@{$self->{pending_restart}}, $job_i);
			} else {
				print STDERR "job $job_i dead - giving up\n";
				push(@{$self->{done}}, $job_i);
				$self->{exit_stat}->[$job_i] = "failed";
			}
			splice(@$running, $i, 1);
			$i--;
		}
	}

	while($#{$self->{pending_restart}} != -1
	   && $#{$self->{running}} < $max_slots) {
		my($job_i) = shift @{$self->{pending_restart}};
		$self->launch_qsub($job_i);
	}
	while($#{$self->{pending}} != -1 
	   && $#{$self->{running}} < $max_slots) {
		my($job_i) = shift @{$self->{pending}};
		$self->launch_qsub($job_i);
	}
}

sub launch_qsub {
	my($self, $job_i) = @_;


	my($work_dir) = $self->{work_dir};
	my($scr) = "$work_dir/qscripts/do$job_i.sh";
	my($queue_list) = $self->{queue_list};
	system("qsub -terse $scr -q $queue_list  >$work_dir/launch.out");
	open(LAUNCH, "$work_dir/launch.out") || die "cannot open launch.out\n";
	my($id);
	$id = <LAUNCH>;
	chop $id;
	if($id !~/\d+/) {
		die "cannot submit job $job_i and get an id via qsub, id was $id\n";
	}
	push(@{$self->{running}}, $job_i);
	$self->{job_id}->[$job_i] = $id;
	$self->{id_to_job_i}->{$id} = $job_i;
	$self->{job_init_time}->[$job_i] = $self->{time};
	print STDERR "launched job $job_i to $id\n";
}

sub query_q {
	my($self) = @_;

	my($work_dir) = $self->{work_dir};
	system("qstat >$work_dir/log.qstat.out");

	open(STAT, "$work_dir/log.qstat.out") || die "cannot open log.qstat.out\n";
	$self->{qstat} = [];
	<STAT>; #header
	<STAT>; #seperator
	while(<STAT>) {	
		chop;
		my(@f) = split(/\s+/, $_);
#		print STDERR "parse qstat line $_, f[0] $f[0]\n";
		my($id) = $f[0];
		my($status) = $f[4];
		my($queue) = $f[7];
		my($job_i) = $self->{id_to_job_i}->{$id};
		if(defined($job_i)) {
			$self->{qstat}->[$job_i] = $status;
			$self->{qstat_queue}->[$job_i] = $queue;
		}
	}
}

sub update_job_status {
	my($self, $job_i) = @_;

	my($id) = $self->{job_id}->[$job_i];
	my($work_dir) = $self->{work_dir};

	my($status) = $self->{qstat}->[$job_i];
	my($queue) = $self->{qstat_queue}->[$job_i];
	my($max_time_gap) = $self->{max_time_gap};

#	print STDERR "job status $job_i id $id is $status queue $queue\n";
	if(!defined($status)) {
		#check if finished ok
		if(-e "$work_dir/status/$job_i") {
			return("done");
		} else {
			my($q) = $self->{qstat_queue}->[$job_i];
			$self->mark_fail($q);
			return("dead");
		}
	} elsif($status eq "r") {
		$self->{job_queue}->[$job_i] = $queue;
		return("running\n");
		#check if expired? 
	} else {
		if($self->{time} - $self->{job_init_time}->[$job_i] > $max_time_gap) {
#			print STDERR "job $job_i id $id is too long in non run state - we kill an reinit\n";
			#too much time without going into run state	
			system("qdel -j $id");
			return("dead");
		}
	}

}

sub mark_fail($$) {
	my($self, $queue) = @_;

	$self->{dead_on_queue}->{$queue}++;
	if($self->{dead_on_queue}->{$queue} > $self->{max_queue_fail}) {
		#black list the queue?
	}
}

sub wrap_up($) {
	my($self) = @_;

	my($proc_lib) = $self->{proc_lib};
	&{\&{$proc_lib."::reduce"}}($self->{params});
}


sub run_on_local($) {
	my($self) = @_;

	die "Not implemented yet\n";

}

sub test_params($) {
	my($self) = @_;

	my($proc_lib) = $self->{proc_lib};
	my($params) = $self->{params};

	eval "use ".($self->{proc_lib});
	my(@required) = eval "$proc_lib"."::required";

	my($miss) = 0;
	my($s);
	foreach $s (@required) {
		if(!exists($params->{$s})) {
			print STDERR "missinig param $s in FastqCountRegexp\n";
			$miss++;
		}
	}
	if($miss > 0) {
		die "cannot operate without these params\n";
	}
}

sub clean($) {
	my($self) = @_;

	if(!exists($self->{params}->{debug})) {
		my($work_dir) = $self->{work_dir};
		my($code_dir) = $self->{code_dir};
		#defer this a bit because nfs handles takes time to relax\n";
		system("perl $code_dir/delayed_rm.pl $work_dir 30 &");
	}
}
