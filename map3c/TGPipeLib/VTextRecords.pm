package TGPipeLib::VTextRecords;

use IO::File;
use strict;

1;

#use diagnostics;

sub new($$$) {
	my($clas, $fread_rec, $fskip_rec) = @_;
	my($self) = {};
	
	bless($self,$clas);

	$self->{fread} = $fread_rec;
	$self->{fskip} = $fskip_rec;

	$self->{scope} = [];
	$self->{scope_start} = [];
	$self->{scope_end} = [];
	$self->{cur} = {};

	return($self);
}

sub scope_size {
	my($self) = @_;
	return(1+$#{$self->{scope}});
}

sub set_scope {
	my($self, $pref, $re, $re2) = @_;

	my(@fns) = <$pref*>;

	my($fn);
	if(defined($re)) {
		foreach $fn (@fns) {
			if($fn =~/$re/
			&& (!defined($re2) 
			  || $re2 eq "NA"
			  || $fn=~/$re2/)) {
				push(@{$self->{scope}}, $fn);
			}
		}
	} else {
		my(@fns) = <$pref*>;
		$self->{scope} = \@fns;
	}
	$self->set_default_subscope();
}

sub set_scope_by_list($$) {
	my($self, $fns) = @_;

	@{$self->{scope}} = @$fns;
	$self->set_default_subscope();
}

sub set_default_subscope($) {
	my($self) = @_;

	my($max_i) = $#{$self->{scope}};
	for(my($i) = 0; $i <= $max_i; $i++) {
		$self->{scope_start}->[$i] = 0;
		$self->{scope_end}->[$i] = 1e+18;
	}
}

sub gen_split_sscope($$$) {
	my($self, $nparts, $out_fn) = @_;


	open(SPLIT, ">$out_fn") || die "cannot write split scop at $out_fn\n";
	#we need to count
	my($scope) = $self->{scope};

	my($tot_l) = 0;
	$self->{lines} = [];
	for(my($i) = 0; $i <= $#$scope; $i++) {
		$self->update_line_count($i);
		$tot_l += $self->{lines}->[$i];
	}

	my($split_sz) = int($tot_l/$nparts)+1;

	my($split_i) = 0;
	my($split_cur_l) = 0;

	for(my($i) = 0; $i <= $#$scope; $i++) {
		my($file_max) = $self->{lines}->[$i];
		my($file_i) = 0;
		my($fn) = $self->{scope}->[$i];
		while($file_i < $file_max) {
			my($file_to_i) = $file_i + ($split_sz - $split_cur_l);
			if($file_to_i > $file_max) {
				$file_to_i = $file_max;
			}
			$split_i < $nparts || die "splitting into too many parts, nparts $nparts and now splitting into $split_i, file $fn, file_i $file_i file_max $file_max\n";
			print SPLIT "$fn\t$file_i\t$file_to_i\t$split_i\n";
			$split_cur_l += ($file_to_i - $file_i);
			if($split_cur_l >= $split_sz) {
				$split_i++;
				$split_cur_l = 0;
			}
			$file_i = $file_to_i;
		}
	}
}

sub set_subscope($$$) {
	my($self, $fn, $foc_ssid) = @_;

	open(SSTAB, $fn) || die "cannot open subscope table $fn\n";

	$self->{scope} = [];
	$self->{scope_start} = [];
	$self->{scope_end} = [];
	while(<SSTAB>) {
		chop;
		my($fn, $fr, $to, $ssid) = split("\t", $_);
		if($ssid == $foc_ssid) {
			push(@{$self->{scope}}, $fn);
			push(@{$self->{scope_start}}, $fr);
			push(@{$self->{scope_end}}, $to);
		}
	}
}

sub first($) {
	my($self) = @_;

	$self->{cur_i} = 0;
	$self->{cur_l} = 0;
	$self->{invalid} = undef;
	$self->init_to_scope(0);

	$self->{fread}->($self->{cur_io}, $self->{cur});
}


sub next($) {
	my($self) = @_;

	$self->{cur_l}++;
	my($cur_i) = $self->{cur_i};
	my($cur_l) = $self->{cur_l};
	my($io) = $self->{cur_io};
	if($io->eof() || $cur_l >= $self->{scope_end}->[$cur_i]) {
		$cur_i++;
		if($cur_i <= $#{$self->{scope}} ){
			$self->{cur_i}++;
			$self->init_to_scope($self->{cur_i});
		} else {
			$self->{invalid} = "OK";
		}
	} else {
		$self->{cur_l}++;
	}
	if(!defined($self->{invalid})) {
		my($status) = $self->{fread}->($self->{cur_io}, $self->{cur});
		if(defined($status)) {
			if($status eq "eof") {
				$self->next();	#valid eof - try going to next file
			} else {
				$self->{invalid} = $status;
			}
		}
	}
}


sub valid($) {
	my($self) = @_;
	return(!defined($self->{invalid}));
}
sub status($) {
	my($self) = @_;
	return($self->{invalid});
}

sub init_to_scope($$) {
	my($self, $i) = @_;

	print STDERR "init scope $i ".($self->{scope}->[$i])."\n";
	$self->{cur_io} = IO::File::new();
	$self->{cur_io}->open($self->{scope}->[$i]);
	if($self->{scope_start} != -1) {
		my($start) = $self->{scope_start}->[$i];
		my($valid)= 1;
		my($l) = 0;
		while($valid && $l < $start) {
			$valid = $self->{fskip}->($self->{cur_io});
			$l++;
		}
		if(!$valid) {
			$self->{invalid} = "Skip error";
		}
		$self->{cur_l} = $start;
	} else {
		$self->{cur_l} = 0;
	}
}

sub update_line_count($$) {
	my($self, $i) = @_;

	my($io) = IO::File::new();
	$io->open($self->{scope}->[$i]);
	my($lc) = 0;
	my($valid)= 1;
	while($valid) {
		$valid = $self->{fskip}->($io);
		$lc++;
	}
	$self->{lines}->[$i] = $lc;
}
