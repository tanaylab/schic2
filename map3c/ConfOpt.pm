use strict;

package ConfOpt;

sub new {
	my($clas, $legit) = @_;


	my($self) = {};
	
	bless($self,$clas);

	$self->{conf} = {};

	if(defined($legit)) {
		$self->{all_opts} = {};
		%{$self->{all_opts}} = %{$legit};
		$self->{force_space} = 1;
	} else {
		$self->{force_space} = 0;
	}

	return($self);
}

sub load_legit($) {


}

sub get_opt {
	my($self, $nm, $default) = @_;
	if(!exists($self->{conf}->{$nm})) {
		if(defined($default)) {
			return($default);
		}
		die "cannot find required config opt $nm\n";
	}
	return($self->{conf}->{$nm});
}

sub read_conf($$)  {
	my($self, $fn) = @_;
	local *CONF;
	open(CONF, $fn) || die "cannot open opt conf file $fn\n";
	while(<CONF>) {

		# include files
		if($_=~/^#include\s*(\S+)/) {
			if (!-e $1) { die "cannot find file $1 to include\n"; }
			$self->read_conf($1);
			next;
		}

		# skipping comment line
		if($_=~/^#/) {
			next;
		}

		my($nm, $val) = /(^[^=]+)\s?=\s?(.*)$/;
		if(defined($nm)) {
			if($self->{force_scope} 
			&& !exists($self->{all_opts}->{$nm})) {
				die "conf file with invalid opt name $nm\n";
			}  
			if ($val =~ /\$ENV\{\"(\w+)\"\}(.+)/) {
				$val = $ENV{$1} . $2;
			}
			if ($val =~ /(.*)\$\{(\S+)\}(.*)/) {
				$val = $1 . $self->{conf}->{$2} . $3;
			}
			$self->{conf}->{$nm} = $val;
		}
	}
}

sub parse_argv_conf {
	my($self, $params) = @_;

	my($i);

	for($i = 0; $i <= $#$params; $i++) {
		if($params->[$i] =~/^\@/) {
			my($opt_fn) = substr($params->[$i],1);
			$self->read_conf($opt_fn);
			next;
		}
		if($params->[$i] !~/^\-/) {
			die "bad parameters format (should be -x x_val -y y_val) at $i\n";
		}
		my($nm) = substr($params->[$i], 1);
		if($self->{force_scope}
		&& !exists($self->{all_opts}->{$nm})) {
			die "unrecogznied argv option(s) $nm\n";
		}
		if($nm =~/=/) {
			die "Bad option name $nm - cannot use \"=\"\n";
		}
		if($i == $#$params || $params->[$i+1] =~/^\-/) {
			$self->{conf}->{$nm} = 1;
		} else {
			$self->{conf}->{$nm} = $params->[$i+1];
			$i++;
		}
	}
}
sub dump_all {
	my($self) = @_;
	my($k);
	foreach $k (keys %{$self->{conf}}) {
		print STDOUT "$k\t".($self->{conf}->{$k})."\n";
	}

}
1;
