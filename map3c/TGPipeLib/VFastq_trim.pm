package TGPipeLib::VFastq;

use strict;

#use diagnostics;
#
#
1;

print "Using vfastq\n";

use parent 'TGPipeLib::VTextRecords';

sub new($) {
	my($clas) = @_;
	
	my($self) = TGPipeLib::VTextRecords::new($clas, \&fastq_read_rec, \&fastq_skip_rec);

	$self->{trim} = -1;
	$self->{min_phred} = -1;
	return($self);
}

sub set_trim($$) {
	my($self, $trim) = @_;
	$self->{trim} = $trim;
}

sub set_min_phred($$) {
	my($self, $phred) = @_;
	$self->{min_phred} = $phred;
}

sub fastq_read_rec($$) {
	my($io, $rec) = @_;

	if($io->eof()) {
		return("eof");
	}
	$rec->{id} = <$io>;
	chop $rec->{id};
	$rec->{seq} = <$io>;
	chop $rec->{seq};
	<$io>;
	if($io->eof()) {
		return("premature eof");
	}
	$rec->{qual} = <$io>;
	chop $rec->{qual};
	if($rec->{trim} != -1) {
		my($trim) = $self->{trim};
		$rec->{seq} = substr($rec->{seq}, 0, $trim);
		$rec->{qual} = substr($rec->{qual}, 0, $trim);
	}
	if($rec->{min_phred} != -1 
	&& $rec->phred_below_thresh($rec->{qual}, $rec->{min_phred}) {
		return($rec->fastq_read_req($rec));
	} else {
		return(undef);
	}
}
sub cur_min_phred($$) {
	my($qual, $thresh) = @_;

	for(my($i) = 0; $i < length($qual); $i++) {
		if(ord(substr($qual, $i, 1)) < $thresh) {
			return(1);
		}
	}
	return(0);
}
sub fastq_skip_rec($$) {
	my($io) = @_;
#there may be a more efficient way
	<$io>;
	<$io>;
	<$io>;
	<$io>;
}

sub cur_id($) {
	my($self) = @_;

	return($self->{cur}->{id});
}
sub cur_seq($) {
	my($self) = @_;

	return($self->{cur}->{seq});
}
sub cur_qual($) {
	my($self) = @_;

	return($self->{cur}->{qual});
}
