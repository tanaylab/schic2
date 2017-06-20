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

	return($self);
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
	return(undef);
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
