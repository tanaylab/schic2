use strict;

package TGPipeLib::FastqSplicer;

1;


sub required() {
	return("re1", "R1", "out_fn1");	#can also  have re2 and re3
}

sub proc_line
{
	my($split_id, $params, $fastq, $fastq2, $fastq3) = @_;

	my($verbose) = 0;

	my($n) = 0;
	my($bad);
	my($code) = "\$bad = 0;\n";
	if(exists($params->{"re1"})) {
#		print STDERR "got paras for re1, ".($params->{"re1"})."\n";
		defined($fastq) || die "Cannot run code with regexp on read1 without fastq 1\n";
		$code .= 
		"\$R1=\$fastq->cur_seq();\n".
		"\$Q1=\$fastq->cur_qual();\n";
		if($params->{"re1"} ne "") {
		   $code .= "if(".($params->{"re1"}).") {\n";
		   $code .="for(my(\$h)=1; \$h <= \$#+; \$h++) {\n".
		     "     push(\@S, substr(\$R1, \$-[\$h], \$+[\$h]-\$-[\$h]));\n".
		     "     push(\@U, substr(\$Q1, \$-[\$h], \$+[\$h]-\$-[\$h]));\n }";
		   $code .= "} else { \$bad=1; }";
		}
	} else {
		$code .= "\n\$R1=\$fastq->cur_seq();\n";
		$code .= "\$Q1=\$fastq->cur_qual()\n;";
	}
	if(exists($params->{"re2"})) {
		defined($fastq2) || die "Cannot run code with regexp on read2 without fastq 2\n";
		$code .= 
		"\$R2=\$fastq2->cur_seq();".
		"\$Q2=\$fastq2->cur_qual();";
		if($params->{"re2"} ne "") {
		   $code .= "if(".($params->{"re2"}).") {\n";
		   $code .="for(my(\$h)=1; \$h <= \$#+; \$h++) {\n".
		     "     push(\@S, substr(\$R2, \$-[\$h], \$+[\$h]-\$-[\$h]));\n".
		     "     push(\@U, substr(\$Q2, \$-[\$h], \$+[\$h]-\$-[\$h])); \n}\n";
		   $code .= "} else { \$bad=1; }\n";
		}
	} elsif(defined($fastq2)) {
		$code .= "\n\$R2=\$fastq2->cur_seq();\n";
		$code .= "\$Q2=\$fastq2->cur_qual()\n;";
	} else {
		print STDERR "fastq2 is undefined??\n";
	}

	if(exists($params->{"re3"})) {
		defined($fastq3) || die "Cannot run code with regexp on read3 without fastq 3\n";
		$code .= 
		"\$R3=\$fastq3->cur_seq();".
		"\$Q3=\$fastq3->cur_qual();";
		if($params->{"re3"} ne "") {
		   $code .= "if(".($params->{"re3"}).") {\n";
		   $code .="for(my(\$h)=1; \$h <= \$#+; \$h++) {".
		     "     push(\@S, substr(\$R3, \$-[\$h], \$+[\$h]-\$-[\$h]));".
		     "     push(\@U, substr(\$Q3, \$-[\$h], \$+[\$h]-\$-[\$h])); }";
		   $code .= "} else { \$bad=1; }";
		}
	} elsif(defined($fastq3)) {
		$code .= "\n\$R3=\$fastq3->cur_seq();";
		$code .= "\$Q3=\$fastq3->cur_qual();";
	}
	my($o1, $o2, $o3);

	if(exists($params->{"ID"})) {
		if($params->{"ID"} =~/^\$ID=/) {
			die "ID line should start with \$ID=";
		}
		$code .= $params->{"ID"};
	}
	if(!exists($params->{"O1"})) {
		die "missing options O1, splicer should write at least one file\n";
	}
	$code .= "\nif(!\$bad) {\n";
	if($params->{"O1"} !~/^\$O1=/) {
		die "O1 line should start with \$O1=";
	}
	if(!exists($params->{"out_fn1"})) {
		die "missing out_fn1 parameter\n";
	}
	my($O1) = $params->{"O1"};
	my($O2);
	my($O3);
	my($OQ1) = $O1;
	$OQ1 =~ s/R/Q/g;
	$OQ1 =~ s/S/U/g;
	$OQ1 =~ s/O1=/Q1=/g;
	my($fn) = $params->{"out_fn1"};
	$code .= "\n$O1;\n";
	$code .= "$OQ1;\n";
	$code .= "\$o1->print(\"\$ID\\n\");";
	$code .= "\$o1->print(\"\$O1\\n\");";
	$code .= "\$o1->print(\"+\\n\");";
	$code .= "\$o1->print(\"\$Q1\\n\");";
	$o1 = IO::File::new();
	$o1->open(">$fn") || die "cannot open output file 1\n";

	if(exists($params->{"O2"})) {
		if($params->{"O2"} !~/^\$O2=/) {
			die "O2 line should start with \$O2=";
		}
		if(!exists($params->{"out_fn2"})) {
			die "missing out_fn2 parameter\n";
		}
		$O2 = $params->{"O2"};
		my($OQ2) = $O2;
		$OQ2 =~ s/R/Q/g;
		$OQ2 =~ s/S/U/g;
		$OQ2 =~ s/O2=/Q2=/g;
		my($fn) = $params->{"out_fn2"};
		$code .= "\n$O2;\n";
		$code .= "$OQ2;\n";
		$code .= "\$o2->print(\"\$ID\\n\");";
		$code .= "\$o2->print(\"\$O2\\n\");";
		$code .= "\$o2->print(\"+\\n\");";
		$code .= "\$o2->print(\"\$Q2\\n\");";
		$o2 = IO::File::new();
		$o2->open(">$fn") || die "cannot open output file 2\n";
	}
	if(exists($params->{"O3"})) {
		$O2 = $params->{"O3"};
		my($OQ3) = $O3;
		$OQ3 =~ s/R/Q/g;
		$OQ3 =~ s/S/U/g;
		$OQ3 =~ s/O2=/Q3=/g;
		my($fn) = $params->{"out_fn3"};
		$code .= "\n$O3;\n";
		$code .= "$OQ3;\n";
		$code .= "\$o3->print(\"\$ID\\n\");";
		$code .= "\$o3->print(\"\$O3\\n\");";
		$code .= "\$o3->print(\"+\\n\");";
		$code .= "\$o3->print(\"\$Q3\\n\");";
		$o3 = IO::File::new();
		$o3->open(">$fn") || die "cannot open output file 3\n";
	}
	$code .= "\n}";

	if($verbose) { print STDERR "code is $code\n\n"; }

	if($code eq "") {
		die "FastqSplicer: no re1, re2, re3, O1, O2 or O3 options, so nothing to splice!\n";
	}

	my($OUT1, $OUT2, $OUT3);

	no strict "vars";
	my(%count);
	$fastq->first();
	if(!defined($fastq2)) {
		my($wrap) = 
		"while(\$fastq->valid()) {\n ".
			"my(\@S, \@U);\n ".
			"\$fastq->next(); \n".
			"my(\$ID) = \$fastq->cur_id(); \n".
			"$code; \n".
		"}";
		eval $wrap;
		die "Gen code failed $@\n" if $@;
	} else {
		$fastq2->first();
		if(!defined($fastq3)) {
			my($wrap) =
			"while(\$fastq->valid()) {\n ".
				"my(\@S, \@U); \n".
				"\$fastq->next(); \n".
				"\$fastq2->next(); \n".
				"my(\$ID) = \$fastq->cur_id(); \n".
				"$code;\n ".
			"}";
			eval $wrap;
			die "Gen code failed $@\n" if $@;
		} else {
			$fastq3->first();
			my($wrap) =
			"while(\$fastq->valid()) {\n ".
				"my(\@S, \@U); \n".
				"\$fastq->next(); \n".
				"\$fastq2->next(); \n".
				"\$fastq3->next(); \n".
				"my(\$ID) = \$fastq->cur_id(); \n".
				"$code;\n ".
			"}";
			eval $wrap;
			die "Gen code failed $@\n" if $@;
		}
	}
	use strict "vars";
}

sub reduce($) {
	my($params) = @_;

	print STDERR "nothing to reduce here\n";
}
