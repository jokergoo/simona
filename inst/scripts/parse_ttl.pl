use strict;

my $file = shift(@ARGV);
my @relation_types = @ARGV;
if($file =~/\.gz$/) {
	open FILE, "gzip -d -c $file |" or die "cannot open $file.";
} else {
	open FILE, $file or die "cannot open $file.";
}

my $line;
my $section = {};
my $id;
my $i_record = 0;
while(my $line = <FILE>) {
	if($line =~/^\@prefix /) {
		next;
	}

	if($line =~/^\s*$/) {
		next;
	} else {
		if($line =~/owl:Class/ and $line !~/\/STY\//) {
			$line =~/<(.*?)>/;
			$id = $1;
			$section->{$id} = {};
			$i_record ++;

			# if($i_record % 10000 == 0) {
			# 	print "$i_record finished...\n";
			# }

			while($line = <FILE>) {
				if($line =~/skos:prefLabel/) {
					$line =~/"""(.*?)"""/;
					$section->{$id}->{prefLabel} = $1;
					$section->{$id}->{prefLabel} =~s/"/``/g;
				}
				if($line =~/skos:notation/) {
					$line =~/"""(.*?)"""/;
					$section->{$id}->{notation} = $1;
				}
				if($line =~/skos:definition/) {
					$line =~/"""(.*?)"""/;
					$section->{$id}->{definition} = $1;
					$section->{$id}->{definition} =~s/"/``/g;
				}
				if($line =~/rdfs:subClassOf/ or $line =~/\/is_?a/i) {
					if($line =~/<([^<]+?)> ;/) {
						if(!defined($section->{$id}->{parent})) {
							$section->{$id}->{parent} = {};
							$section->{$id}->{parent}->{$1} = 1;
							$section->{$id}->{relation_type} = {};
							$section->{$id}->{relation_type}->{$1} = "is_a";

						} else {
							$section->{$id}->{parent}->{$1} = 1;
							$section->{$id}->{relation_type}->{$1} = "is_a";
						}
					}
				}
				foreach my $type (@relation_types) {
					if($line =~/\/$type/i) {
						if($line =~/<([^<]+?)> ;/) {
							if(!defined($section->{$id}->{parent})) {
								$section->{$id}->{parent} = {};
								$section->{$id}->{parent}->{$1} = 1;
								$section->{$id}->{relation_type} = {};
								$section->{$id}->{relation_type}->{$1} = $type;
							} else {
								$section->{$id}->{parent}->{$1} = 1;
								$section->{$id}->{relation_type}->{$1} = $type;
							}
						}
					}
				}

				if($line =~/\.$/) {
					last;
				}
			}
		} else {
			while($line = <FILE>) {
				if($line =~/\.$/) {
					last;
				}
			}
		}
	}
}

if($i_record == 0) {
	die "cannot find any object of 'owl:Class'.";
}

print "\"id\",\"prefLabel\",\"notation\",\"definition\",\"parent\",\"relation_type\"\n";

foreach $id (sort keys %$section) {
	print "\"$id\"";
	if(!defined($section->{$id}->{prefLabel})) {
		print ",\"\"";
	} else {
		print ",\"$section->{$id}->{prefLabel}\"";
	}
	if(!defined($section->{$id}->{notation})) {
		print ",\"\"";
	} else {
		print ",\"$section->{$id}->{notation}\"";
	}
	if(!defined($section->{$id}->{definition})) {
		print ",\"\"";
	} else {
		print ",\"$section->{$id}->{definition}\"";
	}
	if(!defined($section->{$id}->{parent})) {
		print ",\"\"";
	} else {
		print ",\"".join(",", keys %{$section->{$id}->{relation_type}})."\"";
	}
	if(!defined($section->{$id}->{relation_type})) {
		print ",\"\"";
	} else {
		print ",\"".join(",", values %{$section->{$id}->{relation_type}})."\"";
	}
	print "\n";
}

