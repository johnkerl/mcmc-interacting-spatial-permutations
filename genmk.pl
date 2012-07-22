#!/usr/bin/perl -Wall

# ================================================================
# John Kerl
# kerl.john.r@gmail.com
# 2008-03-27
# ================================================================

@mks = ();
while ($line=<>) {
	chomp $line;

	# Strip out comments.
	$line =~ s/#.*//;

	# Strip out whitespace.
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;

	# Skip blank/comment lines.
	next if ($line =~ m/^\s*$/);

	if (! ($line =~ /.*\.mki$/)) {
		die "$0:  Unexpected line \"$line\" in input; expected .mki filenames one per line.\n";
	}
	$line =~ s/mki$/mk/;
	push @mks, $line;
}

print "# This file was automatically generated from Makefile.in by genmk.";
print "";

print "opt:";
for my $mk (@mks) {
	print "\texport OPTCFLAGS=\"-O3\" OPTLFLAGS=\"\"; make -ef $mk";
}
print "";

print "profile:";
for my $mk (@mks) {
	print "\texport OPTCFLAGS=\"-g -O3 -pg\" OPTLFLAGS=\"-g -pg\"; make -ef $mk";
}
print "";

print "gcov:";
for my $mk (@mks) {
	print "\texport OPTCFLAGS=\"-g -fprofile-arcs -ftest-coverage\" OPTLFLAGS=\"-g -fprofile-arcs -ftest-coverage\"; make -ef $mk";
}
print "";

print "debug:";
for my $mk (@mks) {
	print "\texport OPTCFLAGS=\"-g\" OPTLFLAGS=\"-g\"; make -ef $mk";
}
print "";

print "build:";
for my $mk (@mks) {
	print "\tmake -f $mk";
}
print "";

print "mk:";
for my $mk (@mks) {
	print "\tyamm ${mk}i";
}
print "";

print "install:";
for my $mk (@mks) {
	printf "\tmake -f %-30s install\n", $mk;
}
print "";

print "clean:";
for my $mk (@mks) {
	printf "\tmake -f %-30s clean\n", $mk;
}
print "";

print "tags: .PHONY";
print "\tctags *.[ch]";
print "";
print ".PHONY:";
print "";
print "over: clean mk build";
