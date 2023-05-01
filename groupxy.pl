#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (s=>.7, c=>.3, r=>.9);
getopts('s:c:r:', \%opts);
die("Usage: sexchr.pl [-s $opts{s}] [-c $opts{c}] [-r $opts{r}] in.sexchr\n") if @ARGV == 0;

# read data
my @a = ();
while (<>) {
	chomp;
	my @t = split("\t");
	next if $t[0] ne 'S';
	push(@a, \@t);
}

# assign sex chromosome for individual contigs
my @c = (0, 0, 0, 0);
for (@a) {
	next if $_->[5] < $_->[4] * $opts{s};
	next if $_->[6] + $_->[7] < $_->[5] * $opts{c};
	$_->[3] = $_->[6] > ($_->[6] + $_->[7]) * $opts{r}? 3 : $_->[7] > ($_->[6] + $_->[7]) * $opts{r}? 4 : 0;
	next if $_->[3] == 0;
	my $hap = $_->[2] - 1;
	$c[$hap<<1|0] += $_->[6];
	$c[$hap<<1|1] += $_->[7];
}

# determine which partition correspond to sexchr1 or sexchr2
my $max_chr = $c[0] + $c[2] > $c[1] + $c[3]? 0 : 1;
my $type = ($c[0<<1|$max_chr] > $c[1<<1|$max_chr]? 0 : 1) ^ $max_chr;

# update partition
for (@a) {
	if ($_->[3] >= 3) {
		$_->[3] -= 2;
	} else {
		$_->[3] = $type == 0? $_->[2] : 3 - $_->[2];
	}
}

# print out
for (@a) {
	print join("\t", @$_), "\n";
}
