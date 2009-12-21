#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# should be run on the .1.1.ll files from the desired analyses
# generates shell commands to call add_ll.pl to actually load
# the files

my ($species,$version);
GetOptions("species=s"=>\$species,
	   "version=s"=>\$version);

my $addll = $0;
$addll =~ s/_old//;

foreach my $base (@ARGV) {
  $base =~ s/\.1\.1\.ll$//;
  $base =~ s/^(histone\.)|(tf\.)|(pol2\.)//;
  $base = "Sc $base";
  my $expt = $base;
  my $analysis = $base;
  my $files = "$base.*.ll";
  $files =~ s/ /\\ /g;
  print "$addll --species '$species' --version '$version' --analysis '$analysis' --expt '$expt' $files\n";
}
