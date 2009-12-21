#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# takes a directory as the argument

my $addbayes = $0;
$addbayes =~ s/old_alan/bayesian/;
my ($species,$version,$step,$genomeVersion, $exptVersion);
GetOptions("species=s"=>\$species,
	   "genomeVersion=s"=>\$genomeVersion,
	   "step=s"=>\$step,
	   "exptVersion=s"=>\$exptVersion,
	   "version=s"=>\$version);

my $speciesprefix;
if ($species eq 'Homo sapiens') {
  $speciesprefix = 'Hs';
} elsif ($species eq 'Mus musculus') {
  $speciesprefix = 'Mm';
} elsif ($species eq 'Saccharomyces cerevisiae') {
  $speciesprefix = 'Sc';
} else {
  die "Unknown species for prefix $species";
}

my %expts;
foreach my $dir (@ARGV) {
  opendir(DIR,$dir) or die "can't open dir $dir : $!";
  foreach (readdir(DIR)) {
    next unless ($_ =~ /posterior.txt$/);
    if ($_ =~ /^\w+\.([\w\s]*)\./) {
      my $expt = $1;
      my $prefix = "${dir}/${&}";
      my ($f1,$f2,$cond) = ($expt =~ /(.*)\svs\s(.*)\sin\s(.*)/);
      $cond =~ s/hES.h9/unknown/;
      $expts{$prefix} = {dir=>$dir,
			 pattern=>"$expt\\.",
			 expt=>"$speciesprefix $f1 vs $f2 in $cond"};
    } else {
      print STDERR "$$_ didn't match\n";
    }
  }
}
foreach (keys %expts) {
  print "$addbayes --step $step --species '$species' --genomeVersion '$genomeVersion' --exptVersion '$exptVersion' --version '$version' " .
    "--analysis '$expts{$_}->{expt}' --expt '$expts{$_}->{expt}' --dir '$expts{$_}->{dir}' --pattern '$expts{$_}->{pattern}'\n";
}
