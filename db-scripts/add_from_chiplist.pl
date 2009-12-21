#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $addpl = "java edu.mit.csail.cgs.tools.chipchip.AddDataFile ";


my ($expt,$flags,$version);
GetOptions("expt=s"=>\$expt,
	   "flags=s"=>\$flags,
	   "version=s"=>\$version);
warn "Didn't supply --expt, using default naming scheme" unless ($expt);
warn "No --flags supplied: no normalizations, version etc?" unless ($flags);
$flags ||= '';

while (<STDIN>) {
  chomp;
  next if (/^\s*\#/);
  my ($fname,$design,$species,$fragdistname,$fragdistversion,$cells1,$cells2,$cond1,$cond2,$fact1,$fact2,$replicate) = split("\t",$_);
  my $exptname = undef;
  if ($expt) { 
    $exptname =$expt;
  } else {
    my $spec = $species;
    if ($spec eq 'Saccharomyces cerevisiae') {
      $exptname = "Sc ${fact1}:${cells1}:${cond1} vs ${fact2}:${cells2}:${cond2}";
    } elsif ($spec eq 'Mus musculus') {
      $exptname = "Mm ${fact1}:${cells1}:${cond1} vs ${fact2}:${cells2}:${cond2}";
    } elsif ($spec eq 'Homo sapiens') {
      $exptname = "Hs ${fact1}:${cells1}:${cond1} vs ${fact2}:${cells2}:${cond2}";
    } elsif ($spec eq 'Drosophila melanogaster') {
      $exptname = "Dm ${fact1}:${cells1}:${cond1} vs ${fact2}:${cells2}:${cond2}";
    }
  }
  my $cmd = "$addpl --expt \"$exptname;$version;$replicate\" --designname \"$design\" --species \"$species\" ".
    "--fragdistname \"$fragdistname\" --fragdistversion \"$fragdistversion\" --cellsone \"$cells1\" ".
    "--conditionone \"$cond1\" --factorone \"$fact1\" --cellstwo \"$cells2\" --conditiontwo \"$cond2\" ".
    "--factortwo \"$fact2\" $flags --file \"$fname\"";
  print "$cmd\n";
  #  system("$cmd");
}
