#!/usr/bin/perl

use strict;
use warnings;

my ($baseversion) = shift(@ARGV);
print STDERR "Using version base $baseversion\n";

foreach my $dir (@ARGV) {
  open(ARGS,"${dir}/args.txt") or (warn "Can't open ${dir}/args.txt : $!",next);
  my $args = <ARGS>;
  close ARGS;
  my $exptnameversion;
  if ($args =~ /expt\s+[\'\"]([^\'\"]+)[\'\"]\s/) {
    $exptnameversion = $1;
  } else {
    warn "Can't get exptnameversion from args.txt in $dir : $args.";
    next;
  }
  my ($expt,$version,$rep) = split(';',$exptnameversion);
  my $thisversion = '';
  if ($dir =~ /bijective/) {
    $thisversion = "${baseversion} bijective";
  } elsif ($dir =~ /tenarray/) {
    $thisversion = "${baseversion} tenarray";
  } elsif ($dir =~ /all/) {
    $thisversion = "${baseversion} all";
  } else {
    warn "Can't figure out a version for $dir";
    $thisversion = $baseversion;
  }
  print "cd ${dir} && cat args.txt | xargs ~/projects/database/add_bayesian.pl --chromindex 1 --analysis '${expt};${thisversion}'\n";
}
