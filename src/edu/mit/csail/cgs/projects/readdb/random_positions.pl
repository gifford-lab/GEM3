#!/usr/bin/perl

use strict;
use warnings;

my $n = $ARGV[0];

while ($n-- > 0) {
  my $chr = int(rand() * 15) + 1;
  my $start = int(rand() * 200000);
  my $end = $start + 10000;
  print "chr${chr}:${start}-${end}:-\n";
}
