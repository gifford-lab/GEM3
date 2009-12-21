#!/usr/bin/perl

use strict;
use warnings;
my @binOffsets = (512+64+8+1,
		  64+8+1,
		  8+1,
		  1,
		  0);
my $binFirstShift = 17;
my $binNextShift = 3;

# start and stop are the 0-based column indices of the
# start pos and stop pos

my ($startcol,$stopcol);
use Getopt::Long;
GetOptions("start=s"=>\$startcol,
	   "stop=s"=>\$stopcol);

while (<STDIN>) {
  chomp;
  my @l = split(/\t/,$_);
  my $bin = bin_from_range($l[$startcol],$l[$stopcol]);
  print "${bin}\t$_\n";
}

sub bin_from_range {
  my ($start,$end) = @_;
  my $startbin = $start;
  my $endbin = $end - 1;
  my $i;
  $startbin >>= $binFirstShift;
  $endbin >>= $binFirstShift;
  for ($i = 0; $i < @binOffsets; ++$i) {
    if ($startbin == $endbin) {
      return $binOffsets[$i] + $startbin;     
    }
    $startbin >>= $binNextShift;
    $endbin >>= $binNextShift;
  }
  warn "Can't get a bin for $start, $end : $startbin, $endbin";
  return 0;
}
