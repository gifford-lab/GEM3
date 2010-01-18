#!/usr/bin/perl


use strict;
use warnings;
use ReadDBClient;

my $client = new ReadDBClient();
my $cmd = shift(@ARGV);

if ($cmd eq 'getChroms') {
  my $align = shift(@ARGV);
  my @c = $client->getChroms($align);
  print "chroms for $align are @c\n";
} elsif ($cmd eq 'getCount') {
  my $align = shift(@ARGV);
  my $c = $client->getCount($align);
  print "count for $align is $c\n";
} elsif ($cmd eq 'getHits') {
  my ($align, $chrom, $start, $stop) = @ARGV;
  my $hits;
  if ($start and $stop) {
    $hits = $client->getHitsRange($align,$chrom,$start,$stop);
  } else {
    $hits = $client->getHits($align,$chrom);
  }
  foreach (@$hits) {
    print "$_\n";
  }
} elsif ($cmd eq 'getWeights') {
  my ($align, $chrom, $start, $stop) = @ARGV;
  my $hits;
  if ($start and $stop) {
    $hits = $client->getWeightsRange($align,$chrom,$start,$stop);
  } else {
    $hits = $client->getWeights($align,$chrom);
  }
  foreach (@$hits) {
    print "$_\n";
  }
}
