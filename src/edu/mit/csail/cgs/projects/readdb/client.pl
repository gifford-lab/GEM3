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
}
