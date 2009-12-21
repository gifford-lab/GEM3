#!/usr/bin/perl

# This script dumps all of the chip-chip probes for a specified species.  It's meant
# to be used in conjunction with add_probe_locations.pl to map existing array designs
# to a new genome release.
#   The output on STDOUT is FASTA formatted.  The sequence IDs are the internal probe ids.
# An option --arraydesign parameter can be used to limit the output to probes
# from a single array design.

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Database;
my ($species, $designname);
GetOptions("species=s"=>\$species,
	   "arraydesign=s"=>\$designname);

my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');
 
my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
unless (@results) {
  die "No such species : $species";
}
my $speciesid = $results[0];
my $sth = $dbhcore->prepare("select id from genome where species = ?");
$sth->execute($speciesid);
my @genomes;
while (my @results = $sth->fetchrow_array()) {
  push(@genomes,$results[0]);
}

my $designid = undef;
if ($designname) {
  @results = $dbhchip->selectrow_array("select id from arraydesign where name = '$designname'");
  unless (@results) {
    die "No such array design $designname";
  }
  $designid = $results[0];
}

$sth = $dbhchip->prepare("select pd.id, pd.sequence from probedesign pd, arraydesign ad where ad.genome in (" .
			 join(',',@genomes) . ") and ad.id = pd.arraydesign " . 
			 ($designid ? " and pd.arraydesign = $designid" : ""));
$sth->execute();
while (@results = $sth->fetchrow_array()) {
  next unless ($results[1]);
  print ">$results[0]\n$results[1]\n";
}
