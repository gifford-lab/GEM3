#!/usr/bin/perl

# adds a chromosome to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');

my ($species,$version,$chrom);
GetOptions("species=s"=>\$species,
	   "version=s"=>\$version,
	   "chrom=s"=>\$chrom);
my @results = $dbh->selectrow_array("select id from species where name = '$species'");
unless (@results) {
  die "Can't find species $species";
}
my $speciesid = $results[0];
@results = $dbh->selectrow_array("select id from genome where species = $speciesid and version = '$version'");
unless (@results) {
  die "Couldn't get genome from $species ($speciesid) , $version";
}
my $genomeid = $results[0];
$dbh->do("insert into chromosome values (".PSRG::Database::AIinsertText($dbh,'chromosome_id').", '$chrom', $genomeid)");
