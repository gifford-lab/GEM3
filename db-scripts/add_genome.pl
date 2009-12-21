#!/usr/bin/perl

# adds a genome version to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');

my ($species,$version,$desc);
GetOptions("species=s"=>\$species,
	   "version=s"=>\$version,
	   "description=s"=>\$desc);
my @results = $dbh->selectrow_array("select id from species where name = '$species'");
my $speciesid = $results[0];
$dbh->do("insert into genome (id,species,version,description) values (".PSRG::Database::AIinsertText($dbh,'genome_id').",$speciesid,'$version','$desc')");
