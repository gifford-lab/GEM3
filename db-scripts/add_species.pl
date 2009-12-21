#!/usr/bin/perl

# adds a species to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;

my $dbh = PSRG::Database::handleForRole('core');

my ($species);
GetOptions("species=s"=>\$species);
$dbh->do("insert into species values(".PSRG::Database::AIinsertText($dbh,'species_id').",'$species')");

