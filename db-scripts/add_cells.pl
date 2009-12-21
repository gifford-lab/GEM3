#!/usr/bin/perl

# adds a cells to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');

my ($cells);
GetOptions("cells=s"=>\$cells);
$dbh->do("insert into cells values(".PSRG::Database::AIinsertText($dbh,'cells_id').",'$cells')");



