#!/usr/bin/perl

# adds a factor to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');

my ($factor);
GetOptions("factor=s"=>\$factor);
$dbh->do("insert into factors values(".PSRG::Database::AIinsertText($dbh,'factors_id').",'$factor')");



