#!/usr/bin/perl

# adds a species to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');

my ($condition);
GetOptions("condition=s"=>\$condition);
$dbh->do("insert into condition values(condition_id.nextval,'$condition')");



