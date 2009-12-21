#!/usr/bin/perl

use strict;
use warnings;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');
my $schema = $dbh->{Username} || die "no username in $dbh";  

my $sth = $dbh->prepare("select table_name from sys.user_tables where table_name like '\%TMP\%' or table_name like '\%TEMP%\'");
$sth->execute();
while (my @r = $sth->fetchrow_array()) {
  print "$r[0]\n";
  $dbh->do("drop table ${schema}.$r[0]");
}
