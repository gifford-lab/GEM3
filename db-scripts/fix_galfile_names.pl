#!/usr/bin/perl

use PSRG::Database;
use strict;
use warnings;
use DBI;
my $dbh = PSRG::Database::handleForRole('chipchip');

my $fetch = $dbh->prepare("select id, name from galfiles");
my $fix = $dbh->prepare("update galfiles set name = ? where id = ?");

$fetch->execute();
while (my ($id, $name) = $fetch->fetchrow_array()) {
  if ($name =~ /^(\d+)_D/) {
    my $newname = "$1.tdt";
    print STDERR "$name -> $newname\n";
    $fix->execute($newname,$id);
  }
}
