#!/usr/bin/perl

use strict;
use warnings;
use DBI;
my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
my $passwd = `cat ~/.oracle_passwd`;
chomp($passwd);
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
$dbh->do("alter session set current_schema=arolfe");

my $start = 15741763;
my $stop = 16164356;
$dbh->{AutoCommit} = 0;
$dbh->{RaiseError} = 1;

my $j = $start;
my $i = $start + 1000;
my $sth = $dbh->prepare('delete from newprobedesign where id >= ? and id <= ?');
while ($i <= $stop) {
  $sth->execute($j,$i);
  $j = $i;
  $i += 1000;
  $dbh->commit();
  print STDERR "$i ";
  if ($i - $start % 100000 == 0) {print STDERR "\n";}
}
$sth->$dbh->execute($j,$stop);
$dbh->commit();
$dbh->disconnect();
