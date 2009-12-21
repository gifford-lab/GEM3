#!/usr/bin/perl

#use strict;
#use warnings;
use DBI;
use DBD::Oracle qw(:ora_types);

my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
my $passwd = `cat ~/.oracle_passwd`;
chomp($passwd);
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
$dbh->do("alter session set current_schema=arolfe");

foreach my $table (@ARGV) {
  my $filename = "${table}.txt";
  print STDERR "Opening $filename\n";
  open(FILE,">$filename") or die "Can't open $filename : $!";
  my $sth = $dbh->prepare("select * from $table");
  $sth->execute();
  while (my @r = $sth->fetchrow_array()) {
    print FILE join("\t",@r) . "\n";
  }
  close FILE;
}
