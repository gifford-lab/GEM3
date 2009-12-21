#!/usr/bin/perl

# generates .data files for the specified species and experiment

use strict;
use warnings;
use DBI;
use Getopt::Long;

my $dbtype = $ENV{PSRGDBTYPE} || 'oracle';
my $dbh;
if ($dbtype eq 'oracle') {
  my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
  my $passwd = `cat ~/.oracle_passwd`;
  chomp($passwd);
  $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
  $dbh->do("alter session set current_schema=arolfe");
  
} elsif ($dbtype eq 'mysql') {
  my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
  if (`hostname` =~ /opteron.csail.mit.edu/) {
    $host = 'localhost';
  }
  my $passwd = `cat ~/.mysql_passwd`;
  chomp($passwd);
  $dbh = DBI->connect("DBI:mysql:host=${host};database=${sid}", $user, $passwd) or
    die "Can't connect as $user/$passwd: $!";
} else {
  die "Unknown database type $dbtype";
}

my $gen = $0;
$gen =~ s/_all//;

my ($version, $species, $args);
GetOptions("version=s"=>\$version,
	   "species=s"=>\$species,
	   "args=s"=>\$args);

my ($speciesid);
my @results = $dbh->selectrow_array("select id from species where name = '$species'");
if (@results == 0) {
  die "Can't find species $species";
}
$speciesid = $results[0];
my $sth = $dbh->prepare("select name from experiment where species = $speciesid and version = '$version'");
$sth->execute();
while (my @row = $sth->fetchrow_array()) {
  print "$gen --species '$species' --exptVersion '$version' --expt '$row[0]' $args\n";
}
