#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Utils;

die "Needs to be updated for new database connection system";

my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
my $passwd = `cat ~/.oracle_passwd`;
chomp($passwd);
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
$dbh->do("alter session set current_schema=arolfe");
my ($analysis,$species,$version);
GetOptions("analysis=s"=>\$analysis,
	   "species=s"=>\$species);

$dbh->{AutoCommit} = 0;
($analysis,$version) = split(';',$analysis);
my ($speciesname, $speciesversion) = split(';',$species);
my @results = $dbh->selectrow_array("select id from species where name = '$speciesname'");
if (@results == 0) {
  die "Can't find species $species";
}
my $speciesid = $results[0];
@results = $dbh->selectrow_array("select id from genome where species = $speciesid and version = '$speciesversion'");
if (@results == 0) {
  die "Can't find genome $speciesname with version $speciesversion";
}
my $genomeid = $results[0];
@results = $dbh->selectrow_array("select id from mleanalysis where name = '$analysis' and version = '$version' and species = $speciesid");
my $analysisid;
if (@results == 0) {
  $dbh->do("insert into mleanalysis (id,species,name,version,active) values (".PSRG::Database::AIinsertText($dbh,'analysis_id').", $speciesid, '$analysis', '$version',1)");
  @results = $dbh->selectrow_array("select ".PSRG::Database::AIfetchValue($dbh,'analysis_id'));
} else {
  warn "Analysis already present : $analysis, $version";
}
$analysisid = $results[0];

my %chrmap;
my $sth = $dbh->prepare("select name, id from chromosome where genome = $genomeid");
$sth->execute();
while (my @row = $sth->fetchrow_array) {
  $chrmap{$row[0]} = $row[1];
}

my $insert = $dbh->prepare("insert into mleresults (analysis, chromosome, position, b_i) values ($analysisid,?,?,?)");
foreach my $fname (@ARGV) {
  open(FILE,$fname) or die "can't open $fname : $!";
  while (<FILE>) {
    chomp;
    my @line = split('\t',$_);
    $line[0] =~ s/chr//;
    $insert->execute($chrmap{$line[0]},int(($line[4] + $line[3])/2),$line[5]);
  }
}
$dbh->commit();
$dbh->disconnect();
