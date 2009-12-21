#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use DBD::Oracle qw(:ora_types);
use Getopt::Long;
use PSRG::Utils;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhannot = PSRG::Database::handleForRole('annotations');

my ($species,$version);

GetOptions("species=s"=>\$species);

($species,$version) = split(';',$species);
my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
if (@results == 0) {
  die "Can't find species $species";
}
my $speciesid = $results[0];
@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$version'");
if (@results == 0) {
  die "Can't find genome for $speciesid, $version";
}
my $genomeid = $results[0] or die "Can't get genomeid for $speciesid and $version";

my %chrommap = ();
my $sth = $dbhcore->prepare("select id, name from chromosome where genome = $genomeid");
$sth->execute();
while (@results = $sth->fetchrow_array()) {
  $chrommap{$results[1]} = $results[0];
}

$dbhannot->{AutoCommit} = 0;
$dbhannot->{RaiseError} = 1;

my $insert = $dbhannot->prepare("insert into gff (startdate,id,chromosome,source,type,startpos,stoppos,score,strand,frame,attributes) values " .
				" (DATE '2005-01-25',".PSRG::Database::AIinsertText($dbhannot,'gff_id').",?,?,?,?,?,?,?,?,?)");

my $c = 1;
while (<STDIN>) {
  chomp;
  s/^chr//;
  my @line = split(/\t/,$_);
  if (@line != 9) {
    warn "Wrong number of fields on line $c : @line";
    next;
  }
  $line[0] = $chrommap{$line[0]} || (warn "Can't find chromosome ID for $line[0]" && next);
  if ($line[5] eq '.') {
    $line[5] = 0;
  }
  $insert->execute(@line);
  $c++;
  if ($c % 1000 == 0) {
    $dbh->commit();
    print ".";
    if ($c % 20000 == 0) {
      print "  $c \n";
    }
  }
}
$dbhannot->commit();
$dbhannot->disconnect();

