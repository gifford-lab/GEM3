#!/usr/bin/perl

# this script updates the exptToGenome table to indicate that all of the chipchip experiments
# from one genome version are now available for another.  This happens because you've remapped 
# all of the probes for the relevant species to the new genome.  See dump_probes_for_species.pl and 
# add_probe_locations.pl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Database;
my ($oldspecies,$oldgenome,$newspecies,$newgenome);
GetOptions("old=s"=>\$oldspecies,
	   "new=s"=>\$newspecies);
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');
($oldspecies,$oldgenome) = split(';',$oldspecies);
($newspecies,$newgenome) = split(';',$newspecies);

my @results = $dbhcore->selectrow_array("select id from species where name = '$oldspecies'");
unless (@results) {
  die "No such species : $oldspecies";
}
my $oldspeciesid = $results[0];

@results = $dbhcore->selectrow_array("select id from genome where species = $oldspeciesid and version = '$oldgenome'");
unless (@results) {
  die "No such genome : $oldspecies, $oldgenome";
}
my $oldgenomeid = $results[0];

@results = $dbhcore->selectrow_array("select id from species where name = '$newspecies'");
unless (@results) {
  die "No such species : $newspecies";
}
my $newspeciesid = $results[0];

@results = $dbhcore->selectrow_array("select id from genome where species = $newspeciesid and version = '$newgenome'");
unless (@results) {
  die "No such genome : $newspecies, $newgenome";
}
my $newgenomeid = $results[0];

$dbhchip->{AutoCommit} = 0;
$dbhchip->{RaiseError} = 1;

my $select = $dbhchip->prepare("select experiment, genome from exptToGenome where genome = ?");
my $insert = $dbhchip->prepare("insert into exptToGenome(experiment,genome) values(?,?)");
$select->execute($oldgenomeid);
while (@results = $select->fetchrow_array()) {
  $insert->execute($results[0],$newgenomeid);
}
$dbhchip->commit();
$dbhchip->disconnect();
