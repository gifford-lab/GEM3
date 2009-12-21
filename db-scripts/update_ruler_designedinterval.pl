#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Database;
use PSRG::ArrayDesign;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhrulers = PSRG::Database::handleForRole('rulers');
my ($probetsv, $design, $species, $genomeversion);
GetOptions("species=s"=>\$species,
	   "genomeversion=s"=>\$genomeversion,
	   "design=s"=>\$design,
	   "probetsv=s"=>\$probetsv);
if (not $genomeversion) {
  ($species,$genomeversion) = split(';',$species);
}
my $speciesid = PSRG::Database::getSpeciesID($dbhcore,$species);
my $genomeid = PSRG::Database::getGenomeID($dbhcore,$speciesid,$genomeversion);
print STDERR "GenomeID is $genomeid\n";
my $designid = PSRG::ArrayDesign::ensureArrayDesign($dbhrulers,$design,$genomeid);
print STDERR "DesignID is $designid\n";
$dbhrulers->{AutoCommit} = 0;
$dbhrulers->{RaiseError} = 1;

my $stmt = $dbhrulers->prepare("update probedesign set designedinterval = ?, fullsequence = ? where probeid = ? and arraydesign = ?");
open(TSV,$probetsv) or die "Can't open $probetsv : $!";
while (<TSV>) {
  chomp;
  my @l=split(/\t/,$_);
  next unless ($l[1] =~ /^\d+$/);
  $stmt->execute($l[1],$l[-1],$l[0],$designid);
}
close TSV;
$dbhrulers->commit();
$dbhrulers->disconnect();
