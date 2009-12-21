#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhannot = PSRG::Database::handleForRole('annotation');


my ($species,$genomeversion);
GetOptions("species=s"=>\$species);
if ($species =~ /;/) {
  ($species,$genomeversion) = split(/;/,$species);
}
my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
if (@results == 0) {
  die "Can't find species $species";
}
my $speciesid = $results[0];
@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$genomeversion'");
if (@results == 0) {
  die "Can't find genome for $speciesid, $genomeversion";
}
my $genomeid = $results[0] or die "Can't get genomeid for $speciesid and $genomeversion";


my %chrommap = ();
my $sth = $dbhcore->prepare("select id, name from chromosome where genome = $genomeid");
$sth->execute();
while (@results = $sth->fetchrow_array()) {
  $chrommap{$results[1]} = $results[0];
}

my $insert = $dbhannot->prepare("insert into liver_conservation_motifs(chromosome,factor,position,score) values (?,?,?,?)");

foreach (@ARGV) {
  my @name = split('\.',$_);
  my $factor = $name[0];
  open(SCAN,$_) or die "Can't open $_ : $!";
  <SCAN>;<SCAN>;
  while (my $line = <SCAN>) {
    chomp $line;
    my ($gene,@positions) = split(/\s+/,$line);
    my ($name,$chrom,$tss,$rstart,$rstop) = split(/:/,$gene);
    my $chromid = $chrommap{$chrom};
    foreach my $pos (@positions) {
      my ($motifpos,$score) = split(/:/,$pos);
      eval {
	$insert->execute($chromid,$factor,$motifpos,$score);
      };
      if ($@) {
	if ($@ !~ /ORA-00001/) {
	  die $@;
	}
      }
    }
  }
  close SCAN;
}
