#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhannot = PSRG::Database::handleForRole('annotation');


my $species = 'Saccharomyces cerevisiae';
my $genomeVersion = '1-25-05';
my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
if (@results == 0) {
  die "Can't find species $species";
}
my $speciesid = $results[0];
@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$genomeVersion'");
unless (@results) {
  die "Couldn't get genome from $species, $genomeVersion";
}
my $genomeid = $results[0];

my $sth = $dbhcore->prepare("select name, id from chromosome where genome = $genomeid");
$sth->execute();
my %chroms = ();
while (@results = $sth->fetchrow_array()) {
  $chroms{$results[0]} = $results[1];
}


my $probefname = $ARGV[0];
my $callfname = $ARGV[1];

open(PROBES,$probefname) or die "Can't open $probefname : $!";
my %probes = ();
<PROBES>;
while (<PROBES>) {
  chomp;
  my @line = split(/\t/,$_);
  $probes{$line[0]} = {chrom=>$line[1],
		       start=>$line[2],
		       stop=>$line[3]};  
}
close PROBES;

my @fields = qw(startdate stopdate id chromosome source type startpos stoppos score strand frame);
my @values = ("DATE '2005-07-06'","DATE '2000-01-01'",'gff_id.nextval','?',"'Yuan/Rando 2005'","'Nucleosome'",'?','?',
	      '?',"' '","'.'");
$sth = $dbhannot->prepare("insert into gff (" . join(",",@fields) . ") values (" . join(",",@values) .")");

open(CALLS,$callfname) or die "Can't open $callfname : $!";
<CALLS>;
while (<CALLS>) {
  chomp;
  my @line = split(/\t/,$_);
  my $probeinfo = $probes{$line[0]} || next;
  my $score = ($line[3] + $line[5] + $line[7])/3;
  $sth->execute($chroms{$probeinfo->{chrom}},
		$probeinfo->{start},
		$probeinfo->{stop},
		$score);
}
