#!/usr/bin/perl

# adds a chromosome to the database

use strict;
use warnings;
use DBI;
use DBD::mysql;
use Getopt::Long;
use Bio::SeqIO;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');

# if --chrom is given, then read only the first sequence from the fasta file and name it accordingly.
# otherwise, read the entire fasta file and use the names it contains.

my ($species,$version,$chrom,$fasta);
GetOptions("species=s"=>\$species,
	   "version=s"=>\$version,
	   "chrom=s"=>\$chrom,
	   "fasta=s"=>\$fasta);
my @results = $dbh->selectrow_array("select id from species where name = '$species'");
unless (@results) {
  die "Can't find species $species";
}
my $speciesid = $results[0];
@results = $dbh->selectrow_array("select id from genome where species = $speciesid and version = '$version'");
unless (@results) {
  die "Couldn't get genome from $species, $version";
}
my $genomeid = $results[0];

print STDERR "Getting sequence from $fasta\n";
$fasta = new Bio::SeqIO(-file=>$fasta,
			-format=>'fasta');
my $sth = $dbh->prepare("insert into chromsequence values (?, ?)");
if ($chrom) {
  @results = $dbh->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
  unless (@results) {
    die "Couldn't find chromosome $chrom for $species, $version";
  }
  my $chromid = $results[0];
  my $seq = $fasta->next_seq();
  $sth->execute($chromid,$seq->seq());
} else {
  while (my $seq = $fasta->next_seq()) {
    $chrom = $seq->id();
    @results = $dbh->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
    unless (@results) {
      die "Couldn't find chromosome $chrom for $species, $version";
    }
    my $chromid = $results[0];    
    print STDERR "Adding $chrom, $chromid\n";
    $sth->execute($chromid,$seq->seq());
  }
}
