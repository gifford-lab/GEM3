#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use Bio::SeqIO;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('core');

my ($species,$version,$chrom,$fasta,$expand);
$expand = 0;
GetOptions("species=s"=>\$species,
	   "version=s"=>\$version,
	   "expand=s"=>\$expand,
	   "chrom=s"=>\$chrom,
	   "fasta=s"=>\$fasta);
$fasta ||= '-';
if ($species =~ /;/) {
  ($species,$version) = split(';',$species);
}
$dbh->{LongReadLen} = 1000000;
$fasta = new Bio::SeqIO(-file=>">$fasta",
			-format=>'fasta');
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
if ($chrom) {
  @results = $dbh->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
  unless (@results) {
    die "Couldn't find chromosome $chrom for $species, $version";
  }
  my $chromid = $results[0];
  my $sth = $dbh->prepare("select sequence from chromsequence where id = $chromid");
  $sth->execute();
  my @r = $sth->fetchrow_array();
  my $seq = new Bio::PrimarySeq(-seq=>$r[0],
				-id=>$chrom);
  $fasta->write_seq($seq);
} else {
  while (my $line = <STDIN>) {
    chomp $line;
    $line =~ s/^\s+//;
    my ($chrom,$start,$stop);
    if ($line =~ /(.+):(\d+)-(\d+)/) {
      ($chrom,$start,$stop) = ($1,$2,$3);
    } else {
      ($chrom,$start,$stop) = split(/\s+/,$line);
    }
    @results = $dbh->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
    unless (@results) {
      die "Couldn't find chromosome $chrom for $species, $version";
    }
    $start -= $expand;
    $stop += $expand;
    my $chromid = $results[0];
    my $length = $stop - $start;
    my $sth = $dbh->prepare("select substr(sequence,$start,$length) from chromsequence where id = $chromid");
    $sth->execute();
    my @r = $sth->fetchrow_array();
    my $seq = new Bio::PrimarySeq(-seq=>$r[0],
				  -id=>"$chrom:$start-$stop");
    $fasta->write_seq($seq);
  }
}
$fasta->close();

