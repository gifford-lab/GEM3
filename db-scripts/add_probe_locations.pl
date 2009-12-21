#!/usr/bin/perl

# This script loads probe locations from psl output.  The query sequence name
# *MUST* be the database id of the probe- dump_probes_for_species.pl will create
# a FASTA file with correctly named probe sequences.

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Database;
my ($species, $version,$minscore);
$minscore = 100;
GetOptions("species=s"=>\$species,
	   "genomeversion=s"=>\$version,
	   "minscore=s"=>\$minscore);

my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');

if (not $version) {
  ($species,$version) = split(';',$species);
}

my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
unless (@results) {
  die "No such species : $species";
}
my $speciesid = $results[0];

@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$version'");
unless (@results) {
  die "No such genome : $species, $version";
}
my $genomeid = $results[0];

my %chromcache = ();
my $chromsth = $dbhcore->prepare("select id, name from chromosome where genome = ?");
$chromsth->execute($genomeid);
while (@results = $chromsth->fetchrow_array()) {
  $chromcache{$results[1]} = $results[0];
}

my $addpos = $dbhchip->prepare("insert into probelocation(id,chromosome,startpos,stoppos,strand,loccount,bitscore) values (?,?,?,?,?,?,?)");
my @rows = ();
my $lastprobe = -1;
$dbhchip->{AutoCommit} = 0;
$dbhchip->{RaiseError} = 1;
while (<STDIN>) {
  last if (/\-{30}/);
}
while (<STDIN>) {
  while (<STDIN>) {
    chomp;
    my @line = split(/\t/,$_);
    # 0        1       2       3        4        5      6          7      8     9              19      11      12      13               14     15      16     17         18             19         29
    #match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
    #     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
    my ($qname,$tname,$tstart,$tend,$qstart,$qend,$strand) = @line[(9,13,15,16,11,12,8)];
    my $score = 2 * $line[0] - $line[1];
    unless (defined $score) {
      die "Invalid row : $_";
    }
    next unless ($score >= $minscore);
    ($tstart,$tend) = sort {$a <=> $b} ($tstart,$tend);
    if ($lastprobe > 0 and
	$qname != $lastprobe) {
      insert_rows();
    }
    $tname =~ s/^chr//;
    $tname =~ s/\.fa//;
    $tname =~ s/\.fasta//;
    unless ($chromcache{$tname}) {
      warn "No Chromosome ID for $tname";
      next;
    }
    $lastprobe = $qname;
    push(@rows,[$qname,$chromcache{$tname},$tstart,$tend,$strand,$score]);
  }
}
insert_rows();
$dbhchip->commit();
$dbhchip->disconnect();

sub insert_rows {
  my $count = scalar(@rows);
  foreach my $row (@rows) {
    eval {
      $addpos->execute($row->[0],$row->[1],$row->[2],$row->[3],$row->[4],$count,$row->[5]);
    };
    if ($@ and $@ !~ /duplicate/i and $@ !~ /unique.constraint/) {
      $dbhchip->rollback();
      die "$@";
    }
  }
  @rows = ();
}
