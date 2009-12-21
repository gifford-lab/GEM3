#!/usr/bin/perl

use strict;
use warnings;
use POSIX qw(ceil floor);
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Database;
use PSRG::ArrayDesign;
my $dbhrulers = PSRG::Database::handleForRole('rulers');

# a probeset is a set of probes that already exist.  Each probeset is identified by a name;

my (@probeid, @probesequence, $id, $seq);
my ($name);

GetOptions("name=s"=>\$name,
	   "id"=>\$id,
	   "seq"=>\$seq,
	   "probeid=s"=>\@probeid,
	   "probesequence=s"=>\@probesequence);


$dbhrulers->do("insert into probeset(id,name) values(". PSRG::Database::AIinsertText($dbhrulers,'probeset_id') .
	       ",'$name')");
my @results = $dbhrulers->selectrow_array('select '.PSRG::Database::AIfetchValue($dbhrulers,'probeset_id'));
my $setid = $results[0];

my $insert = $dbhrulers->prepare("insert into probesetentry(pset,probe) values(?,?)");

foreach (@probeid) {
  add_probeid($setid,$_);
}
foreach (@probesequence) {
  add_probesequence($setid,$_);
}
if ($seq) {
  add_probesequence($setid,'-');
}
if ($id) {
  add_probeid($setid,'-');
}

# uses the first column of fname as the probeid
sub add_probeid {
  my ($setid,$fname) = @_;
  open(FILE,$fname) or die "Can't open $fname : $!";
  my $get = $dbhrulers->prepare("select id from probedesign where probeid = ?");
  my %seen = ();
  while (<FILE>) {
    chomp;
    my ($probeid) = split(/\t/,$_);
    next if ($seen{$probeid});
    $seen{$probeid} = 1;
    $get->execute($probeid);    
    while (my ($id) = $get->fetchrow_array()) {
      $insert->execute($setid,$id);
    }
  }
  close FILE;
}

# uses the first column of fname as a probe sequence.  
# compares this to fullsequence
sub add_probesequence {
  my ($setid,$fname) = @_;
  open(FILE,$fname) or die "Can't open $fname : $!";
  my $get = $dbhrulers->prepare("select id from probedesign where fullsequence = ?");
  my %seen = ();
  while (<FILE>) {
    chomp;
    my ($seq) = split(/\t/,$_);
    $get->execute($seq);
    next if ($seen{$seq});
    $seen{$seq} = 1;
    while (my ($id) = $get->fetchrow_array()) {
      $insert->execute($setid,$id);
    }
  }
  close FILE;

}
