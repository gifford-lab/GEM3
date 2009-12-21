#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhannot = PSRG::Database::handleForRole('annotations');

# species1 is that query genome and species2 is the target genome

my ($scanname, $scanversion,$species1,$species2,$genome1,$genome2);
GetOptions("scan=s"=>\$scanname,
	   "species1=s"=>\$species1,
	   "species2=s"=>\$species2);

print STDERR "Using $scanname, $scanversion.  $species1, $genome1.  $species2, $genome2\n";
($scanname,$scanversion) = split(';',$scanname);
($species1,$genome1) = split(';',$species1);
($species2,$genome2) = split(';',$species2);

my @results = $dbhcore->selectrow_array("select id from species where name = '$species1'");
if (@results == 0) {
  die "Can't find species ${species1}";
}
my $speciesid1 = $results[0];
@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid1 and version = '$genome1'");
if (@results == 0) {
  die "Can't find genome $genome1 for $species1";
} 
my $genomeid1 = $results[0];

@results = $dbhcore->selectrow_array("select id from species where name = '$species2'");
if (@results == 0) {
  die "Can't find species $species2";
}
my $speciesid2 = $results[0];
@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid2 and version = '$genome2'");
if (@results == 0) {
  die "Can't find genome $genome2 for $species2";
} 
my $genomeid2 = $results[0];


@results = $dbhannot->selectrow_array("select id from syntenyscan where name = '$scanname' and version = '$scanversion'");
my $scanid;
if (@results == 0) {
  $dbhannot->do("insert into syntenyscan(id,name,version) values (".PSRG::Database::AIinsertText($dbhannot,'syntenyscan_id') 
		.",'$scanname','$scanversion')");
  @results = $dbhannot->selectrow_array("select " . PSRG::Database::AIfetchValue($dbhannot,'syntenyscan_id'));
} 
$scanid = $results[0];

@results = $dbhannot->selectrow_array("select count(*) from syntenygenome where scan = $scanid and genome1 = $genomeid1 and genome2 = $genomeid2");
if ($results[0] == 0) {
  $dbhannot->do("insert into syntenygenome(scan,genome1,genome2) values ($scanid,$genomeid1,$genomeid2)");
}

my %chroms1 = ();
my %chroms2 = ();
my $getchroms = $dbhcore->prepare("select name,id from chromosome where genome = ?");
$getchroms->execute($genomeid1);
while (@results = $getchroms->fetchrow_array()) {
  $chroms1{$results[0]} = $results[1];
}
$getchroms->execute($genomeid2);
while (@results = $getchroms->fetchrow_array()) {
  $chroms2{$results[0]} = $results[1];
}

my $insert = $dbhannot->prepare("insert into syntenyblock(scan,genome1,chrom1,start1,stop1,genome2,chrom2,start2,stop2,strand,score) "
				."values(?,?,?,?,?,?,?,?,?,?,?)");

$dbhchip->{AutoCommit} = 0;
$dbhchip->{RaiseError} = 1;
foreach (@ARGV) {
  open(PSL,$_) or die "can't open $_ : $!";
  while (<PSL>) {
    last if (/\-{30}/);
  }
  while (<PSL>) {
    chomp;
    my @l = split(/\t/,$_);
    my ($match,$mismatch,$strand,$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend) =
      @l[0,1,8,9,10,11,12,13,14,15,16];
    my $chrom1 = $chroms1{$qname} || next;
    my $chrom2 = $chroms2{$tname} || next;
    my @vals = ($scanid,$genomeid1,$chrom1,$qstart,$qend,$genomeid2,$chrom2,$tstart,$tend,$strand,$match * 2 - $mismatch);
#    print STDERR "Inserting @vals\n";
    eval {
      $insert->execute(@vals);
    };
    if ($@ and $@ !~ /duplicate/i and $@ !~ /unique.constraint/) {
      $dbhannot->rollback();
      die "$@";
    }
  }
}

$dbhannot->commit();
$dbhannot->disconnect();
