#!/usr/bin/perl

# adds a fragment size distribution to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use Statistics::Lite;
use PSRG::projects::Agilent::ErrorModel;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('chipchip');


my ($distname,$distversion,$description);
GetOptions("name=s"=>\$distname,
	   "version=s"=>\$distversion,
	   "description=s"=>\$description);

my @results = $dbh->selectrow_array("select count(*) from fragdist where name = '$distname' and version = '$distversion'");
if ($results[0] > 0) {
  warn "Fragdist already exists.  Not adding.";
  exit;
}
$dbh->do("insert into fragdist values (".PSRG::Database::AIinsertText($dbh,'fragdist_id').",'$distname','$distversion','$description')");
@results = $dbh->selectrow_array("select id from fragdist where name = '$distname' and version = '$distversion'");
my $distid = $results[0];
my $maxdist = PSRG::projects::Agilent::ErrorModel::getMaximumFragmentLength();
my $settozero = .005;
my @points = ();
foreach my $fname (@ARGV) {
  open(FILE,$fname) or die "can't open $fname : $!";
  my $zeroint = undef;
  while (<FILE>) {
    chomp;
    my ($length,$int) = split(/\s+/,$_);
    if ($length == 0) {
      $zeroint = $int;
    }
    push(@{$points[$length]},$int/$zeroint);
  }
}
for (my $i = $#points; $i > 0; $i--) {
  my $m = Statistics::Lite::mean(@{$points[$i]});
  if ($m < $settozero) {
    $maxdist = $i;
  }
}
for (my $i = 0; $i < $maxdist; $i++) {
  my $val = Statistics::Lite::mean(@{$points[$i]});
  last if ($val < $settozero);
  $val = int($val * 1000)/1000;
  $dbh->do("insert into fragdistentry values ($distid,$i,$val)");
}
