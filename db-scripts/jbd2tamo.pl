#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use DBD::Oracle qw(:ora_types);
use Getopt::Long;
use PSRG::Utils;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');

# species and genome version
my ($species,$genomeversion);
# analysis name and version
my ($analysis,$analysisversion);
# any point with a posterior probability of binding greater than this value
# is considered in a bound region
my $posteriorbackground = .1;
# a bound region must have at least one point with a posterior greater than this
# value to be included in the output
my $posteriorthresh = .2;
# a bound region must have a total "size" (sum of posterior * strength) greater than
# this value to be included in the output
my $sizethresh = 2;
# pull this much sequence around each binding event
# if windowsize is <= 1, then it's interpreted as a posterior and the window
# around the binding event in which the posterior is greater than this value is kept.
my $windowsize = 150;
# base filename for output.  We'll produce .bound.fasta and .bound.score files based on this name
my $basename;
# maximum distance between two points in the bayes results that will
# still be considered part of the same region.  The default is 
# either taken from the bayesparams table (the unit parameter)
# or, if that's not present, 30.
my $unit;
# if defined, randomly select this many binding events to dump out.  Otherwise, prints 
# all the ones that meet the criteria
my $selectrandom = undef;

GetOptions("analysis=s"=>\$analysis,
	   "species=s"=>\$species,
	   "posteriorbackground=s"=>\$posteriorbackground,
	   "posteriorthresh=s"=>\$posteriorthresh,
	   "unit=s"=>\$unit,
	   "sizethresh=s"=>\$sizethresh,
	   "windowsize=s"=>\$windowsize,
	   "basename=s"=>\$basename,
	   "selectrandom=s"=>\$selectrandom);

$dbhcore->{LongReadLen} = 1000000;

($species,$genomeversion) = split(';',$species);
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

($analysis,$analysisversion) = split(';',$analysis);
@results = $dbhchip->selectrow_array("select id from bayesanalysis where name = '$analysis' and version = '$analysisversion' and species = $speciesid");
my $analysisid;
if (@results == 0) {
  die "Can't find analysis $analysis, $analysisversion";
} else {
  $analysisid = $results[0];
}

@results = $dbhchip->selectrow_array("select value from bayesparameters where analysis = $analysisid and name = 'unit'");
my $analysisunit = $results[0] || 30;
if (@results && not defined($unit)) {
  $unit = $analysisunit;
}

if ($windowsize < $posteriorbackground) {
  warn "Setting windowsize from $windowsize to $posteriorbackground";
  $windowsize = $posteriorbackground;
}

my %chroms = ();
my $sth = $dbhcore->prepare("select id, name from chromosome where genome = $genomeid");
$sth->execute();
while (my @row = $sth->fetchrow_array()) {
  $chroms{$row[0]} = $row[1];
  print STDERR "Got chrom $row[0] -> $row[1]\n";
}
my $chromstring = '(' . join(',',keys %chroms) . ')';

my $fastaout = $basename . '.bound.fasta';
my $scoreout = $basename . '.bound.score';

my $seqsth = $dbhcore->prepare("select substr(sequence,?,?) from chromsequence where id = ?");
my $lensth = $dbhcore->prepare("select length(sequence) from chromsequence where id = ?");

$sth = $dbhchip->prepare("select chromosome, position, posterior, strength from bayesresults where "
			 . "analysis = ? and posterior > ? and chromosome in $chromstring order by chromosome, position");
my $allsth = $dbhchip->prepare("select chromosome, position, posterior, strength from bayesresults where "
			       . "analysis = ? and chromosome = ? and position >= ? and position <= ?");
$sth->execute($analysisid,$posteriorbackground);
my ($lastchrom, $lastpos, $count) = (-1,-1,0);
my @toprint = ();
my $elevated = [];
while (my @r = $sth->fetchrow_array()) {
  if ($r[0] != $lastchrom) {
#    print STDERR "CHROM is $r[0], $chroms{$r[0]}\n";
  }
  if ($r[0] != $lastchrom or
      $r[1] > $lastpos + $unit) {
#    print STDERR "Skip from $lastpos to $r[1]\n";
    push(@toprint,$elevated);
    $elevated = [];
  }
#  print STDERR "  saw @r\n";
  push(@$elevated,\@r);
  ($lastchrom,$lastpos) = @r;
}
if (@$elevated) {
  push(@toprint,$elevated);
}
if (defined $selectrandom) {
  @toprint = sort { rand() <=> rand()} @toprint;
} else {
  $selectrandom = scalar(@toprint);
}
open(FASTA,">$fastaout") or die "can't open $fastaout for writing : $!";
open(SCORE,">$scoreout") or die "can't open $scoreout for writing : $!";
my $printed = 0;
while (@toprint and ($printed <= $selectrandom)) {
  my $region = shift(@toprint);
  $printed += printbinding(@{$region});
}
close FASTA;
close SCORE;

# elements of @_ are arrayrefs with fields
# 0: chrom
# 1: pos
# 2: posterior
# 3: strength
sub printbinding {
  my $sum = 0;
  my $possum = 0;
  my $postsum = 0;
  my $maxposterior = 0;
  my $chromid = undef;
  return 0 unless (@_);
  foreach (@_) {
    unless ($chromid) {
      $chromid = $_->[0];
    }
    if ($_->[2] > $maxposterior) {
      $maxposterior = $_->[2];
    }
    $sum += $_->[2] * $_->[3];
    $possum += $_->[2] * $_->[3] * $_->[1];
  }
  return 0 unless  ($maxposterior >= $posteriorthresh and
		  $sum >= $sizethresh);
  my $center = int($possum / $sum);
  my ($min,$max);
  # if windowsize is in bases, just compute max and min
  # based on the center and the windowsize
  if ($windowsize > 1) {    
    $min = int($center - $windowsize/2);
    $max = int($center + $windowsize/2) - 1;
  } else {
    # if windowsize is < 1, then it's a posterior value.
    # max and min are the coordinates of the region around the center
    # in which the posterior is greater than this value.
    my $centerindex = 0;
    for (my $i = 0; $i < @_; $i++) {
      if (abs($center - $_[$i][1]) < 
	  abs($center - $_[$centerindex][1])) {
	$centerindex = $i;
      }
    }
#    print STDERR "Datapoints go from $_[0][1] to $_[-1][1]\n";
#    print STDERR "Centerindex is $centerindex.  Center is $center\n";    
    my $i = $centerindex;
    while ($_[$i][2] >= $windowsize and $i < $#_) {
      $i++;
    }
    $max = $_[$i][1];
    $i = $centerindex;
    while ($_[$i][2] >= $windowsize and $i > 0) {
      $i--;
    }
    $min = $_[$i][1];
    if (abs($max - $min) < 2 * $analysisunit) {
      $max += $analysisunit;
      $min -= $analysisunit;
    }
  }
  if ($min < 1) {$min = 1;}
  $lensth->execute($chromid);
  my ($strlen) = $lensth->fetchrow_array();
  if ($max > $strlen) {$max = $strlen;}
  unless ($chroms{$chromid}) {die "No such chromid : $chromid";}
  my $outname = $chroms{$chromid} . ':' . $min . '-' . $max;
  $seqsth->execute($min,($max - $min + 1),$chromid);
  my ($outstring) = $seqsth->fetchrow_array();
  unless (length($outstring) == ($max - $min + 1)) {
    print STDERR "Length is " . length($outstring) . "\n";
    die "expected $min to $max";
  }
  my @l = ();
  $allsth->execute($analysisid, $chromid,$min,$max);
  while (my @r = $allsth->fetchrow_array()) {
    push(@l,\@r);
  }
  return 0 unless (@l);
  print FASTA ">$outname\n";
  while ($outstring) {
    my $substr = substr($outstring,0,60,'');
    print FASTA "$substr\n";
  }
  print SCORE ">$outname\n";
  my @score = ();
#  print STDERR "Window is from $min to $max\n";
  for (my $i = $min; $i < $l[0]->[1]; $i++) {
#    print STDERR "Adding zero for $i\n";
    push(@score,0);
  }
  for (my $j = 0; $j < $#l; $j++) {
    for (my $i = $l[$j][1]; $i < $l[$j+1][1]; $i++) {
      push(@score,$l[$j][2] + ($l[$j+1][2] - $l[$j][2]) *
	   (($i - $l[$j][1])/($l[$j+1][1] - $l[$j][1])));
#      print STDERR "$i (at $l[$j][1]) : $score[-1]\n";
    }
  }
  push(@score,$l[-1][2]);
  for (my $i = $l[-1][1] + 1; $i <= $max; $i++) {
#    print STDERR "Adding zero for $i\n";
    push(@score,0);
  }
  $sum = 0;
  foreach (@score) {$sum += $_;}
  @score = map {sprintf("%0.6f",$_ * $maxposterior / $sum)} @score;
  if (@score != ($max - $min + 1)) {
    print STDERR "scores length is " . scalar(@score) . " but length of outstring is " . length($outstring) . "\n";
    print STDERR "posteriors are " . join(" ",map {$_->[2]} @l) . "\n";
    print STDERR "strengths are " . join(" ",map {$_->[3]} @l) . "\n";
    die "Got wrong number of bytes";
  }
  while (@score) {
    my @subscore = splice(@score,0,60);
    print SCORE join(" ",@subscore) . " \n";
  }
  return 1;
}

