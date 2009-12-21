#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');


my ($species, $genomeversion, $analysis);

GetOptions("analysis=s"=>\$analysis,
	   "species=s"=>\$species);

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
my $genomeid = $results[0];

my ($analysisname, $analysisversion) = split(';',$analysis);
@results = $dbhchip->selectrow_array("select id from rosettaanalysis where name = '$analysisname' and version = '$analysisversion' and species = $speciesid");
if (@results == 0) {
  $dbhchip->do("insert into rosettaanalysis (id, species, name, version, active) values (".PSRG::Database::AIinsertText($dbhchip,'rosettaanalysis_id').", $speciesid, '$analysisname', '$analysisversion',1)");
  @results = $dbhchip->selectrow_array("select ".PSRG::Database::AIfetchValue($dbhchip,'rosettaanalysis_id'));
} else {
  warn "Using existing analysis record";
}
my $analysisid = $results[0];

@results = $dbhchip->selectrow_array("select count(*) from rosettaToGenome where analysis = $analysisid and genome = $genomeid");
if ($results[0] == 0) {
  $dbhchip->do("insert into rosettaToGenome(analysis,genome) values($analysisid,$genomeid)");
}

my %chrmap;
my $sth = $dbhcore->prepare("select name, id from chromosome where genome = $genomeid");
$sth->execute();
while (my @row = $sth->fetchrow_array) {
  $chrmap{$row[0]} = $row[1];
}

my %mspchr;
if ($species eq 'Homo sapiens') {
  %mspchr = (23=>'X',
	     24=>'Y',
	     25=>'mt');
} elsif ($species eq 'Mus musculus') {
  %mspchr = (20=>'X',
	     21=>'Y',
	     22=>'mt');
} elsif ($species eq 'Danio rerio') {
  %mspchr = (26=>'Un',
	     27=>'NA');
} elsif ($species eq 'Drosophila melanogaster') {
  %mspchr = (1=>'2L',
	      2=>'2R',
	      3=>'3L',
	      4=>4,
	      5=>'3R',
	      6=>'X',
	      7=>'Yh');
  foreach (keys %chrmap) {
    next if (/^chr/i);
    $mspchr{"chr${_}"} = $_;
  }
}


$dbhchip->{AutoCommit} = 0;
$dbhchip->{RaiseError} = 1;
$sth = $dbhchip->prepare("insert into rosettaresults(analysis,chromosome,position,ratio,x,pval,pval3,red,green,medianofratios) " .
		     "values ($analysisid,?,?,?,?,?,?,?,?,?)");
my @wantcols = qw(chr pos ratio x pval1 pval3 red green medianofratios);
$wantcols[3] = "x'";

foreach my $filename (@ARGV) {
  print STDERR "Reading $filename\n";
  open(FILE,$filename) or die "Can't open $filename : $!";
  my $header = <FILE>;
  chomp($header);
  my @cols = split(/\t/,$header);
  my %cols = ();
  for (my $i = 0; $i <= $#cols; $i++) {
    $cols{lc($cols[$i])} = $i;
  }
  foreach (@wantcols) {
    die "Can't find column $_ in @cols" unless (defined $cols{$_});
  }
  while (my $line = <FILE>) {
    chomp $line;
    next if ($line =~ /FLAG/);
    my @line = split(/\t/,$line);
    if (exists $mspchr{$line[$cols{chr}]}) {
      $line[$cols{chr}] = $mspchr{$line[$cols{chr}]};
    }
    $line[$cols{chr}] = $chrmap{$line[$cols{chr}]} || (warn "No chrom for $line[$cols{chr}]", next);
#    for the H3 depletion file that stuart had
#    $line[$cols{ratio}] = 1/$line[$cols{ratio}];
    eval {
      $sth->execute(map {$line[$cols{$_}]} @wantcols);
    };
    if ($@) {
      if ($@ !~ /ORA-00001/) {
	die $@;
      }
    }
  }

  $dbhchip->commit();
}
$dbhchip->disconnect();
