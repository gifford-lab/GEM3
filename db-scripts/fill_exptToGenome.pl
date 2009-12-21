#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use PSRG::Database;

my $dbhcore = PSRG::Database::handleForRole('core');
my $dbh = PSRG::Database::handleForRole('chipchip');

my $sth = $dbhcore->prepare("select id, genome from chromosome");
$sth->execute();
my %chrommap = ();
while (my @r = $sth->fetchrow_array()) {
  $chrommap{$r[0]} = $r[1];
}

#fixExpt();
#fixAnalysis('mle');
fixAnalysis('bayes');
fixAnalysis('rosetta');

sub fixExpt {
  my $sth = $dbh->prepare("select id from experiment");
  my $gh = $dbh->prepare("select unique(pl.chromosome) from probelocation pl, data d where  d.probe = pl.id and d.experiment = ?");
  my $ih = $dbh->prepare("insert into exptToGenome(experiment,genome) values(?,?)");
  
  $sth->execute();
  while (my @r = $sth->fetchrow_array()) {
    my $exptid = $r[0];
    $gh->execute($exptid);
    my %seen = ();
    while (my @g = $gh->fetchrow_array()) {
      my $chrom = $g[0];
      my $genomeid = $chrommap{$chrom};
      next if ($seen{$genomeid});
      $seen{$genomeid} = 1;
      $ih->execute($exptid,$genomeid);
    }
  }
}

sub fixAnalysis {
  my ($analysis) = @_;
  my $rtable = "${analysis}results";
  my $sth = $dbh->prepare("select id from ${analysis}analysis");
  my $gh = $dbh->prepare("select unique(chromosome) from ${rtable} r where analysis = ?");
  my $ih = $dbh->prepare("insert into ${analysis}ToGenome(analysis,genome) values(?,?)");
  $sth->execute();
  while (my @r = $sth->fetchrow_array()) {
    my $exptid = $r[0];    
    $gh->execute($exptid);
    my %seen = ();
    while (my @g = $gh->fetchrow_array()) {      
      my $chrom = $g[0];
      my $genomeid = $chrommap{$chrom};
      next if ($seen{$genomeid});
      $seen{$genomeid} = 1;
      $ih->execute($exptid,$genomeid);
    }
  }
}
