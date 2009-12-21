#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use PSRG;
use PSRG::Database;

my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhannotations = PSRG::Database::handleForRole('annotations');

my %chrom2genome = ();
my %scan2genome = ();

my $getgenomes = $dbhcore->prepare("select id, genome from chromosome");
$getgenomes->execute();
while (my @r = $getgenomes->fetchrow_array()) {
  $chrom2genome{$r[0]} = $r[1];
}


my $getchroms = $dbhannotations->prepare("select unique chromosome, scan from wms_hits");
$getchroms->execute();
while (my @r = $getchroms->fetchrow_array()) {
  $scan2genome{$r[1]}{$chrom2genome{$r[0]}} = 1;
}

my $add = $dbhannotations->prepare("insert into wms_scanned_genomes (scan,genome) values (?,?)");
foreach my $scan (keys %scan2genome) {
  for my $genome (keys %{$scan2genome{$scan}}) {
    $add->execute($scan,$genome);
  }
}
