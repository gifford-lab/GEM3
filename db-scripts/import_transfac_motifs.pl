#!/usr/bin/perl

use strict;
use warnings;

use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhannot = PSRG::Database::handleForRole('annotations');
$dbhannot->{AutoCommit} = 0;

my %background = ('Homo sapiens'=>{A=>.295,C=>.205,G=>.205,T=>.295},
		  'Mus musculus'=>{A=>.291,C=>.209,G=>.209,T=>.291},
		  'Drosophila melanogaster'=>{A=>.29,C=>.21,G=>.21,T=>.29},
		  'Saccharomyces cerevisiae'=>{A=>.310,C=>.190,G=>.190,T=>.310});

# TRANSFAC refers to species, so cache a list
my %species;
my $sth = $dbhcore->prepare("select id, name from species");
$sth->execute();
while (my @r = $sth->fetchrow_array()) {
  $species{$r[1]} = $r[0];
}

my $addmatrix = $dbhannot->prepare("insert into weightmatrix(id,species,name,version,type) values (" .
				  PSRG::Database::AIinsertText($dbhannot,"weightmatrix_id") .
				  ",?,?,?,?)");
my $addcolumn = $dbhannot->prepare("insert into weightmatrixcols (weightmatrix,position,letter,weight) values".
				  "(?,?,?,?)");
my $transfac_version = undef;

my ($name,$acc,$species, $version);
my $inmatrix = 0;
my @matrix;
while (<STDIN>) {
  chomp;
  if (/^VV.*Release\s+([\d\.]+)\s/) {
    $transfac_version = $1;
    next;
  }
  if (/^\/\//) {
    addmatrix();
    ($name,$species,$version) = (undef,undef,undef,undef);
    $inmatrix = 0;
    @matrix = ();
    next;    
  }
  if (/^NA\s+(.*)/) {
    $name = $1;
    next;
  }
  if (/^AC\s+(.*)/) {
    $version = "TRANSFAC ${transfac_version}, $1";
    next;
  }
  if (/^BF\s+.*:(.*)\./ and
      not $species) {
    $species = $1;
    $species =~ s/^.*,\s*//;
    next;
  }
  if (/^P0/) {
    unless (/P0\s*A\s*C\s*G\s*T/) {
      die "Someone switched the base order: $_";
    }
    $inmatrix = 1;
  }
  if (/^\d+\s*/ and $inmatrix) {
    my @line = split(/\s+/,$');
#    print STDERR "PUSHING @line\n";
    push(@matrix,{A=>$line[0],
		  C=>$line[1],
		  G=>$line[2],
		  T=>$line[3]});
  }
  
}
addmatrix();

sub addmatrix {
  unless (@matrix and $name and $species and $version) {
    return;
  }
  unless ($species{$species}) {
    return;
  }
  $addmatrix->execute($species{$species}, $name, $version, 'TRANSFAC');
  my ($matrixid) = $dbhannot->selectrow_array("select ".PSRG::Database::AIfetchValue($dbhannot,'weightmatrix_id'));
#  print STDERR "$name, $version\n";
  my $count = 0;
  foreach my $letter (keys %{$matrix[0]}) {
    $count += $matrix[0]{$letter};
  }
#  print STDERR "Count is $count\n";
  my %bg = %{$background{$species} || $background{'Homo sapiens'}};
  for (my $i = 0; $i < @matrix; $i++) {
    foreach my $letter (keys %{$matrix[$i]}) {
      my $prob = ($matrix[$i]{$letter}) / $count;
      if ($prob < .0001) {
	$prob = .0001;
      }
      my $logodds = log($prob/$bg{$letter})  / log(2);
      $addcolumn->execute($matrixid,$i,$letter,$logodds);
#      print STDERR "$i, $letter => $prob, $logodds\n";
    }
  }
#  $dbhannot->rollback();
#  exit;
  $dbhannot->commit();
}
