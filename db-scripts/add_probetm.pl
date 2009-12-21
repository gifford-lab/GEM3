#!/usr/bin/perl

BEGIN {
  # for Tiling.pm
  use lib "$ENV{HOME}/projects/probedesign";
}

use strict;
use warnings;
use DBI;
use Bio::SeqIO;
use Getopt::Long;
use Tiling;
use PSRG::Database;
my $dbh = PSRG::Database::handleForRole('chipchip');



# blat file in psl format
my ($blat,$fasta);
GetOptions("blat=s"=>\$blat,
	   "fasta=s"=>\$fasta);

my $seqs = new Bio::SeqIO(-format=>'fasta',
			  -file=>$fasta);

open(BLAT,$blat) or die "Can't open $blat for reading : $!";
while (<BLAT>) {
  last if (/\-{30}/);
}

my $getprobeidsql = 'select id from probedesign where probeid = ?';
my $addprobetmsql = 'insert into probetm(id,tm) values(?,?)';
my $getsth = $dbh->prepare($getprobeidsql);
my $addsth = $dbh->prepare($addprobetmsql);

my $blatline = <BLAT>;
chomp $blatline;
my @blatline = split(/\t/,$blatline);
while (my $seq = $seqs->next_seq()) {
  my $letters = $seq->seq();
  my $name = $seq->id();
  my ($bestscore,$beststart,$beststop) = (0,0,0);  
  do {
    if ($blatline[9] eq $name) {
      my $score = abs($blatline[11] - $blatline[12]);
      if ($score > $bestscore) {
	($bestscore,$beststart,$beststop) = ($score,$blatline[11],$blatline[12]);
      }
      $blatline = <BLAT>;
      if ($blatline) {
	chomp $blatline;
	@blatline = split(/\t/,$blatline);
      }
    }
  } while ($blatline and $blatline[9] eq $name);
  if ($bestscore > 0) {
    my $string = substr($letters,$beststart,$beststop);
    my $tm = Tiling::getTM($string);
    print STDERR "$letters ($beststart:$beststop)\t$string\t$tm\n";
#    $getsth->execute($name);
#    my ($id) = $getsth->fetchrow_array();
#    $addsth->execute($id,$tm);
#    print STDERR "Going to insert $id\t$tm\n";
  }
}


