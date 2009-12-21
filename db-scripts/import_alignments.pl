#!/usr/bin/perl

# imports alignment information into ucsc multiz style tables.
# this program reads blasttab or blastz formated output
# Input is on STDIN

use strict;
use warnings;
use Getopt::Long;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');

my ($genomequery,$genometarget,
    $queryid, $targetid, $format, $prefix);

GetOptions("query=s"=>\$genomequery,
	   "target=s"=>\$genometarget,
	   "format=s"=>\$format,
	   "prefix=s"=>\$prefix);

unless ($format eq 'blasttab' or
	$format eq 'blastz') {
  die "Unknown format : $format.  I know blasttab and blastz";
}
unless ($prefix) {
  print STDERR "Must supply --prefix.  The standard ucsc blastz prefix is chain.";
  print STDERR "We've also used blast.";
  die "bye";
}
$prefix = '_' . $prefix;

my $dbh = PSRG::Database::handleForRole('ucsc_' . $genometarget);
$dbh->{AutoCommit} = 1;
$dbh->{RaiseError} = 1;
my $tabletemplate = <<END;
CREATE TABLE TABLENAME (
  bin smallint(5) unsigned NOT NULL default '0',
  score double NOT NULL default '0',
  tName varchar(255) NOT NULL default '',
  tSize int(10) unsigned NOT NULL default '0',
  tStart int(10) unsigned NOT NULL default '0',
  tEnd int(10) unsigned NOT NULL default '0',
  qName varchar(255) NOT NULL default '',
  qSize int(10) unsigned NOT NULL default '0',
  qStrand char(1) NOT NULL default '',
  qStart int(10) unsigned NOT NULL default '0',
  qEnd int(10) unsigned NOT NULL default '0',
  id int(10) unsigned NOT NULL default '0',
  KEY bin (bin),
  KEY tStart (tStart),
  KEY tEnd (tEnd),
  KEY id (id)
) TYPE=MyISAM;
END
my %created = ();
my @binOffsets = (512+64+8+1,
		  64+8+1,
		  8+1,
		  1,
		  0);
my $binFirstShift = 17;
my $binNextShift = 3;
my %chromsizes = ();

&cacheChromSizes();

my %insert;
# blasttab is 
# query id, subject id, percent identity, alignment length, mismatches, gap openings, query start, query end, subject start, subject end, e-value, and bit score
my $count = 0;
while (<STDIN>) {
  chomp;
  my ($bin,$tsize,$qstrand, $qsize,$qid,$tid,$percentident,$alignlen,$mismatch,$gaps,$qstart,$qend,$tstart,$tend,$eval,$score);
  if ($format eq 'blasttab') {
    ($qid,$tid,$percentident,$alignlen,$mismatch,$gaps,$qstart,$qend,$tstart,$tend,$eval,$score) = split(/\t/,$_);
    die "Invalid line at $. : $_" unless (defined $score);
    if ($qid eq 'chrM') {$qid = 'chrmt';}
    if ($tid eq 'chrM') {$tid = 'chrmt';}
    $tsize = getChromSize($genometarget,$tid);
    $qsize = getChromSize($genomequery,$qid);
    $qstrand = '+';
  } elsif ($format eq 'blastz') {
    ($score,$tid,$tsize,$tstart,$tend,$qid,$qsize,$qstrand,$qstart,$qend) = split(/\t/,$_);
    if ($qid eq 'chrM') {$qid = 'chrmt';}
    if ($tid eq 'chrM') {$tid = 'chrmt';}
    die "Invalid line at $. : $_" unless (defined $qend);
  }
  unless ($tid =~ /^chr/ or $tid =~ /^scaffold/) {
    $tid = 'chr' . $tid;
  }
  unless ($qid =~ /^chr/ or $qid =~ /^scaffold/) {
    $qid = 'chr' . $qid;
  }
  my $tablename = $tid . $prefix . $genomequery;
  $tablename =~ s/[^\w_]+/_/g;
  check_and_create_table($tablename);
  if (!$insert{$tablename}) {
    $insert{$tablename} = $dbh->prepare("insert into $tablename values(?,?,?,?,?,?,?,?,?,?,?,?)");
  }
  if ($tstart > $tend) {
    $qstrand = '-';
    ($tstart,$tend) = ($tend,$tstart);
  }
  $bin = bin_from_range($tstart,$tend);
  my $sth = $insert{$tablename};
  $sth->execute($bin,
		$score,
		$tid,
		$tsize,
		$tstart,
		$tend,
		$qid,
		$qsize,
		$qstrand,
		$qstart,
		$qend,
		$count++);
}
$dbh->disconnect();


sub cacheChromSizes() {
  my $chromsth = $dbhcore->prepare("select id, name from chromosome where genome = ?");
  my $sizesth = $dbhcore->prepare("select length(sequence) from chromsequence where id = ?");
  foreach my $genome (($genomequery,$genometarget)) {
    my @results = $dbhcore->selectrow_array("select count(*) from genome where version = '$genome'");
    if ($results[0] > 1) {
      die "You lose.  Multiple genomes named $genome";
    }
    @results = $dbhcore->selectrow_array("select id from genome where version = '$genome'");
    my $genomeid = $results[0];
    $chromsth->execute($genomeid);
    while (@results = $chromsth->fetchrow_array()) {
      $sizesth->execute($results[0]);
      my @size = $sizesth->fetchrow_array();
      $chromsizes{$genome}{$results[1]} = $size[0];
    }
  }
}
sub getChromSize {
  my ($genome,$chrom) = @_;
  if (exists $chromsizes{$genome}{$chrom}) {
    return $chromsizes{$genome}{$chrom};
  }
  if ($chrom =~ /^\d/) {
    $chrom = 'chr' . $chrom;
    if (exists $chromsizes{$genome}{$chrom}) {
      return $chromsizes{$genome}{$chrom};
    }
  } elsif ($chrom =~ /^chr/) {
    $chrom =~ s/^chr//;
    if (exists $chromsizes{$genome}{$chrom}) {
      return $chromsizes{$genome}{$chrom};
    }
  } elsif ($chrom eq 'mito') {
    $chrom = 'mt';
    if (exists $chromsizes{$genome}{$chrom}) {
      return $chromsizes{$genome}{$chrom};
    }
  }

  die "No size for $chrom in $genome";
}

sub bin_from_range {
  my ($start,$end) = @_;
  my $startbin = $start;
  my $endbin = $end - 1;
  my $i;
  $startbin >>= $binFirstShift;
  $endbin >>= $binFirstShift;
  for ($i = 0; $i < @binOffsets; ++$i) {
    if ($startbin == $endbin) {
      return $binOffsets[$i] + $startbin;     
    }
    $startbin >>= $binNextShift;
    $endbin >>= $binNextShift;
  }
  warn "Can't get a bin for $start, $end : $startbin, $endbin";
  return 0;
}


#creates the alignment table if it doesn't exist;
sub check_and_create_table {
  my ($tablename) = @_;  
  $tablename =~ s/[^\w_]+/_/g;
  return if ($created{$tablename});
  my @results;
  eval {  # this will die if there's no table
    @results = $dbh->selectrow_array("describe $tablename");
  };
  unless (@results) {
    print STDERR "Creating $tablename\n";
    my $sql = $tabletemplate;
    $sql =~ s/TABLENAME/$tablename/;
#    print STDERR "$sql\n";
    $dbh->do($sql);
    $created{$tablename} = 1;
  } else {
    print STDERR "Clearing $tablename\n";
    $dbh->do("delete from $tablename");
    $created{$tablename} = 1;
  }
}

