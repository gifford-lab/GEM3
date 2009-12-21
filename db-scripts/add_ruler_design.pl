#!/usr/bin/perl

use strict;
use warnings;
use POSIX qw(ceil floor);
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::GALFile;
use PSRG::Utils;
use Bio::SeqIO;
use PSRG::Database;
use PSRG::ArrayDesign;
use PSRG::AddRulerGALHandler;

my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhrulers = PSRG::Database::handleForRole('rulers');
my $minbitscore = 32;
my ($probetsv,$blasttab, $design, $species, $genomeversion,$allnew);
GetOptions("species=s"=>\$species,
	   "genomeversion=s"=>\$genomeversion,
	   "design=s"=>\$design,
	   "probetsv=s"=>\$probetsv,
	   "blasttab=s"=>\$blasttab,
	   "allnew"=>\$allnew,
	   "minbitscore=s"=>\$minbitscore);

if (not $genomeversion) {
  ($species,$genomeversion) = split(';',$species);
}
unless (-e $probetsv) {
  die "No such file $probetsv for probe information";
}
unless (-e $blasttab) {
  die "No such file $blasttab";
}
my $speciesid = PSRG::Database::getSpeciesID($dbhcore,$species);
my $genomeid = PSRG::Database::getGenomeID($dbhcore,$speciesid,$genomeversion);
print STDERR "GenomeID is $genomeid\n";
my $designid = PSRG::ArrayDesign::ensureArrayDesign($dbhrulers,$design,$genomeid);
print STDERR "DesignID is $designid\n";

$dbhrulers->{AutoCommit} = 0;
$dbhrulers->{RaiseError} = 1;

my %probes = ();
# probetsv has 4 columns: probeid, left part of probe, right part of probe, full sequence.
# if the probe isn't an interval probe, right part is empty.
# left and right parts should only contain active sequence and should not include stilt.
open(TSV,$probetsv) or die "Can't open $probetsv : $!";
while (<TSV>) {
  chomp;
  my @l=split(/\t/,$_);
  $probes{$l[0]} = {interval=>$l[1],
		    left=>$l[2], 
		    right=>$l[3],
		    full=>$l[4]};
}
print STDERR "Loaded probe sequences from $probetsv\n";

my ($galfileid,%pids);
my $designstmt = $dbhrulers->prepare("insert into probedesign (id, arraydesign, blockno, colno, rowno, galfile, probename, probeid, sequenceleft, sequenceright, lengthleft, lengthright, fullsequence, designedinterval) ".
				     "values(".PSRG::Database::AIinsertText($dbhrulers,'probe_id').",?,?,?,?,?,?,?,?,?,?,?,?,?)");


foreach (@ARGV) {
  unless (-e $_) {
    warn "$_ doesn't exist";
    next;
  }
  my $galfilename = $_;
  $galfilename =~ s/.*\///;
  if ($galfilename =~ /^(\d+)_D/) {
    $galfilename = "$1.tdt";
  }
  $galfileid = PSRG::ArrayDesign::ensureGalfile($dbhrulers,$galfilename);
  print STDERR "GALFILE id is $galfileid\n";
  %pids = %{PSRG::ArrayDesign::getProbeIDs($dbhrulers,$galfileid,$designid)};
  my $handler = new PSRG::AddRulerGALHandler(species=>$speciesid,
					     genome=>$genomeid,
					     design=>$designid,
					     galfile=>$galfileid,
					     pids=>\%pids,
					     dbh=>$dbhrulers,
					     insert=>$designstmt,
					     seqs=>\%probes);
  my $parser;
  if (/gal$/) {
    $parser = new PSRG::GALFile($_,$handler);
  } elsif (/ndf$/i) {
    $parser = new PSRG::NDFFile($_,$handler);
  } elsif (/tdt$/i) {
    $parser = new PSRG::TDTFile($_,$handler);
  }
  print STDERR "Finished loading galfile $_\n";
  $dbhrulers->commit();
}
%probes = ();
my $sth = $dbhrulers->prepare("select id, probeid from probedesign where arraydesign = $designid");
$sth->execute();
while (my @r = $sth->fetchrow_array()) {
  push(@{$probes{$r[1]}},$r[0]);
}

my %chromcache = ();
$sth = $dbhcore->prepare("select id, name from chromosome where genome = ?");
$sth->execute($genomeid);
while (my @r = $sth->fetchrow_array()) {
  $chromcache{$r[1]} = $r[0];
}
my $locstmt = $dbhrulers->prepare("insert into probelocation (id, chromosome, startpos, stoppos, probeside) values (?,?,?,?,?)");
my $intstmt = $dbhrulers->prepare("insert into probeinterval (id, chromosome, startpos, stoppos, orientation, leftstart, leftstop, rightstart, rightstop) " .
				  "values (?,?,?,?,?,?,?,?,?)");

open(BLASTTAB,$blasttab) or die "Can't open $blasttab : $!";
my %lefts = ();
my %rights = ();
my $lastprobe = ''; my $lastchrom;
# query id, subject id, percent identity, alignment length, mismatches, gap openings, query start, query end, subject start, subject end, e-value, and bit score
my ($qid,$tid,$percent,$alignlen,$mismatch,$gap,$qstart,$qend,$tstart,$tend,$eval,$score) = (0,1,2,3,4,5,6,7,8,9,10,11);
print STDERR "Reading $blasttab\n";
while (<BLASTTAB>) {
  chomp;
  my @line = split(/\t/,$_);
  next unless ($line[$score] >= $minbitscore and $line[$gap] == 0);
  my ($probeid,$type) = ($line[$qid] =~ /(.*)_(.*)/);
  
  if ($lastprobe and $probeid ne $lastprobe) {
    flushintervals();
  }


  my $chrom = $line[$tid];
  my $chromid;
  $chrom =~ s/^chr0?//;
  $chrom =~ s/\.fs?a//;
  if ($chrom eq 'Mito' or $chrom eq 'M') {
    $chrom = 'mt';
  }
  unless ($chromid = $chromcache{$chrom}) {
    die "Can't get a chromosome for $chrom : $chromid";
  }

  my $strand;
  if ($line[$tstart] >= $line[$tend]) {
    $strand = '+';
  } else {
    $strand = '-';
    ($line[$tstart],$line[$tend]) = ($line[$tend],$line[$tstart]);
  }
  if ($type eq 'left') {
    push(@{$lefts{$chromid}},{score=>$line[$score],
			      strand=>$strand,
			      qstart=>$line[$qstart],
			      qend=>$line[$qend],
			      tstart=>$line[$tstart],
			      tend=>$line[$tend]});
    foreach my $pid (@{$probes{$probeid}}) {
      my @p = ($pid,$chromid,$line[$tstart],$line[$tend],'l');
      eval {
	$locstmt->execute(@p);
      };
      if ($@ and $@ !~ /duplicate/i and $@ !~ /ORA.00001/) {
	$dbhrulers->rollback();
	die "$@";
      }
    } 
  } elsif ($type eq 'right') {
    push(@{$rights{$chromid}},{score=>$line[$score],
			      strand=>$strand,
			      qstart=>$line[$qstart],
			      qend=>$line[$qend],
			      tstart=>$line[$tstart],
			      tend=>$line[$tend]});
    foreach my $pid (@{$probes{$probeid}}) {
      my @p = ($pid,$chromid,$line[$tstart],$line[$tend],'r');
      eval {
	$locstmt->execute(@p);
      };
      if ($@ and $@ !~ /duplicate/i and $@ !~ /ORA.00001/) {
	$dbhrulers->rollback();
	die "$@";
      }
    } 
  } elsif ($type eq 'all') {
    foreach my $pid (@{$probes{$probeid}}) {
      my @p = ($pid,$chromid,$line[$tstart],$line[$tend],'a');
      eval {
	$locstmt->execute(@p);
      };
      if ($@ and $@ !~ /duplicate/i and $@ !~ /ORA.00001/) {
	$dbhrulers->rollback();
	die "$@";
      }
    } 
    my $qmidlower = floor(($line[$qstart] + $line[$qend]) / 2.0);
    my $qmidupper = ceil(($line[$qstart] + $line[$qend]) / 2.0);
    my $tmidlower = floor(($line[$tstart] + $line[$tend]) / 2.0);
    my $tmidupper = ceil(($line[$tstart] + $line[$tend]) / 2.0);
    if ($qmidlower == $qmidupper) {
      $qmidupper++;
    }
    if ($tmidlower == $tmidupper) {
      $tmidupper++;
    }
    if ($strand eq '+') {
      push(@{$lefts{$chromid}},{score=>$line[$score],
				strand=>$strand,
				qstart=>$line[$qstart],
				qend=>$qmidlower,
				tstart=>$line[$tstart],
				tend=>$tmidlower});
      push(@{$rights{$chromid}},{score=>$line[$score],
				 strand=>$strand,
				 qstart=>$qmidupper,
				 qend=>$line[$qend],
				 tstart=>$tmidupper,
				 tend=>$line[$tend]});
    } else {
      push(@{$rights{$chromid}},{score=>$line[$score],
				 strand=>$strand,
				 qstart=>$line[$qend],
				 qend=>$qmidupper,
				 tstart=>$line[$tstart],
				 tend=>$tmidlower});
      push(@{$lefts{$chromid}},{score=>$line[$score],
				strand=>$strand,
				qstart=>$qmidlower,
				qend=>$line[$qstart],
				tstart=>$tmidupper,
				tend=>$line[$tend]});
    }
  } else {
    warn "Unknown type from $probeid";
  }
  $lastprobe = $probeid;
  $lastchrom = $chromid;
}
flushintervals();
print STDERR "Added probe locations\n";
$dbhrulers->commit();
$dbhrulers->disconnect();

sub flushintervals {
#  print STDERR "Adding intervals for $lastprobe: @{$probes{$lastprobe}}\n";
  foreach my $c (keys %lefts) {
#    print STDERR "chrom $c\n";
    foreach my $leftside (@{$lefts{$c}}) {
#      print STDERR "LEFT $leftside\n";
      foreach my $rightside (@{$rights{$c}}) {
#	print STDERR "  $c $leftside->{tstart}, $rightside->{tstart}, $leftside->{strand}, $rightside->{strand}\n";
	next unless (abs($leftside->{tstart} - $rightside->{tstart}) < 100000);
	next unless ($leftside->{strand} eq $rightside->{strand});
	my @p = sort {$a <=> $b} ($leftside->{tstart}, $leftside->{tend},
				  $rightside->{tstart}, $rightside->{tend});
#	print STDERR "\t$c: @p\n";

	my @v = ($c,
		 $p[1],$p[2],
		 $leftside->{tstart} > $rightside->{tstart} ? '-' : '+',
		 $leftside->{qstart}, $leftside->{qend},
		 $rightside->{qstart}, $rightside->{qend});
	foreach my $pid (@{$probes{$lastprobe}}) {
	  eval {
	    $intstmt->execute($pid,@v);
	  };
	  if ($@ and $@ !~ /duplicate/i and $@ !~ /ORA.00001/) {
	    $dbhrulers->rollback();
	    die "$@";
	  }
	}
      }
    }
  }
  %lefts = ();
  %rights = ();
}
