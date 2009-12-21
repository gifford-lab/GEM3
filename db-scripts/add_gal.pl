#!/usr/bin/perl

# adds a gal file or tdt file to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::GALFile;
use PSRG::NDFFile;
use PSRG::Utils;
use PSRG::dataset;
use Bio::SearchIO;
use PSRG::Database;
use PSRG::ArrayDesign;
use PSRG::AddGALHandler;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');
$dbhchip->{RowCacheSize} = 100000;
 
my ($species,$version,$design,@blastname,$hashname,@galname);
GetOptions("species=s"=>\$species,
	   "genomeversion=s"=>\$version,
	   "design=s"=>\$design,
	   "blast=s"=>\@blastname,
	   "galfile=s"=>\@galname);

if (not $version) {
  ($species,$version) = split(';',$species);
}
my @hitsfiles = ();
foreach my $blastname (@blastname) {
  if (-e $blastname) {
    print STDERR "Using tophits : $blastname\n";
  } else {  
    $blastname = $blastname . '.tophits';
    unless (-e $blastname) {
      die "Need .tophits file $blastname";
    }
  }
  push(@hitsfiles,$blastname);
}

my $speciesid = PSRG::Database::getSpeciesID($dbhcore,$species);
my $genomeid = PSRG::Database::getGenomeID($dbhcore,$speciesid,$version);
print STDERR "GenomeID is $genomeid\n";
my $designid = PSRG::ArrayDesign::ensureArrayDesign($dbhchip,$design,$genomeid);
print STDERR "DesignID is $designid\n";
$dbhchip->{AutoCommit} = 0;
$dbhchip->{RaiseError} = 1;
my ($galfileid,$pids);
my $designstmt = $dbhchip->prepare("insert into probedesign (id, arraydesign, blockno, colno, rowno, galfile, probename, probeid, type, sequence) ".
				   "values(".PSRG::Database::AIinsertText($dbhchip,'probedesign_id').",?,?,?,?,?,?,?,?,?)");
foreach (@galname) {
  unless (-e $_) {
    warn "$_ doesn't exist";
    next;
  }
  my $galfilename = $_;
  $galfilename =~ s/.*\///;
  if ($galfilename =~ /^(\d+)_D/) {
    $galfilename = "$1.tdt";
  }
  $galfileid = PSRG::ArrayDesign::ensureGalfile($dbhchip,$galfilename);
  
  print STDERR "GALFILE id is $galfileid\n";
  $pids = PSRG::ArrayDesign::getProbeIDs($dbhchip,$galfileid,$designid);

  my $handler = new PSRG::AddGALHandler(species=>$speciesid,
					genome=>$genomeid,
					design=>$designid,
					galfile=>$galfileid,
					pids=>$pids,
					dbh=>$dbhchip,
					insert=>$designstmt);
  my $parser;
  if (/gal$/) {
    $parser = new PSRG::GALFile($_,$handler);
  } elsif (/ndf$/i) {
    $parser = new PSRG::NDFFile($_,$handler);
  } elsif (/tdt$/i) {
    $parser = new PSRG::TDTFile($_,$handler);
  }
  print STDERR "Finished loading galfile $_\n";
  $dbhchip->commit();
  $pids = undef;
  $handler = undef;
  $parser = undef;
  print STDERR "Building probeid -> dbid mapping\n";
  my %probes = ();
  my $sth = $dbhchip->prepare("select id, probeid from probedesign where arraydesign = ? and galfile = ?");
  $sth->execute($designid,$galfileid);
  my $c = 0;
  while (my @r = $sth->fetchrow_array()) {
    push(@{$probes{$r[1]}},$r[0]);
    if ($c++ % 100000 == 0) {
      print STDERR time() . "\n";
    }
  }
  
  $dbhchip->{PrintError} = 0;
  $dbhchip->{PrintWarn} = 0;
  my %chromcache = ();
  my $locstmt = $dbhchip->prepare("insert into probelocation (id, chromosome, startpos, stoppos, loccount, bitscore) values (?,?,?,?,?,?)");
  my %seen = ();
  my %warned = ();
  my $hitsfile = shift(@hitsfiles);
  print STDERR "Reading hitsfile $hitsfile\n";
  open(LOCS,$hitsfile) or die "Can't open $hitsfile : $!";
  my $lastprobe = '';
  my @lines = ();
  my @line;
  my ($chrom,$chromid, $start,$stop,$pos);
  while (<LOCS>) {
    chomp;
    if ($. % 100000 == 0) {
      $dbhchip->commit();
    }
    @line = split(/\t/,$_);
    if ($line[0] ne $lastprobe and @lines and not $seen{$lastprobe}) {
      foreach my $r (@lines) {
	$pos = $r->[0];
	($chrom,$start,$stop) = ($pos =~ /(.+)\:(\d+)\-(\d+)/);
	$chrom =~ s/\.fa$//;
	$chrom =~ s/\.fsa$//;
	$chrom =~ s/\.fasta$//;
	$chrom =~ s/^chr//;
	if ($chrom eq 'Mito') {
	  $chrom = 'mt';
	}
      
	if ($chromcache{$chrom}) {
	  $chromid = $chromcache{$chrom};
	} else {
	  my @results = $dbhcore->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
	  $chromid = $results[0];
	  $chromcache{$chrom} = $chromid;
	  unless ($chromid) {
	    unless ($warned{$chrom}) {
	      warn "Can't find a chromosome for $genomeid, $lastprobe: $chrom : $pos";
	      $warned{$chrom} = 1;
	    }
	    next;
	  }
	}
	unless ($chromid) {
	  warn "Can't get a chromosome for $chrom : $chromid";
	}
	foreach my $pid (@{$probes{$lastprobe}}) {
	  eval {
	    $locstmt->execute($pid,
			      $chromid,
			      $start,
			      $stop,
			      scalar(@lines),
			      $r->[1]);
	  };
	  if ($@ and $@ !~ /duplicate/i and $@ !~ /unique.constraint/) {
	    $dbhchip->rollback();
	    die "$@";
	  }
	}
      }
      $seen{$lastprobe} = 1;
      @lines = ();
    }
    next unless ($probes{$line[0]});
    $lastprobe = $line[0];
    push(@lines,[$line[1],$line[3]]);  
  }
  my ($seen,$unseen) = 0;
  foreach (keys %probes) {
    if ($seen{$_}) {$seen++;} else {$unseen++;}
  }
  $dbhchip->{PrintError} = 1;
  $dbhchip->{PrintWarn} = 1;
  $dbhchip->commit();
  print STDERR "Added probe locations.  $seen had locations and $unseen had no locations.\n";
}

$dbhchip->commit();
$dbhchip->disconnect();

