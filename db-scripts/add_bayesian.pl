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

my ($analysis,$species,$version,@expt,$exptVersion,$genomeVersion, $pattern,$chromindex,$paramfile);
my @dir;


# analysis is the name and version of the results object. 
#  expt is the name of the experiment that was the input.
$chromindex = 1;
GetOptions("analysis=s"=>\$analysis,
	   "species=s"=>\$species,
	   "version=s"=>\$version,
# the filename has compoonents separated by '.'
# chromindex is the zero-based index of the chromosome name in the filename
	   "chromindex=s"=>\$chromindex,
	   "expt=s"=>\@expt,
	   "paramfile=s"=>\$paramfile,
	   "dir=s"=>\@dir,
	   "pattern=s"=>\$pattern);
my ($exptid, $designid,$speciesid,$genomeid);
my %dir;
map {$dir{$_} = 1;} @dir;
@dir = keys %dir;

# make a guess for the parameters file name
if (-e "$dir[0]/bayes.params") {
  $paramfile = "$dir[0]/bayes.params";
  warn "Using default parameters file : $paramfile";
}

($species,$genomeVersion) = split(';',$species);
my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
if (@results == 0) {
  die "Can't find species $species";
}
$speciesid = $results[0];

if ($version and not $analysis) {
  my ($e,$v) = split(';',$expt[0]);
  $analysis = $e;
  warn "Guessing analysis as $analysis ; $version";
} else {
  ($analysis,$version) = split(';',$analysis);
}

@results = $dbhchip->selectrow_array("select id from bayesanalysis where name = '$analysis' and version = '$version' and species = $speciesid");
my $analysisid;
if (@results == 0) {
  $dbhchip->do("insert into bayesanalysis (id,species,name,version,active) values (".PSRG::Database::AIinsertText($dbhchip,'analysis_id').", $speciesid, '$analysis', '$version',1)");
  my ($newexptid) = $dbhchip->selectrow_array("select ".PSRG::Database::AIfetchValue($dbhchip,'analysis_id'));

  foreach my $expt (@expt) {
    my ($exptname,$exptversion,$exptreplicate) = split(';',$expt);
    if ($exptreplicate) {
      @results = $dbhchip->selectrow_array("select id from experiment where name = '$exptname' and species = $speciesid and version='$exptversion' and replicate = '$exptreplicate'");
      if (@results == 0) {
	die "Can't get ID for experiment $expt, $speciesid, $exptVersion";
      }
      $exptid = $results[0];
      $dbhchip->do("insert into bayesanalysisinputs values ($newexptid, $exptid)");
    } else {
      my $sth = $dbhchip->prepare("select id from experiment where name = '$exptname' and species = $speciesid and version='$exptversion'");
      $sth->execute();
      while (@results = $sth->fetchrow_array()) {
	$exptid = $results[0];
	$dbhchip->do("insert into bayesanalysisinputs values ($newexptid, $exptid)");
      }
    }
  }
  $paramfile ||= $analysis . ".params";
  if (-e $paramfile) {
    my $params = PSRG::Utils::readHash(file=>$paramfile);
    my $sth = $dbhchip->prepare("insert into bayesparameters values ($newexptid, ?, ?)");
    foreach (keys %$params) {
      $sth->execute($_,$params->{$_});
    }
  }
  @results = ($newexptid);
} else {
#  die "Analysis already present : $analysis, $version";
  warn "Using existing analysis record : $analysis, $version -> $results[0]";
}
$analysisid = $results[0] or die "Can't get analysis id";

@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$genomeVersion'");
if (@results == 0) {
  die "Can't find genome for $speciesid, $genomeVersion";
}
$genomeid = $results[0] or die "Can't get genomeid for $speciesid and $genomeVersion";

@results = $dbhchip->selectrow_array("select count(*) from bayesToGenome where analysis = $analysisid and genome = $genomeid");
if ($results[0] == 0) {
  $dbhchip->do("insert into bayesToGenome(analysis,genome) values($analysisid,$genomeid)");
}

my $sth = $dbhchip->prepare('insert into bayesresults(analysis,chromosome,position,posterior,posteriorstd,strength,strengthstd) values (?,?,?,?,?,?,?)');

if ($pattern and @dir) {
  @ARGV = ();
  foreach my $dir (@dir) {
    processdir($dir);
  }
}

$dbhchip->{AutoCommit} = 0;
foreach my $posterior (@ARGV) {
  next unless ($posterior =~ /\.jbd$/);
  my $fname = $posterior;
  $fname =~ s/^[^\/]*\///;
  my @name = split(/\./,$fname);
  my $chrom = $name[$chromindex];
  print "FNAME is $fname.  CHROM is $chrom\n";

  @results = $dbhcore->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
  my ($chromid);
  if (@results) {
    $chromid = $results[0];
  } else {
    warn "Can't find chromosome $chrom for genome $genomeid";
    next;
  }

  open(POST,$posterior) or die "Can't open $posterior : $!";
  while (<POST>) {
    chomp;
#    my ($pos,$bm,$bv,$sm,$sv) = split(/\t/,$_);
#    $sth->bind_param(1, $analysisid);
#    $sth->bind_param(2,$chromid);
#    $sth->bind_param(3,$pos);
#    $sth->bind_param(4,$bm);
#    $sth->bind_param(5,$bv);
#    $sth->bind_param(6,$sm);
#    $sth->bind_param(7,$sv);
#    $sth->execute();
    my @l = split(/\t/,$_);
    unless (@l == 5) {
      warn ("Too few fields " . scalar(@l));
    }
    $sth->execute($analysisid,$chromid,@l);
  }
  $dbhchip->commit();
}
$dbhchip->commit();
$dbhchip->disconnect();

sub processdir {
  my ($dir) = @_;
  opendir(DIR,$dir) or die "Can't open dir $dir : $!";
  while ($_ = readdir(DIR)) {
    next if ($_ =~ /^\./);
    if (/$pattern/ and /\.jbd$/) {
      push(@ARGV,"${dir}/${_}");
      #      print STDERR "KEEPING $_\n";
    } elsif (-d "${dir}/${_}") {
      processdir("${dir}/${_}");
    }
  }
}
