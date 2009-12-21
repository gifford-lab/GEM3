#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG::Utils;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');

my ($analysis,$species,@expt,$version,$chromnumkey,@dir,$pattern,$paramfile);

# analysis is the name of the results object. 
#  expt is the name of the experiment that was the input.
#  right now, this only handles a single experiment as input to an
#  analysis and assumes that the experiment version and the
# analysis version are the same
GetOptions("chromindex=s"=>\$chromindex,
	   "analysis=s"=>\$analysis,
	   "species=s"=>\$species,
	   "paramfile=s"=>\$paramfile,
	   "expt=s"=>\@expt,
	   "dir=s"=>\@dir,
	   "pattern=s"=>\$pattern);

$dbhchip->{AutoCommit} = 0;
# make a guess for the parameters file name
if (-e "$dir[0]/mle.params") {
  $paramfile = "$dir[0]/mle.params";
  warn "Using default parameters file : $paramfile";
}

($analysis,$version) = split(';',$analysis);
my ($exptid, $designid,$speciesid,$genomeid);
my ($speciesname, $speciesversion) = split(';',$species);
my @results = $dbhcore->selectrow_array("select id from species where name = '$speciesname'");
if (@results == 0) {
  die "Can't find species $species";
}
$speciesid = $results[0];
@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$speciesversion'");
if (@results == 0) {
  die "Can't find genome $speciesname with version $speciesversion";
}
$genomeid = $results[0];

@results = $dbhchip->selectrow_array("select id from mleanalysis where name = '$analysis' and version = '$version' and species = $speciesid");
my $analysisid;
if (@results == 0) {
  $dbhchip->do("insert into mleanalysis (id,species,name,version,active) values (".PSRG::Database::AIinsertText($dbhchip,'analysis_id').", $speciesid, '$analysis', '$version',1)");
  my ($newexptid) = $dbhchip->selectrow_array("select ".PSRG::Database::AIfetchValue($dbhchip,'analysis_id'));
  foreach my $expt (@expt) {
    my ($exptname,$exptversion,$exptreplicate) = split(';',$expt);
    @results = $dbhchip->selectrow_array("select id from experiment where name = '$exptname' and species = $speciesid and version='$exptversion' and replicate = '$exptreplicate'");
    if (@results == 0) {
      die "Can't get ID for experiment $exptname, $exptversion, $exptreplicate";
    }
    $exptid = $results[0];
    $dbhchip->do("insert into mleanalysisinputs values ($newexptid, $exptid)");
  }
  my $params = PSRG::Utils::readHash(file=>$paramfile);
  my $sth = $dbhchip->prepare("insert into mleparameters values ($newexptid, ?, ?)");
  foreach (keys %$params) {
    $sth->execute($_,$params->{$_});
  }
  @results = ($newexptid);
} else {
  warn "Analysis already present : $analysis, $version";
}
$analysisid = $results[0];

@results = $dbhchip->selectrow_array("select count(*) from mleToGenome where analysis = $analysisid and genome = $genomeid");
if ($results[0] == 0) {
  $dbhchip->do("insert into mleToGenome(analysis,genome) values($analysisid,$genomeid)");
}

if ($pattern and @dir) {
  @ARGV = ();
  foreach my $dir (@dir) {
    opendir(DIR,$dir) or die "Can't open dir $dir : $!";
    while ($_ = readdir(DIR)) {
      if (/$pattern/ and /\.mle$/) {
	push(@ARGV,"${dir}/${_}");
      }
    }
  }
}

foreach my $fname (@ARGV) {
  $fname =~ s/^[^\/]*\///;
  my @fname = split(/\./,$fname);
  my $chrom = $fname[$chromindex];
  @results = $dbhcore->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
  if (@results == 0) {
    warn "Can't find chrom for $chrom in genome $genomeid.  Skipping file $fname";
    next;
  }
  my $chromid = $results[0];
  open(MLE,$fname) or die "Can't open $fname for reading : $!";
  my $sth=$dbhchip->prepare("insert into mleresults values ($analysisid, $chromid, ?, trunc(?,4), trunc(?,2), trunc(?,2), trunc(?,2), trunc(?,4))");
  while (<MLE>) {
    chomp;
    my @line = split("\t",$_);
    next if ($line[1] eq '-');
    if (@line == 2) {
      $sth->execute(@line,-1,-1,-1,-1);
    } else {
      $sth->execute(@line);
    }
  }  
}
$dbhchip->commit();
$dbhchip->disconnect();
