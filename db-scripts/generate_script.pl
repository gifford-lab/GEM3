#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;

my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
my $passwd = `cat ~/.oracle_passwd`;
chomp($passwd);
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
$dbh->do("alter session set current_schema=arolfe");
my $thisprog = $0;
my $gendata = $thisprog;
$gendata =~ s/script\.pl/datafiles.pl/;
$thisprog =~ /^(.*\/)/;
my $agilentem = $1;
my $loadll = $thisprog;
$loadll =~ s/generate_script/add_ll/;

my ($sqlfilter, $beta, $step, $mean, $var, $maxval, $type, $species, $data);
GetOptions("sqlfilter=s"=>\$sqlfilter,
	   "species=s"=>\$species,
	   "data"=>\$data,
	   "beta=s"=>\$beta,
	   "step=s"=>\$step,
	   "mean=s"=>\$mean,
	   "var=s"=>\$var,
	   "maxval=s"=>\$maxval,
	   "type=s"=>\$type);

my @results = $dbh->selectrow_array("select id from species where name = '$species'");
if (@results == 0) {
  die "Can't find species $species";
}
my $speciesid = $results[0];
if ($sqlfilter) {
  $sqlfilter .= " and species=$speciesid";
} else {
  $sqlfilter = " species=$speciesid";
}
my $sth = $dbh->prepare("select name from experiment where $sqlfilter");
$sth->execute();
while (my ($expt) = $sth->fetchrow_array()) {
  open(SCRIPT,">${expt}.sh") or die "can't open script for $expt : $!";
  if ($data) {
    SCRIPT print "perl $gendata --expt '$expt' --species '$species'";
  }
  my $matname = $expt;
  $namename =~ s/[\.\s]/_/g;
  print SCRIPT "${agilentem}/bin/analyze_chrom --beta $beta --step $step --${type} --mean $mean --var $var --maxval $maxval --coeffs '${expt}.coeffs' --paramsfile '${expt}.params' '${expt}'.*.data > ${matname}.m\n";
  print SCRIPT "matlab -nosplash -nojvm -r ${matname}\n";
  print SCRIPT "perl ${agilentem}/compareModels.pl ${expt}.peak\n";
  print SCRIPT "perl $loadll --expt '${expt}' --analysis '${expt}' --version '${expt}'*.ll\n";
}
