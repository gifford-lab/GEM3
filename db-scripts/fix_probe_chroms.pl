#!/usr/bin/perl

# adds a cells to the database

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::GALFile;
use PSRG::Utils;
use PSRG::dataset;
use Bio::SeqIO;
use Bio::SearchIO;

my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
my $passwd = `cat ~/.oracle_passwd`;
chomp($passwd);
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
$dbh->do("alter session set current_schema=arolfe");
my ($species,$version,$design,$blastname);
GetOptions("species=s"=>\$species,
	   "genomeversion=s"=>\$version,
	   "design=s"=>\$design,
	   "blast=s"=>\$blastname);

my %hits = ();
my $hashname = $blastname . '.tophits';
if (-e $hashname) {
  print STDERR "Parsing tophits\n";
  my $h = PSRG::Utils::readHash(file=>$hashname);
  %hits = %$h;
  $h = undef;
} else {  
  my $blast = new Bio::SearchIO(-format=>'blast',
				-file=>$blastname);
  while( my $result = $blast->next_result ) {
    my $probename = $result->query_name() . $result->query_description();
    my $hit = $result->next_hit();
    my $hsp = $hit->next_hsp();
    if ($hits{$probename}) {
      warn "Duplicate hit for $probename";
      if ($hits{$probename} ne $hit->name() . ':' . $hsp->start('hit') . '-' . $hsp->end('hit')) {
	die "Hit Moves : $probename";
      }
    }
    $hits{$probename} = $hit->name() . ':' . $hsp->start('hit') . '-' . $hsp->end('hit');
#    print STDERR "Got Name as $probename\n";
  }
  PSRG::Utils::saveHash(file=>$hashname,hash=>\%hits);  
}

my @results = $dbh->selectrow_array("select id from species where name = '$species'");
unless (@results) {
  die "No such species : $species";
}
my $speciesid = $results[0];
@results = $dbh->selectrow_array("select id from genome where species = $speciesid and version = '$version'");
unless (@results) {
  die "No such genome : $species, $version";
}
my $genomeid = $results[0];
@results = $dbh->selectrow_array("select id from arraydesign where name='$design'");
if (@results == 0) {
  die "Can't find $design";
}
my $designid = $results[0];
@results = $dbh->selectrow_array("select id from chromosome where genome = $genomeid and name = 'X'");
if (@results == 0) {
  die "Can't find chromosome X";
}
my $Xid = $results[0];
@results = $dbh->selectrow_array("select id from chromosome where genome = $genomeid and name = '10'");
if (@results == 0) {
  die "Can't find chromosome 10";
}
my $tenid = $results[0];
my $check = $dbh->prepare("select count(*) from probedesign where arraydesign = $designid and chromosome = $tenid and probeid = ?");
my $fix = $dbh->prepare("update probedesign set chromosome = $Xid where arraydesign = $designid and chromosome = $tenid and probeid = ?");
foreach my $probe (keys %hits) {
  next unless ($hits{$probe} =~ /^chrX/);
#  print STDERR "Fixing $probe : $hits{$probe}\n";
#  $check->execute($probe);
#  @results = $check->fetchrow_array();
  #  if ($results[0] >= 1) {
    $fix->execute($probe);
#  } 
}

