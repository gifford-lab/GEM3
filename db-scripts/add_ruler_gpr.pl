#!/usr/bin/perl

# cy3 = green
# cy5 = red

use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Utils;
use PSRG::AFEFile;
use PSRG::AddRulerAFEHandler;

use POSIX;
use PSRG::Database;

my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhruler = PSRG::Database::handleForRole('rulers');

my ($rexpt, $rversion, $rreplicate, $rfragdistname, $rfragdistversion, $rspecies, $rcells, $rcondition);
my ($gexpt, $gversion, $greplicate, $gfragdistname, $gfragdistversion, $gspecies, $gcells, $gcondition);
my ($designname, $nofar, $blockoffsetparam);
GetOptions("designname=s"=>\$designname,
	   "nofar"=>\$nofar,
	   "rexpt=s"=>\$rexpt,
	   "rfragdist=s"=>\$rfragdistname,
	   "rspecies=s"=>\$rspecies,
	   "rcells=s"=>\$rcells,
	   "rcondition=s"=>\$rcondition,
	   "gexpt=s"=>\$gexpt,
	   "gfragdist=s"=>\$gfragdistname,
	   "gspecies=s"=>\$gspecies,
	   "gcells=s"=>\$gcells,
	   "gcondition=s"=>\$gcondition,
	   "blockoffset=s"=>\$blockoffsetparam);

# in case we got a genome version too
$rspecies =~ s/;.*//;
$gspecies =~ s/;.*//;

my $rspeciesid = PSRG::Database::getSpeciesID($dbhcore,$rspecies);
my $gspeciesid = PSRG::Database::getSpeciesID($dbhcore,$gspecies);
my @results = $dbhruler->selectrow_array("select id from arraydesign where name = '$designname'");
die "Can't find design $designname" if (@results == 0);
my $designid = $results[0];


($rexpt,$rversion,$rreplicate) = split(';',$rexpt);
($gexpt,$gversion,$greplicate) = split(';',$gexpt);
($rfragdistname,$rfragdistversion) = split(';',$rfragdistname);
($gfragdistname,$gfragdistversion) = split(';',$gfragdistname);
my $rexptid = get_expt($rexpt,$rspeciesid,$rversion,$rreplicate,$rfragdistname,$rfragdistversion,
		      $rcells,$rcondition,'');
my $gexptid = get_expt($gexpt,$gspeciesid,$gversion,$greplicate,$gfragdistname,$gfragdistversion,
		      $gcells,$gcondition,'');
print STDERR "Experiments are $rexptid, $gexptid\n";
$dbhruler->{AutoCommit} = 0;
$dbhruler->{RaiseError} = 1;

my $stmt = $dbhruler->prepare("insert into datatemp(experiment,probeid,blockno,rowno,colno,medianval,bgsubmedian) values (?,?,?,?,?,?,?)");

foreach my $fname (@ARGV) {
  my ($gprfile,$handler);
  print STDERR "ADDING FILE $fname\n";
  $dbhruler->do("delete from datatemp");
  if ($fname =~ /afe$/ or $fname =~ /txt$/) {
    $handler = new PSRG::AddRulerAFEHandler(insert=>$stmt,
					    redexpt=>$rexptid,
					    greenexpt=>$gexptid,
					    dbh=>$dbhruler);
    $gprfile = new PSRG::AFEFile($fname,$handler);
  } else {
    die "Unknown type for $fname";
  }
  my ($maxcol,$maxrow,$blockrow)=(0,0,0);
  my $galfileid = $handler->{galfileid} or die "No galfileid";
  ($maxcol,$maxrow) = $dbhruler->selectrow_array("select max(colno), max(rowno) from probedesign where galfile = $galfileid and arraydesign = $designid");
  my $blockoffset = $blockoffsetparam || $maxrow;
  my $sth = $dbhruler->prepare("select id, blockno, colno, rowno, probeid from probedesign where arraydesign = $designid and galfile = $galfileid");
  print STDERR "Populating probecache for $galfileid\n";
  $sth->execute();
  my @probecache = ();
  while (@results = $sth->fetchrow_array()) {
    my ($id,$block,$col,$row,$probeid) = @results;
    $probecache[$block][$col][$row] = {probeid=>$probeid,
				       id=>$id};
  }
  my ($afe,$far, $half, $side);
  ($afe,$far, $half, $side, $blockoffsetparam) = PSRG::ArrayDesign::determine_geometry($dbhruler,\@probecache,$maxrow,$maxcol,$blockoffsetparam);
  PSRG::ArrayDesign::fix_geometry($dbhruler,$galfileid,$maxrow,$maxcol,$blockoffset,$afe,$far,$half,$side,$blockoffsetparam);
  $dbhruler->do("update datatemp set id = (select id from probedesign pd where arraydesign = $designid and " .
		"galfile = $galfileid and pd.blockno = datatemp.blockno and pd.colno = datatemp.colno " .
		"and pd.rowno = datatemp.rowno and pd.probeid = datatemp.probeid)");
  @results = $dbhruler->selectrow_array("select count(*) from datatemp");
  print STDERR "There are $results[0] rows\n";
  @results = $dbhruler->selectrow_array("select count(*) from datatemp where id is null");
  print STDERR "There were $results[0] probes that I will delete because they aren't in the probedesign table\n";
  $dbhruler->do("delete from datatemp where id is null");
  $dbhruler->do("insert into data(experiment,probe,medianval, bgsubmedian) " .
	       "(select experiment, id, medianval, bgsubmedian from datatemp)");
  $dbhruler->do("delete from datatemp");
  $dbhruler->commit();
  print STDERR "Finished loading $fname\n";
}
$dbhruler->commit();
$dbhruler->disconnect();





sub get_expt {
  my ($expt,$speciesid,$version,$replicate,$fragdistname,$fragdistversion,$cells,$condition,$normalization) = @_;
  my @results = $dbhruler->selectrow_array("select id from experiment where name = '$expt' and species = $speciesid and version = '$version' and replicate = '$replicate'");
  if (@results) {
    print STDERR "Got id for $expt, $speciesid as $results[0]\n";
    return $results[0];
  } else {
    @results = $dbhruler->selectrow_array("select id from fragdist where name='$fragdistname' and version='$fragdistversion'");
    if (@results == 0) {
      die "Can't find fragdist $fragdistname, $fragdistversion";
    }
    my $fragdistid = $results[0];
    
    # Make sure both cells types are there and get there IDs
    @results= $dbhcore->selectrow_array("select id from cells where name='$cells'");
    if (@results == 0) {
      die "Can't find cells $cells";
    }
    my $cellsid = $results[0];
    # get the condition IDs
    @results= $dbhcore->selectrow_array("select id from condition where name='$condition'");
    if (@results == 0) {
      die "Can't find condition $condition";
    }
    my $conditionid = $results[0];
    my @fields = qw( id name version replicate fragdist species cells condition normalization active );
    my @values = (PSRG::Database::AIinsertText($dbhruler,'experiment_id'),"'$expt'","'$version'","'$replicate'",$fragdistid,
		  $speciesid,$cellsid, $conditionid,
		  "'$normalization'",1);
    $dbhruler->do("insert into experiment (".join(",",@fields).") values (".join(",",@values).")");
    @results = $dbhruler->selectrow_array("select ".PSRG::Database::AIfetchValue($dbhruler,'experiment_id'));
    print "Created experiment id $results[0] for $expt\n";
    return $results[0];
  }
}
  
