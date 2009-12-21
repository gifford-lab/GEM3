#!/usr/bin/perl

# adds data from a GPR file to the database
# subexpt is a colon-separated list of cells, condition, factor
# for the two conditions OR the control experiment id.
#
# either --exptid must be supplied or all of  designname, species, fragdistname, fragdistversion,
# cellsone, conditionone, factorone, cellstwo, conditiontwo, and factortwo
# must be supplied.
#
# --median: makes the median intensity of both channels equal
# --linefit : do linear regression on each GPR at then rotate the log-transformed data such
#     that it lies along the ip=wce line.
#
# --flippedchannels : normally assume that cy5 (635nm) is the IP channel and that cy3 (530nm) is WCE.
#                     this option flips the two channels.  It works *after* the file has been read but
#                     before any normalization happens

use strict;
use warnings;
use DBI;
use DBI::Const::GetInfoType;
use Getopt::Long;
use PSRG;
use PSRG::Utils;
use PSRG::GPRFile;
use PSRG::AFEFile;
use PSRG::PairFile;
use PSRG::CELFile;
use PSRG::AddGPRHandler;
use PSRG::AddAFEHandler;
use PSRG::AddPairHandler;
use PSRG::AddCELHandler;
use POSIX;
use PSRG::Database;
use Statistics::Regression;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');

my ($expt,$exptid,$designid,$sth,@results);
my ($designname,$species,$genomeversion,$cellsone,$conditionone,$factorone,
    $cellstwo, $conditiontwo,$factortwo,$normalization,$version);
my ($fragdistname,$fragdistversion,$median,$logratio, $replicate);
my ($mean,$linefit,$nofar,$noafe,$blockoffsetparam,$blockoffset,@celip, @celwce, @galoverride);
my $verbose = 0;
my $flippedchannels = 0;
GetOptions("expt=s"=>\$expt,   
	   "designname=s"=>\$designname,  
	   "species=s"=>\$species,
	   "fragdistname=s"=>\$fragdistname,
	   "fragdistversion=s"=>\$fragdistversion,
	   "cellsone=s"=>\$cellsone,
	   "conditionone=s"=>\$conditionone,
	   "factorone=s"=>\$factorone,
	   "cellstwo=s"=>\$cellstwo,
	   "conditiontwo=s"=>\$conditiontwo,
	   "factortwo=s"=>\$factortwo,
	   "celip=s"=>\@celip,
	   "celwce=s"=>\@celwce,
	   "gal=s"=>\@galoverride,
	   "normalization=s"=>\$normalization,
	   "median"=>\$median,
	   "linefit"=>\$linefit, 
	   "logratio"=>\$logratio,
	   "nofar"=>\$nofar,
	   "noafe"=>\$noafe,
	   "verbose"=>\$verbose,
	   "flippedchannels"=>\$flippedchannels,
	   "blockoffset=s"=>\$blockoffsetparam);

if (not $genomeversion) {
  ($species,$genomeversion) = split(';',$species);
}
my $speciesid = PSRG::Database::getSpeciesID($dbhcore,$species);
@results = $dbhchip->selectrow_array("select id from arraydesign where name = '$designname'");
die "Can't find design $designname" unless (@results);
$designid = $results[0];

($expt,$version,$replicate) = split(';',$expt);
my ($datatemptable, $galtemptable, $datatempindex, $galtempindex) = create_temp_tables($dbhchip);
$dbhchip->do("lock table ${datatemptable} in exclusive mode");
$dbhchip->do("lock table ${galtemptable} in exclusive mode");


# Get the experiment id, or create it if it doesn't exist yet
@results = $dbhchip->selectrow_array("select id from experiment where name = '$expt' and species = $speciesid and version = '$version' and replicate = '$replicate'");
if (@results) {
  $exptid = $results[0];
  print STDERR "Got id for $expt, $species as $exptid\n";
} else {
  @results = $dbhchip->selectrow_array("select id from fragdist where name='$fragdistname' and version='$fragdistversion'");
  if (@results == 0) {
    die "Can't find fragdist $fragdistname, $fragdistversion";
  }
  my $fragdistid = $results[0];
  
  # Make sure both cells types are there and get there IDs
  @results= $dbhcore->selectrow_array("select id from cells where name='$cellsone'");
  if (@results == 0) {
    die "Can't find cells $cellsone";
  }
  my $cellsoneid = $results[0];
  @results = $dbhcore->selectrow_array("select id from cells where name='$cellstwo'");
  if (@results == 0) {
    die "Can't find cells $cellstwo";
  }
  my $cellstwoid = $results[0];
  # get the condition IDs
  @results= $dbhcore->selectrow_array("select id from condition where name='$conditionone'");
  if (@results == 0) {
    die "Can't find condition $conditionone";
  }
  my $conditiononeid = $results[0];
  @results = $dbhcore->selectrow_array("select id from condition where name='$conditiontwo'");
  if (@results == 0) {
    die "Can't find condition $conditiontwo";
  }  
  my $conditiontwoid = $results[0];
  @results = $dbhcore->selectrow_array("select id from factors where name = '$factorone'");
  if (@results == 0) {
    die "Can't find factor $factorone";
  }
  my $factoroneid = $results[0];
  @results = $dbhcore->selectrow_array("select id from factors where name = '$factortwo'");
  if (@results == 0) {
    die "Can't find factor $factortwo";
  }
  my $factortwoid = $results[0];
  my @fields = qw(id name fragdist species cellsone conditionone factorone cellstwo conditiontwo
		  factortwo normalization version replicate active);
  my @values = (PSRG::Database::AIinsertText($dbhchip,'experiment_id'),"'$expt'",$fragdistid,$speciesid,$cellsoneid,
		$conditiononeid,$factoroneid,$cellstwoid,$conditiontwoid,$factortwoid,
 		"'$normalization'","'$version'","'$replicate'",1);
  $dbhchip->do("insert into experiment (".join(",",@fields).") values (".join(",",@values).")");
  @results = $dbhchip->selectrow_array("select ".PSRG::Database::AIfetchValue($dbhchip,'experiment_id'));
  $exptid = $results[0];
  print "Created experiment id $exptid for $expt\n";
}

# add mapping for this genome
my @genomes = ();
$sth = $dbhchip->prepare("select unique(chromosome) from probelocation where id in (select id from probedesign where arraydesign = $designid and rownum < 10000)");
my @chroms;
$sth->execute();
while (@results = $sth->fetchrow_array()) {
  push(@chroms,$results[0]);
}
$sth = $dbhcore->prepare("select unique(genome) from chromosome where id in (" . join(',',@chroms) . ")");
$sth->execute();
while (@results = $sth->fetchrow_array()) {
  push(@genomes,$results[0]);
}
foreach my $genomeid (@genomes) {
  @results = $dbhchip->selectrow_array("select count(*) from exptToGenome where experiment = $exptid and genome = $genomeid");
  if ($results[0] == 0) {
    $dbhchip->do("insert into exptToGenome(experiment,genome) values($exptid,$genomeid)");
  }
}


my ($stmt,$ipbs,$wcebs, $ipcol,$wcecol,$idcol,$blockcol,$rowcol,$colcol,$morcol);

$dbhchip->do("delete from ${datatemptable}");
$stmt = $dbhchip->prepare("insert /*+ append */ into ${datatemptable}(experiment,id,probeid,blockno,rowno,colno,channelone,channeltwo,channelratio,mor) values ($exptid,NULL,?,?,?,?,?,?,?,?)");
$dbhchip->{AutoCommit} = 0;
$dbhchip->{RaiseError} = 1;

my $usingcel = 0;
if (@celip) {
  die "Unequal ip and wce files" unless (@celip == @celwce);
  print STDERR "Using --celip and --celwce and ignoring other file names\n";
  @ARGV = @celip;
  $usingcel = 1;
}

foreach my $fname (@ARGV) {
  if ((not -e "/afs/csail.mit.edu/group/psrg/datasets/gpr_files/$fname") and 
      -e "/afs/csail.mit.edu/group/psrg/datasets/gpr_files/unsorted/$fname") {
    $fname = "unsorted/$fname";
  }

  print STDERR "ADDING FILE $fname\n";  
  my $prefix;
  if ($fname !~ /^\// and not -e $fname) {
    $prefix = '/afs/csail.mit.edu/group/psrg/datasets/gpr_files/';
  }
  print STDERR "Reading file to ${datatemptable} " . time() . "\n";
  my ($gprfile,$handler);
  if ($fname =~ /gpr$/i and not $usingcel) {
    $handler = new PSRG::AddGPRHandler(design=>$designid,
				       expt=>$exptid,
				       dbh=>$dbhchip,
				       insert=>$stmt);
    $gprfile = new PSRG::GPRFile($prefix.$fname,$handler);
  } elsif (($fname =~ /pair.txt/ or $fname =~ /\.pair$/)  and not $usingcel){
    $handler = new PSRG::AddPairHandler(design=>$designid,
					expt=>$exptid,
					dbh=>$dbhchip,
					insert=>$stmt);
    $gprfile = new PSRG::PairFile($prefix.$fname,$handler);
  } elsif (($fname =~ /afe$/i or $fname =~ /txt$/)  and not $usingcel){
    $handler = new PSRG::AddAFEHandler(design=>$designid,
				       expt=>$exptid,
				       dbh=>$dbhchip,
				       insert=>$stmt);   
    $gprfile = new PSRG::AFEFile($prefix.$fname,$handler);
  } elsif (($fname =~ /\.cel/i or $fname =~ /\.cel\.txt/i) and $usingcel) {
    $handler = new PSRG::AddCELHandler(design=>$designid,
				       expt=>$exptid,
				       dbh=>$dbhchip,
				       insert=>$stmt,
				       galfile=>shift(@galoverride));
    $gprfile = new PSRG::CELFile($prefix.$fname,$prefix.shift(@celwce),$handler);
  }
  my $galfileid = $handler->{galfileid} or die "No galfileid";

  if ($flippedchannels) {
    print STDERR "Flipping channels\n";
    $dbhchip->do("update ${datatemptable} set channelone = channeltwo, channeltwo = channelone, channelratio = 1.0/channelratio, mor = 1.0/mor");
  }

  $dbhchip->commit();
# all of the Add*Handlers now do the regularization when they read the data
#  print STDERR "Regularizing " . time() . "\n";
#  $dbhchip->do("update ${datatemptable} set channelone = greatest(1,channelone), channeltwo = greatest(1,channeltwo)");
#  $dbhchip->commit();
  print STDERR "Doing id assignments " . time() . "\n";
  if ($usingcel) {
    $dbhchip->{RaiseError} = 0;
    $dbhchip->do("drop index ${galtempindex}");
    $dbhchip->{RaiseError} = 1;
    $dbhchip->do("delete from ${galtemptable}");
    $dbhchip->do("insert /*+ append */ into ${galtemptable}(id,probeid,blockno,rowno,colno) (select id, probeid, blockno,rowno,colno from ".
		 "probedesign where arraydesign = $designid and galfile = $galfileid)");
    $dbhchip->commit();
    $dbhchip->do("create index ${galtempindex} on ${galtemptable}(blockno,rowno,colno) tablespace scratch nologging");
    $dbhchip->do("update ${datatemptable} set id = (select id from ${galtemptable} where " .
		 "${galtemptable}.blockno = ${datatemptable}.blockno and ${galtemptable}.colno = ${datatemptable}.colno " .
		 "and ${galtemptable}.rowno = ${datatemptable}.rowno)");
  } else {
    my ($maxcol,$maxrow,$blockrow)=(0,0,0);
    # lookup the probe
    $dbhchip->{RaiseError} = 0;
    $dbhchip->do("drop index ${galtempindex}");
    $dbhchip->{RaiseError} = 1;
    $dbhchip->do("delete from ${galtemptable}");
    $dbhchip->do("insert /*+ append */ into ${galtemptable}(id,probeid,blockno,rowno,colno) (select id, probeid, blockno,rowno,colno from ".
		 "probedesign where arraydesign = $designid and galfile = $galfileid)");
    $dbhchip->commit();
    my ($count) = $dbhchip->selectrow_array("select count(*) from ${galtemptable}");
    print STDERR "Populating probecache for $galfileid, $designid : count is $count\n";
    ($maxcol,$maxrow) = $dbhchip->selectrow_array("select max(colno), max(rowno) from ${galtemptable}");
    unless (defined $maxcol and defined $maxrow) {
      die "Couldn't get maxcol ($maxcol) and maxrow ($maxrow)";
    }
    $blockoffset = $blockoffsetparam || $maxrow;
    my $sth = $dbhchip->prepare("select id, blockno, colno, rowno, probeid from ${galtemptable}");
    $sth->execute();
    my @probecache = ();
    while (@results = $sth->fetchrow_array()) {
      my ($id,$block,$col,$row,$probeid) = @results;
      $probecache[$block][$col][$row] = {probeid=>$probeid,
					 id=>$id};
    }
    if ($gprfile->isa('PSRG::AFEFile') or $gprfile->isa('PSRG::GPRFile') or $gprfile->isa('PSRG::NDFFile')) {
      my ($afe,$far, $half, $side);
      ($afe,$far, $half, $side,$blockoffsetparam) = PSRG::ArrayDesign::determine_geometry($dbhchip,$datatemptable,\@probecache,$maxrow,$maxcol,$blockoffsetparam);
      PSRG::ArrayDesign::fix_geometry($dbhchip,$datatemptable,$galfileid,$maxrow,$maxcol,$blockoffset,$afe,$far,$half,$side,$blockoffsetparam);
    }
    $dbhchip->do("create index ${galtempindex} on ${galtemptable}(blockno,rowno,colno)  tablespace scratch nologging");
    $dbhchip->do("update ${datatemptable} set id = (select id from ${galtemptable} where " .
		 "${galtemptable}.blockno = ${datatemptable}.blockno and ${galtemptable}.colno = ${datatemptable}.colno " .
		 "and ${galtemptable}.rowno = ${datatemptable}.rowno and ${galtemptable}.probeid = ${datatemptable}.probeid)");
  }
  $dbhchip->commit();
  @results = $dbhchip->selectrow_array("select count(*) from ${datatemptable} where id is null");
  print STDERR "There were $results[0] probes that I will delete because they aren't in the probedesign table\n";
  $dbhchip->do("delete from ${datatemptable} where id is null");
  normalize();
  print STDERR "inserting into data " . time() . "\n";
  if ($logratio) {
    $dbhchip->do("insert /*+ append */ into data(experiment,probe,channelone,channeltwo,mor,channelratio,ratio) " .
		 "(select experiment, id, channelone, channeltwo, mor, channelratio, log(channelone/channeltwo) from ${datatemptable})");
  } else {
    $dbhchip->do("insert /*+ append */ into data(experiment,probe,channelone,channeltwo,mor,channelratio,ratio) " .
		 "(select experiment, id, channelone, channeltwo, mor, channelratio, (channelone/channeltwo) from ${datatemptable})");
  }

  $dbhchip->do("delete from ${datatemptable}");
  $dbhchip->commit();
  print STDERR "Finished loading $fname " . time() . "\n";
}
destroy_temp_tables($dbhchip,$datatemptable,$galtemptable);
$dbhchip->commit();
$dbhchip->disconnect();

sub create_temp_tables {
  my $dbh = shift(@_);
  my $schema = $dbh->{Username} || die "no username in $dbh";  
  print STDERR "Using schema $schema for temporary tables\n";
  my $createsql = <<DATATEMP;
  create table SCHEMA.datatempPID (
				   experiment number(10) not null,
				   id number(10),
				   probeid varchar(200),
				   blockno number(10),
				   rowno number(10),
				   colno number(10),
				   channelone binary_float,
				   channeltwo binary_float,
				   mor binary_float,
				   -- ROM from channelone / channeltwo
				   channelratio binary_float,
				   -- final output ratio
				   ratio binary_float,
				   controlratio binary_float) tablespace scratch nologging
DATATEMP
    $createsql =~ s/SCHEMA/$schema/;
  $createsql =~ s/PID/$$/;
  $dbh->do($createsql);
  $createsql = <<GALTEMP;
  create table SCHEMA.galtempPID (
			id number(10),
			probeid varchar2(200),
			arraydesign number(10),
			galfile number(10),
			blockno number(10),
			colno number(10),
			rowno number(10)
		       ) tablespace scratch nologging
GALTEMP
    $createsql =~ s/SCHEMA/$schema/;
  $createsql =~ s/PID/$$/;
  $dbh->do($createsql);
  return ("${schema}.datatemp$$","${schema}.galtemp$$",
	  "ix_datatemp$$","ix_galtemp$$");
  
}

sub destroy_temp_tables {
  my $dbh = shift(@_);
  my @tables = @_;
  foreach my $table (@tables) {
    $dbh->do("drop table $table");
  }
}

sub normalize {
  if ($median or $mean) {
    my ($ip,$wce) = (1,1);
    if ($median) {
      ($ip,$wce) = $dbhchip->selectrow_array("select median(channelone), median(channeltwo) from ${datatemptable}");
    } elsif ($mean) {
      ($ip,$wce) = $dbhchip->selectrow_array("select mean(channelone), mean(channeltwo) from ${datatemptable}");
    } 
    my $diff = ($wce+1)/($ip+1);
    $dbhchip->do("update ${datatemptable} set channelone = channelone * $diff");
  }
  if ($linefit) {
    $dbhchip->do("update ${datatemptable} set channelone = log(2,channelone), channeltwo = log(2,channeltwo)");
    my ($factor,$intercept);
    if ($dbhchip->get_info($GetInfoType{SQL_DBMS_NAME}) eq 'Oracle') {
      ($factor,$intercept) = $dbhchip->selectrow_array("select to_number(atan(REGR_SLOPE(channelone, channeltwo)) - atan(1)), to_number(REGR_INTERCEPT(channelone,channeltwo)) from ${datatemptable} " .
						       "where channelone - channeltwo < 3 and channelone - channeltwo > -3");
    } else {
      # we can compute the slope and intercept, but the update part will 
      # not work right.  If you do "update set c1 = c1*c2, c2 = c1*c2",
      # mysql updates c1 and then uses that value in the second clause.
      # oracle correctly reads both values and then does both computations.
      die "Sorry.  I can't do linefit on MySQL";
      my $tmpsth = $dbhchip->prepare("select channelone, channeltwo from ${datatemptable} where channelone - channeltwo < 3 and channelone - channeltwo > -3");
	$tmpsth->execute();
      my $reg = Statistics::Regression->new(2,"linefit",["const","channeltwo"]);
      while (my @r = $tmpsth->fetchrow_array()) {
	$reg->include($r[0],[1,$r[1]]);
      }
      my @theta = $reg->theta();
      ($factor,$intercept) = (POSIX::atan($theta[1]) - POSIX::atan(1),
			      $theta[0]);
    }
    print STDERR "FACTOR is $factor, intercept is $intercept\n";
    unless (defined($factor) and defined($intercept)) {
      print STDERR "Oops!\n";
      destroy_temp_tables($dbhchip,$datatemptable,$galtemptable);
      $dbhchip->commit();
      $dbhchip->disconnect();
      exit;      
    }
    $dbhchip->do("update ${datatemptable} set channelone = channelone - $intercept");    
    $dbhchip->do("update ${datatemptable} set channelone = power(2,sqrt(channelone*channelone + channeltwo*channeltwo) * sin(atan(channelone/channeltwo) - $factor)), " .
		 "channeltwo = power(2,sqrt(channelone*channelone + channeltwo*channeltwo) * cos(atan(channelone/channeltwo) - $factor))");
  }
}


#                if ($$header[$j] eq "rBGSubSignal")     {$$header[$j] = "F635 MEDIAN - B635";}
#                 if ($$header[$j] eq "gBGSubSignal")     {$$header[$j] = "F532 MEDIAN - B532";}
#                 if ($$header[$j] eq "SystematicName")   {$$header[$j] = "NAME";}
#                 if ($$header[$j] eq "ProbeName")                {$$header[$j] = "ID";}
#                 if ($$header[$j] eq "LogRatio")                 {$$header[$j] = "Median of Ratios";} ##wrong
#                 if ($$header[$j] eq "rBGPixSDev")               {$$header[$j] = "B635 SD";}
#                 if ($$header[$j] eq "gBGPixSDev")               {$$header[$j] = "B532 SD";}
#                 if ($$header[$j] eq "IsManualFlag")     {$$header[$j] = "FLAGS";}
  
