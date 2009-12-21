#!/usr/bin/perl

# generates .data files for the specified species and experiment.  These are the input files for JBD
#
# --expt : experiment name;version;replicate.  can be given more than once
# --minstd : each output row is one probe, however, the standard deviation of that probe across all
#            experiments being dumped is computed.  If minstd is supplied, the standard deviation
#            in the output files is never less than this value
# --bychrom : generate a directory per chromosome
# --subdirs : break output up over multiple subdirectories


use strict;
use warnings;
use DBI;
use Getopt::Long;
use PSRG;
use PSRG::Utils;
use PSRG::GPRFile;
use Statistics::Lite;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');

my (@expt, $species, $minstd, $subdirs, $bychrom, $genomeversion, $maxdist,$outname, $design,$coeffsonly,$randomizeprobes,$randomizevalues,$unique);
$minstd = 0;
$maxdist = 0;
my $gff= 0;
$coeffsonly = 0;
$randomizeprobes = 0;
$randomizevalues = 0;
GetOptions("expt=s"=>\@expt,  
	   "species=s"=>\$species,
	   "minstd=s"=>\$minstd,
	   "bychrom"=>\$bychrom,
	   "subdirs=s"=>\$subdirs,
	   "randomizevalues"=>\$randomizevalues,
	   "randomizeprobes"=>\$randomizeprobes,
	   "gff"=>\$gff,
	   "coeffsonly"=>\$coeffsonly,
	   "outname=s"=>\$outname,
	   "maxdist=s"=>\$maxdist,
	   "unique=s"=>\$unique,  # maximum genomic matches for probes to be included in the output
	   "design=s"=>\$design);

# bychrom creates on subdirectory per chromosome
# subdirs sets the number of subdirectories across which
#   the output files should be spread.  It can be combined 
#   with bychrom to subdivide each chromosome.
# If either bychrom or subdirs is set, then the d ata will be split into pieces
#   based on the fragment length distribution size and the gaps between probes. 
#   If neither is set, then there will be exactly one file per chromosome.

unless ($outname) {
  die "Must supply an output filename base with --outname";
}


($species,$genomeversion) = split(';',$species);
# get the experiment ID
my ($speciesid);
my @fragdistid;
my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
if (@results == 0) {
  die "Can't find species $species";
}
$speciesid = $results[0];
if ($design) {
  @results = $dbhchip->selectrow_array("select id from arraydesign where name = '$design'");
  $design = $results[0];
}
my %exptsseen = ();
my @exptid=();
my $repcount = 1;
my %repmap = ();
my @repnames;
foreach my $expt (@expt) {
  my ($exptname,$exptversion,$replicate) = split(';',$expt);
  if ($replicate) {
    @results = $dbhchip->selectrow_array("select id, fragdist from experiment where name = '$exptname' and species = $speciesid and version = '$exptversion' and replicate = '$replicate'");
    if (@results) {
      push(@exptid,$results[0]);
      $repmap{$results[0]} = $repcount;
      push(@fragdistid,$results[1]);
      push(@repnames,$repcount++);
    } else {
      die "Can't get experiment with name $exptname and version $exptversion";
    }
  } else {
    my $sth = $dbhchip->prepare("select id, fragdist, replicate from experiment where name = '$exptname' and species = $speciesid and version = '$exptversion' and active = 1");
    print STDERR "Executing \"select id, fragdist, replicate from experiment where name = '$exptname' and species = $speciesid and version = '$exptversion'\"\n";
    $sth->execute();
    while (my @results = $sth->fetchrow_array()) {
      push(@exptid,$results[0]);
      $repmap{$results[0]} = $repcount;
      push(@repnames,$repcount++);
      print STDERR "got $results[0] as $exptname, $exptversion\n";
      push(@fragdistid,$results[1]);
    }
  }
}
unless (@exptid) {
  die "Didn't find any experiments in the data";
}
my $exptid = "(" . join(" or ", map {' data.experiment = ? '} @exptid) . ")";

@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$genomeversion'");
if (@results == 0) {
  die "Can't get genome for $speciesid, $genomeversion";
}
my $genomeid = $results[0];
foreach my $fragdistid (@fragdistid) {
  unless ($maxdist) {
    @results = $dbhchip->selectrow_array("select max(distance) from fragdistentry where value > 0 and distribution = $fragdistid");
    if (@results == 0) {
      die "Can't get maxdist";
    }
    if ($results[0] > $maxdist) {
      $maxdist = $results[0];
    }
  }
}

my %chroms = (); my %chromids = ();
my $sth = $dbhcore->prepare("select id, name from chromosome where genome = $genomeid");
$sth->execute();
my @chroms = ();
while (my @row = $sth->fetchrow_array()) {
  $chroms{$row[1]} = $row[0];
  $chromids{$row[0]} = $row[1];
#  print STDERR "CHROM $row[0] -> $row[1]\n";
  push(@chroms,$row[0]);
}
my $chromclause = join(',',@chroms);
my $safeexptname = $outname;
$safeexptname =~ s/\s/_/g;

my ($piece,$fname,$dir,$lastpos,$chrom);
my @args = ("--species '${species};${genomeversion}' --pattern '${outname}'");  # args for add_bayesian or add_mle
my @d;
my %directories = ();
my $dumpedthisfile = 0;
unless ($coeffsonly) {
  # this comes up with the list of output directories to use.  keyed on chromosome, 
  # each entry can either be a string (whole chrom in one dir) or an array
  my %outdirs;
  system("mkdir -p '$outname'");
  chdir($outname);
  if ($bychrom) {
    if ($subdirs) {
      foreach my $chrom (keys %chroms) {
	my @d = map {"./${chrom}/${_}"} 1..$subdirs;   
	$outdirs{$chrom} = \@d;
	foreach (@d) {
	  system("mkdir -p '$_'");
	  $directories{$_} = 1;
	}    
      }
    } else {
      foreach (keys %chroms) {
	$outdirs{$_} = "./${_}";
	system("mkdir -p './${_}'");
	$directories{"./$_"} = 1;
      }
    }
  } else {
    if ($subdirs) {
      my @d = map {"./${_}"} 1..$subdirs;   
      foreach my $chrom (keys %chroms) {
	$outdirs{$chrom} = \@d;
      }
      foreach (@d) {
	system("mkdir -p '$_'");
	$directories{$_} = 1;
      }
    } else {
      foreach (keys %chroms) {
	$outdirs{$_} = '.';
	$directories{'.'} = 1;
      }
    }
  }
  foreach (keys %directories) {
    push(@args,"--dir '$_'");
  }
  my (@vals, %vals, $randprobeid);
  if ($randomizevalues or $randomizeprobes) {
    if ($design) {
      $sth = $dbhchip->prepare("select data.probe, data.channelone, data.channeltwo, data.ratio,data.experiment " .
			       " from data, probedesign where $exptid and " .
			       " and data.probe = probedesign.id and probedesign.arraydesign = $design order by data.probe");
      
    } else {
      $sth = $dbhchip->prepare("select data.probe,data.channelone, data.channeltwo, data.ratio,data.experiment " .
			       " from data where $exptid order by data.probe");
      
    }
    $sth->execute(@exptid);
    while (my ($probe,$ip,$wce,$ratio,$rep) = $sth->fetchrow_array()) {
      if ($randomizevalues) {
	push(@vals,[rand(),$ip,$wce,$ratio]);
      } else {
	push(@{$vals{$probe}},[$ip,$wce,$ratio,$repmap{$rep}]);
      }
    }
    if ($randomizevalues) {
      @vals = sort {$a->[0] <=> $b->[0]} @vals;    
    } else {
      @vals = sort {rand() <=> rand()} keys %vals;
    }
    print STDERR "Have $#vals randomized datapoints\n";
  }
  if ($design) {
    $sth = $dbhchip->prepare("select probelocation.chromosome, data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, data.experiment " .
			     " from data, probelocation, probedesign where $exptid and data.probe = probelocation.id " .
			     " and data.probe = probedesign.id and probedesign.arraydesign = $design " . 
			     ($unique ? " and probelocation.loccount <= $unique" : '') .
 			     " order by probelocation.chromosome, probelocation.startpos");
  
  } else {
    $sth = $dbhchip->prepare("select /*+ LEADING (data)*/ /*+ INDEX (ix_data_expt)*/ probelocation.chromosome, data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, data.experiment " .
			     " from data, probelocation where $exptid and data.probe = probelocation.id  " .
			     ($unique ? " and probelocation.loccount <= $unique" : '') . 
			     " order by probelocation.chromosome, probelocation.startpos");
  
  }

  my $lastchrom = '';
  $sth->execute(@exptid);
  while (my ($chromid,$ip,$wce,$ratio,$start,$stop,$rep) = $sth->fetchrow_array()) {
    if ($randomizevalues) {
      my $rv = shift(@vals);
      $ip = $rv->[1];
      $wce = $rv->[2];
      $ratio = $rv->[3];
    }

    $exptsseen{$rep} = 1;
    $chrom = $chromids{$chromid} || next;
    $rep = $repmap{$rep};
    if ($chrom ne $lastchrom) {
      $lastchrom = $chrom;
      $piece = 1;
      if (ref($outdirs{$chrom}) eq 'ARRAY') {
	my $i = $piece % @{$outdirs{$chrom}};
	$dir = $outdirs{$chrom}[$i];
      } else {
	$dir = $outdirs{$chrom};
      }
      if ($gff) {
	$fname = "${dir}/${outname}.${chrom}.${piece}.gff";
      } else {
	$fname = "${dir}/${outname}.${chrom}.${piece}.data";
      }
      open(DATA,">$fname") or die "Can't open $fname to write : $!";
      $dumpedthisfile = 0;
      print STDERR "writing $fname\n";    
      $lastpos = -1;
    }
    my $pos = int(($start + $stop)/2);
    if ($pos != $lastpos and $lastpos != -1) {
      if ($randomizeprobes) {
	$randprobeid = shift @vals;
	@d = ();
	foreach (@{$vals{$randprobeid}}) {
	  push(@d,$_);
	}
      }

      dumpdata(@d);
      @d = ();
      if ($bychrom || $subdirs and ($lastpos != -1) and ($pos > $lastpos + 2*$maxdist) and ($dumpedthisfile > 20000)) {
	close DATA;
	$piece++;
	if (ref($outdirs{$chrom}) eq 'ARRAY') {
	  my $i = $piece % @{$outdirs{$chrom}};
	  $dir = $outdirs{$chrom}[$i];
	} else {
	  $dir = $outdirs{$chrom};
	}
	if ($gff) {
	  $fname = "${dir}/${outname}.${chrom}.${piece}.gff";
	} else {
	  $fname = "${dir}/${outname}.${chrom}.${piece}.data";
	}
	open(DATA,">$fname") or die "Can't open $fname to write : $!";
	$dumpedthisfile = 0;
      }
    }
    push(@d,[$ip,$wce,$ratio,$rep]);
    $lastpos = $pos;
  }
  if (@d) {			# process remaining datapoints
    dumpdata(@d);
    @d = ();
  }
  close DATA;  
}
		       
for (my $i = 0; $i < @fragdistid; $i++) {
  my $fragdistid = $fragdistid[$i];
  $sth = $dbhchip->prepare("select distance, value from fragdistentry where distribution = $fragdistid order by distance");
  $sth->execute();
  my $repname = $repnames[$i];
  my $fname = "./${outname}.${repname}.coeffs";
  open(COEFFS,">$fname") || die "can't open $fname for writing : $!";
  while (my ($dist,$val) = $sth->fetchrow_array()) {
    print COEFFS "$dist\t$val\n";
  }
  close COEFFS;
}

# this section generates partial input arguments to add_bayesian or add_mle
$sth = $dbhchip->prepare('select name, version, replicate from experiment where id = ?');
foreach (keys %exptsseen) {
  $sth->execute($_);
  my ($name,$version,$replicate) = $sth->fetchrow_array();
  print STDERR "$_ -> $name, $version, $replicate\n";
  push(@args,(" --expt '$name;$version;$replicate'"));
}
my $argsfname = "./args.txt";
open(ARGS,">$argsfname") or die "Can't open $argsfname : $!";
print ARGS join(" ",@args) . "\n";
close ARGS;

sub dumpdata {
  my @ratios = map {$_->[2]} @_;
  my $std = Statistics::Lite::stddev(@ratios) + $minstd;
  #      if ($gff) {@d = ($d[0])};
  foreach (@_) {
    $dumpedthisfile++;
    if ($gff) {
      print DATA join("\t",('chr' . $chrom,
			    $safeexptname,
			    'agilentobs',
			    $lastpos-30,
			    $lastpos+30,
			    log($_->[2] > .1 ? $_->[2] : .1),
			    '.',
			    '.',
			    'obs'.$lastpos))."\n";
      
    } else {
      print DATA "$lastpos $_->[0] $_->[1] $_->[2] $std $_->[3]\n";
    }
  }
}
