#!/usr/bin/perl

use strict;
use warnings;

use DBI;
use Getopt::Long;
use PSRG;
use PSRG::GALFile;
use PSRG::Utils;
use Bio::SeqIO;
use PSRG::Database;

my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhrulers = PSRG::Database::handleForRole('rulers');

my ($fastaname,$splitpslfname, $design, $species, $genomeversion,$allnew);
GetOptions("species=s"=>\$species,
	   "genomeversion=s"=>\$genomeversion,
	   "design=s"=>\$design,
	   "fasta=s"=>\$fastaname,
	   "psl=s"=>\$splitpslfname);

if (not $genomeversion) {
  ($species,$genomeversion) = split(';',$species);
}
unless (-e $fastaname) {
  die "No such file $fastaname for FASTA file";
}
unless (-e $splitpslfname) {
  die "No such file $splitpslfname for split probe PSL output";
}
my @results = $dbhcore->selectrow_array("select id from species where name = '$species'");
unless (@results) {
  die "No such species : $species";
}
my $speciesid = $results[0];

@results = $dbhcore->selectrow_array("select id from genome where species = $speciesid and version = '$genomeversion'");
unless (@results) {
  die "No such genome : $species, $genomeversion";
}
my $genomeid = $results[0];

@results = $dbhrulers->selectrow_array("select count(*) from arraydesign where name='$design'");
if ($results[0] == 0) {
  die "Can't find design $design";
} 
@results = $dbhrulers->selectrow_array("select id from arraydesign where name='$design'");
my $designid = $results[0];

my %ids = ();
my $sth = $dbhrulers->prepare("select probeid, id from probedesign where arraydesign = $designid");
$sth->execute();
while (my @r = $sth->fetchrow_array()) {
  if (exists $ids{$r[0]}) {
    die "@r : $ids{$r[0]}";
  }
  $ids{$r[0]} = $r[1];
}

print STDERR "Fixing probe sequences\n";
$sth = $dbhrulers->prepare("update probedesign set sequenceleft = ?, sequenceright = ? where id = ?");

# my $fasta = new Bio::SeqIO(-file=>$fastaname,
# 			   -format=>'fasta');
# while (my $seq = $fasta->next_seq()) {
#   my $probeseq = $seq->seq();
#   my $name = $seq->id();
#   if ($ids{$name}) {
#     $sth->execute(substr($probeseq,0,30),
# 		  substr($probeseq,30,30),
# 		  $ids{$name});
#   }
# }
print STDERR "Done fixing probe sequences\n";

print STDERR "Clearing probelocation and probeinterval\n";
$dbhrulers->do("delete from probelocation");
$dbhrulers->do("delete from probeinterval");

my %chromcache = ();
my $locstmt = $dbhrulers->prepare("insert into probelocation (id, chromosome, startpos, stoppos, probeside) values (?,?,?,?,?)");
my $intstmt = $dbhrulers->prepare("insert into probeinterval (id, chromosome, startpos, stoppos, orientation, leftstart, leftstop, rightstart, rightstop) " .
				  "values (?,?,?,?,?,?,?,?,?)");

print STDERR "Adding positions and intervals\n";
open(PSL,$splitpslfname) or die "Can't open $splitpslfname : $!";
do {
  $_ = <PSL>;
} until (/\-{30}/);

my %lefts = ();
my %rights = ();
my $lastprobe = '';
my $lastchrom;
while (<PSL>) {
  chomp;
  my @line = split(/\t/,$_);
  my ($match,$mismatch,$gapcount,$strand,$qname,$qstart,$qend,$tname,$tstart,$tend) = @line[0,1,6,8,9,11,12,13,15,16];
  next unless (($match - $mismatch) >= 25 and $gapcount == 0);
  my ($probeid,$type) = ($qname =~ /(.*)_(.*)/);
  
  if ($lastprobe and $probeid ne $lastprobe) {
    print STDERR "Adding intervals for $lastprobe: $ids{$lastprobe}\n";
    foreach my $c (keys %lefts) {
      foreach my $leftside (@{$lefts{$c}}) {
	foreach my $rightside (@{$rights{$c}}) {
	  next unless (abs($leftside->{tstart} - $rightside->{tstart}) < 100000);
	  next unless ($leftside->{strand} eq $rightside->{strand});
	  my @p = sort {$a <=> $b} ($leftside->{tstart}, $leftside->{tend},
				    $rightside->{tstart}, $rightside->{tend});
	  print STDERR "\t$c: @p\n";
	  my @v = ($c,
		   $p[1],$p[2],
		   $leftside->{tstart} > $rightside->{tstart} ? '-' : '+',
		   $leftside->{qstart}, $leftside->{qend},
		   $rightside->{qstart}, $rightside->{qend});
	  my $pid = $ids{$lastprobe};
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
    %lefts = ();
    %rights = ();
  }


  my $chrom = $tname;
  my $chromid;
  $chrom =~ s/^chr0?//;
  $chrom =~ s/\.fs?a//;
  if ($chrom eq 'Mito') {
    $chrom = 'mt';
  }
  
  if ($chromcache{$chrom}) {
    $chromid = $chromcache{$chrom};
  } else {
    @results = $dbhcore->selectrow_array("select id from chromosome where genome = $genomeid and name = '$chrom'");
    $chromid = $results[0];
    $chromcache{$chrom} = $chromid;
    unless ($chromid) {
      warn "Can't find a chromosome for $genomeid, $lastprobe: $chrom";
      next;
    }
  }
  unless ($chromid) {
    warn "Can't get a chromosome for $chrom : $chromid";
  }

  my $pid = $ids{$probeid};
  my @p = ($pid,$chromid,$tstart,$tend,$type eq 'left' ? 'l' : 'r');
  eval {
    $locstmt->execute(@p);
  };
  if ($@ and $@ !~ /duplicate/i and $@ !~ /ORA.00001/) {
    $dbhrulers->rollback();
    die "$@";
  }
  if ($type eq 'left') {
    push(@{$lefts{$chromid}},{score=>$match,
			    strand=>$strand,
			    qstart=>$qstart,
			    qend=>$qend,
			    tend=>$tend,
			    tstart=>$tstart});
  } else {
    push(@{$rights{$chromid}},{score=>$match,
			     strand=>$strand,
			     qstart=>$qstart,
			     qend=>$qend,
			     tend=>$tend,
			     tstart=>$tstart});
  }
  $lastprobe = $probeid;
  $lastchrom = $chromid;
}

print STDERR "Added probe locations\n";
$dbhrulers->commit();
$dbhrulers->disconnect();
