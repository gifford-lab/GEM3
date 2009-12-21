#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;

my ($host,$sid,$user,$passwd,$table,$file,$species,$version);
$host = 'opteron.csail.mit.edu';
$sid = 'psrg';
GetOptions("host=s"=>\$host,
	   "sid=s"=>\$sid,
	   "user=s"=>\$user,
	   "passwd=s"=>\$passwd,
	   "table=s"=>\$table,
	   "file=s"=>\$file,	   
	   "species=s"=>\$species);
unless ($table and -r $file) {
  die "Must supply --table and --file";
}
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);

# step one is to learn about the table;
print "load data infile '$file' \"str '\\n'\"\n";
print "into table $table fields terminated by '\\t'\n";
print "trailing nullcols\n";
my $sth = $dbh->column_info(undef,undef,uc($table),undef);
my @cols = ();
my $chromcol = -1;
my $colnum = 0;
while (my @row = $sth->fetchrow_array()) {
   my $type = 'char';
   if ($row[5] =~ /(NUMBER)|(DOUBLE)|(FLOAT)/i) {
     if ($row[8] == 0) {
       $type = "integer external nullif $row[3]='NULL'";
     } else {
       $type = "float external nullif $row[3]='NULL'";
     }
   }
   if ($row[5] eq 'VARCHAR2') {
     $type = "char ($row[6])";
   }
   push(@cols, "$row[3] $type");
   if (lc($row[3]) eq 'chrom') {
     $chromcol = $colnum;
   }
   $colnum++;
}
print "(" . join(",\n",@cols) . ")\n";

if ($chromcol) {
  ($species,$genomeVersion) = split(';',$species);
  my @results = $dbh->selectrow_array("select id from species where name = '$species'");
  if (@results == 0) {
    die "Can't find species $species";
  }
  $speciesid = $results[0];
  @results = $dbh->selectrow_array("select id from genome where species = $speciesid and version = '$genomeVersion'");
  if (@results == 0) {
    die "Can't find genome for $speciesid, $genomeVersion";
  }
  $genomeid = $results[0] or die "Can't get genomeid for $speciesid and $genomeVersion";

  my $sth = $dbh->prepare("select id, name from chromosome where genome = $genomeid");
  my %chrommap;
  while (my @r = $sth->fetchrow_array()) {
    $chrommap{'chr'. $r[1]} = $r[0];
  }
  my $tmpfile = $file . '.bak';
  system("mv '$file' '$tmpfile'");
  open(IN,$tmpfile) or die "Can't open $tmpfile : $!";
  open(OUT,">$file") or die "Can't open $file : $!";
  while (<IN>) {
    chomp;
    my @l = split('\t',$_);
    $l[$chromcol] = $chrommap{$l[$chromcol]};
    print OUT join("\t",@l) . "\n";
  }
}
