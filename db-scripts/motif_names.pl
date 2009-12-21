#!/usr/bin/perl

# fills the gff_motif_name table

use strict;
use warnings;
use DBI;
use PSRG::Utils;
use URI::Escape;
use Getopt::Long;

my $dbtype = $ENV{PSRGDBTYPE} || 'oracle';
my $dbh;
if ($dbtype eq 'oracle') {
  my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
  my $passwd = `cat ~/.oracle_passwd`;
  chomp($passwd);
  $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
  $dbh->do("alter session set current_schema=arolfe");
  
} elsif ($dbtype eq 'mysql') {
  my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
  if (`hostname` =~ /opteron.csail.mit.edu/) {
    $host = 'localhost';
  }
  my $passwd = `cat ~/.mysql_passwd`;
  chomp($passwd);
  $dbh = DBI->connect("DBI:mysql:host=${host};database=${sid}", $user, $passwd) or
    die "Can't connect as $user/$passwd: $!";
} else {
  die "Unknown database type $dbtype";
}

my ($fromscratch,$maxid,$table);
GetOptions("fromscratch"=>\$fromscratch);

if ($fromscratch) {
  $maxid = 0;
  $dbh->do("drop table old_gff_motif_name");
  $dbh->{RaiseError} = 1;
  $dbh->do("create table new_gff_motif_name(id number(11), name varchar2(100))");
  $dbh->do("grant select on new_gff_motif_name to psrgread");
  $dbh->do("grant all on new_gff_motif_name to psrgread");
  $table = 'new_gff_motif_name';
} else {
  my @results =  $dbh->selectrow_array("select max(id) from gff_motif_name");
  $maxid = $results[0] || -1;
  $table = 'gff_motif_name';
}

my $sth = $dbh->prepare("select id, attributes from gff where type = 'Motif' and id > $maxid");
$sth->execute();
my $insert = $dbh->prepare("insert into $table (id, name) values( ?,?)"); 
my $i = 0;
$dbh->{AutoCommit} = 0;
$dbh->{RaiseError} = 1;
while (my @results = $sth->fetchrow_array) {
  my ($id,$attr) = @results;
  my @pairs = split(';',$attr);
  my %attrs;
  my $name;
  foreach (@pairs) {
    my ($k,$v) = split('=',$_);
    my @v = map {uri_unescape($_)} split(',',$v);
    $attrs{$k} = \@v;
  }
  if ($attrs{Site}) {
    $insert->execute($id,$attrs{Site}[0]);
  }
  if ($i % 1000 == 0) {
    $dbh->commit();
  }
}
if ($fromscratch) {
  $dbh->do("drop index ix_new_motif_genename");
  $dbh->do("create index ix_new_motif_genename on new_gff_motif_name(name)");
  $dbh->do("rename gff_motif_name to old_gff_motif_name");
  $dbh->do("rename new_gff_motif_name to gff_motif_name");
}
$dbh->commit();
$dbh->disconnect();
