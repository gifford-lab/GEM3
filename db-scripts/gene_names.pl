#!/usr/bin/perl

# fills the gff_gene_name table by parsing the name from the 
# attribues.
# Note that this must use the same attribute parsing as edu.mit.csail.psrg.Bio.RecordTranslator
# 

use strict;
use warnings;
use DBI;
use PSRG::Utils;
use URI::Escape;

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
my $sth = $dbh->prepare("select id, attributes from gff where type = 'gene'");
$sth->execute();

my $insert = $dbh->prepare("insert into gff_gene_name (id, name) values( ?,?)"); 
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
  foreach my $field (qw(gene ID Alias Entrez Name)) {
    if ($attrs{$field}) {
      foreach (@{$attrs{$field}}) {
	$insert->execute($id,$_);
      }
    } 
  }
}
