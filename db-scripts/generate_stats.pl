#!/usr/bin/perl

use strict;
use warnings;
use DBI;
my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
my $passwd = `cat ~/.oracle_passwd`;
chomp($passwd);
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);

my $sth = $dbh->prepare('select table_name from sys.user_tables');
$sth->execute();
my @tables = ();
while (my @r = $sth->fetchrow_array()) {
  push(@tables,$r[0]);
}
my @indexes = ();
$sth = $dbh->prepare('select index_name from sys.user_indexes');
$sth->execute();
while (my @r = $sth->fetchrow_array()) {
  push(@indexes,$r[0]);
}
foreach (@tables) {
  print "exec dbms_stats.gather_table_stats(  ownname=> '$ENV{USER}',  tabname=> '$_',  granularity=> 'DEFAULT',  block_sample=> FALSE,  cascade=> TRUE,  degree=> DBMS_STATS.DEFAULT_DEGREE,  method_opt=> 'FOR ALL COLUMNS SIZE AUTO');\n";
}
foreach (@indexes) {
  print "exec dbms_stats.gather_index_stats( ownname=> '$ENV{USER}',  indname=> '$_');\n";
}
