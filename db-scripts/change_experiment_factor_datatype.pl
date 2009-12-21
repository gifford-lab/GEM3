#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use DBD::Oracle qw(:ora_types);
use Getopt::Long;
use PSRG::Utils;
use PSRG::Database;
my $dbhcore = PSRG::Database::handleForRole('core');
my $dbhchip = PSRG::Database::handleForRole('chipchip');

# this program changes the datatype of the factorone and factortwo columns
# in the experiment table.  They started out as strings and will be numbers, just
# like cells and condition.
#
# the factors table must be created in the core schema before this program is run
#

$dbhchip->{AutoCommit} = 0;
$dbhchip->{RaiseError} = 1;
$dbhcore->{AutoCommit} = 0;
$dbhcore->{RaiseError} = 1;

$dbhchip->do("alter table experiment rename column factortwo to factortwostring");
$dbhchip->do("alter table experiment rename column factorone to factoronestring");
$dbhchip->do("alter table experiment add (factorone number(10), factortwo number(10))");

my %factors = ();
my $addfactor = $dbhcore->prepare("insert into factors(id,name) values (factors_id.nextval,?)");
my $getfactorid = $dbhcore->prepare("select factors_id.currval from dual");
my $getfactors = $dbhcore->prepare("select id, name from factors");
$getfactors->execute();
while (my @r = $getfactors->fetchrow_array()) {
  $factors{$r[1]} = $r[0];
}
my $updatefactors = $dbhchip->prepare("update experiment set factorone = ?, factortwo = ? where id = ?");

my $getexpts = $dbhchip->prepare("select id, factoronestring, factortwostring from experiment");
$getexpts->execute();
while (my ($id, $f1, $f2) = $getexpts->fetchrow_array()) {
  my $f1val = getFactor($f1);
  my $f2val = getFactor($f2);
  print STDERR "$id : $f1 -> $f1val, $f2 -> $f2val\n";
  $updatefactors->execute($f1val,$f2val,$id);
}

$dbhchip->commit();
$dbhcore->commit();
$dbhchip->disconnect();
$dbhcore->disconnect();

sub getFactor {
  my ($factor) = @_;
  if ($factors{$factor}) {
    return $factors{$factor};
  } else {
    $addfactor->execute($factor);
    $getfactorid->execute();
    my @r = $getfactorid->fetchrow_array();
    my $id = $r[0] || die "Can't get id for $factor";
    $factors{$factor} = $id;
    return $id;
  }
}



