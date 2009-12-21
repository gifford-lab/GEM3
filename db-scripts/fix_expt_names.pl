#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use PSRG::Database;

my $dbh = PSRG::Database::handleForRole('chipchip');


fixTable('experiment','name');
fixTable('bayesanalysis','name');
fixTable('bindingscan','version');

sub fixTable {
  my ($table,$col) = @_;
  my $sth = $dbh->prepare("select id, ${col} from $table");
  my $fix = $dbh->prepare("update $table set ${col} = ? where id = ?");
  $sth->execute();
  while (my ($id,$name) = $sth->fetchrow_array()) {
    my $newname = undef;
    if ($name =~ /Mm.*Olig2.differentiating/) {
      $newname = $name;
      $newname =~ s/Olig2.differentiating/HBG3:Olig2 Stage/g;
    }
    if ($newname) {
      my $space = ' ' x (80 - length($name));
      print "${name}${space}${newname}\n";
      $fix->execute($newname,$id);
    }
  }
}
