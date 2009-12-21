#!/usr/bin/perl

use strict;
use warnings;
use PSRG::Database;
# can operate in two modes:
# If no --role options are provided: 
#     takes a set of SQL filesnames on command line.
#     tries to guess what kind of chromosome, startpos index is needed for
#     each one and emits the corresponding SQL on STDOUT
# if --role is provided, connects to that database.  If --table
#     options are provided, works on those tables.  Otherwise works
#     on all tables.  If it can't find an index on the table, it generates
#     one

my @tables = ();
my @roles = ();
use Getopt::Long;
GetOptions("--role=s"=>\@roles,
	   "--table=s"=>\@tables);
if (@roles) {
  foreach (@roles) {
    do_from_role($_,@tables);
  }
} else {
  do_from_files();
}


sub do_from_role {
  my ($role,@tables) = @_;
  my $dbh = PSRG::Database::handleForRole($role) or
    die "No handle for $role";
  unless (@tables) {
    my $sth = $dbh->prepare('show tables');
    $sth->execute();
    while (my @r = $sth->fetchrow_array()) {
      push(@tables,$r[0]);
    }
    print STDERR "Defaulting to tables @tables\n";
  }
  foreach (@tables) {
    my $sth = $dbh->prepare("show create table $_");
    my $schema = '';
    $sth->execute();
    while (my @r = $sth->fetchrow_array()) {
      $schema .= join("",@r);
    }
    my @l = sql_from_schema($_,$schema);
    foreach my $sql (@l) {
      print STDERR "$sql\n";
      $dbh->do($sql);
    }
  }
  $dbh->disconnect();
}




sub do_from_files {
 TABLE: foreach my $fname (@ARGV) {
    open(FILE,$fname) or die "Can't open $fname : $!";
    my $tname= $fname;
    $tname =~ s/.sql$//;
    my $lines = join('',<FILE>);
    print sql_from_schema($lines);
  }
}

sub sql_from_schema {
  my ($tname,$lines) = @_;
  my $prefix;
  if ($lines =~/\WtxStart\W/) {
    $prefix = 'tx';
  } elsif ($lines =~ /\WtStart\W/) {
    $prefix = 't';
  } elsif ($lines =~ /\WchromStart\W/) {
    $prefix = 'chrom';
  } else {
    warn "Can't find a suitable position column in $tname";
    print STDERR "LINES were $lines\n";
    return ();
  }
  my @cols = ();
  my ($onechrom) = ($tname  =~ /^chr/);
  # see if there's already a suitable index
  if ($onechrom) {
    if ($lines =~ /KEY.*Start/i) {
      return ();
    }
  } else {
    if ($lines =~ /KEY.*chrom.*Start/i) {
      return ();
    }
  }
  if (not $onechrom) {
    if ($lines =~ /\Wchrom\W/) {
      push(@cols,'chrom');
    } elsif ($lines =~ /\WtName\W/) {
      push(@cols,'tName');
    }
  }
  
  push(@cols,"${prefix}Start");
  push(@cols,"${prefix}End");
  @cols = grep {$lines =~ /$_/} @cols;
  my $cols = join(',',@cols);
  return ("alter table ${tname} drop index ix_${tname}_chrompos;\n",
	  "create index ix_${tname}_chrompos on $tname($cols);\n");
}
