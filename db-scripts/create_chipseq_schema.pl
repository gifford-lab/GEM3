#!/usr/bin/perl

# this creates new oracle schemas and users for the chipseq system

use strict;
use warnings;
use DBI;
use DBD::Oracle qw(:ora_session_modes);
use Getopt::Long;
use Crypt::RandPasswd;

my ($syspasswd,$prefix);
my $defaulttablespace = 'CGS-SDB1';
my $chipseqschemafile = '/oracle/scripts/chipseq.oracle';
unless (-e $chipseqschemafile) {
  die "Can't find $chipseqschemafile";
}
my $grouppwfile = '/oracle/scripts/oracle_passwords';
my $superuserrole = 'cgs';
GetOptions("prefix=s"=>\$prefix);
unless ($prefix) {
  die "Must supply --prefix   (eg fink, young, odom)";
}

print STDERR "If this doesn't ask you for the sys password in a few seconds, kill it and try\n";
print STDERR "again.  The password generation sometimes hangs.\n";

my $chipseqpasswd = Crypt::RandPasswd->word(9,12);
my $publicpasswd = Crypt::RandPasswd->word(9,12);
my $pwfilesdir = "/oracle/scripts/${prefix}";
my $chipsequser = $prefix . 'chipseq';
my $publicuser = $prefix . 'public';

unless (-d $pwfilesdir) {
  die "$pwfilesdir doesn't exist. Please create it first.";
}

print "Enter sys passwd:\n";
system("stty -echo");
$syspasswd = <STDIN>;
system("stty echo");
chomp($syspasswd);
print STDERR "Connecting to database\n";
delete $ENV{TWO_TASK};
unless ($ENV{ORACLE_SID}) {die "Must set ORACLE_SID in environment";}
my $dbh = DBI->connect("dbi:Oracle:",
		       "sys",
		       $syspasswd,
		       {ora_session_mode => 2});
die "No connection" unless ($dbh);
$dbh->{RaiseError} = 1;

print STDERR "Creating password files...\n";
open(PWFILE,">${pwfilesdir}/chipseq_passwd") or die "Can't open ${pwfilesdir}/chipseq_passwd : $!";
print PWFILE "dbiconnectstring=dbi:Oracle:host=olig2.csail.mit.edu;sid=cgs\n";
print PWFILE "jdbcconnectstring=jdbc:oracle:thin:\@olig2.csail.mit.edu:1521:cgs\n";
print PWFILE "user=${publicuser}\nschema=${chipsequser}\npasswd=${publicpasswd}\n";
close PWFILE;

print STDERR "Appending to ${grouppwfile}\n";
open(PWFILE,">>${grouppwfile}") or die "Can't open $grouppwfile : $!";
print PWFILE "${chipsequser}\t${chipseqpasswd}\n";
close PWFILE;

print STDERR "Creating users\n";
$dbh->do("create user ${chipsequser} identified by ${chipseqpasswd} " .
	 "default tablespace \"${defaulttablespace}\" quota unlimited on \"${defaulttablespace}\" ");
foreach my $role (("ALTER SESSION","CREATE PROCEDURE","CREATE SEQUENCE","CREATE SESSION","CREATE TABLE","CREATE TRIGGER","CREATE TYPE","CREATE VIEW")) {
  $dbh->do("grant ${role} to ${chipsequser}");
}

print STDERR "Populating schemas\n";
system("sqlplus ${chipsequser}/${chipseqpasswd}\@cgs < ${chipseqschemafile}");

print STDERR "Granting permissions\n";
my $stmt = $dbh->prepare("select owner, table_name from dba_tables where lower(owner) in('${chipsequser}','core')");
$stmt->execute();
while (my @r = $stmt->fetchrow_array()) {
  $dbh->do("grant all on $r[0].$r[1] to cgs");
  $dbh->do("grant select on $r[0].$r[1] to ${publicuser}");
}

$stmt = $dbh->prepare("select sequence_owner, sequence_name from dba_sequences where lower(sequence_owner) in('${chipsequser}')");
$stmt->execute();
while (my @r = $stmt->fetchrow_array()) {
  $dbh->do("grant all on $r[0].$r[1] to cgs");
}
