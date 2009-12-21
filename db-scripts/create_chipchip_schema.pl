#!/usr/bin/perl

# this creates new oracle schemas and users for the chipchip and annotations schemas

use strict;
use warnings;
use DBI;
use DBD::Oracle qw(:ora_session_modes);
use Getopt::Long;
use Crypt::RandPasswd;

my ($syspasswd,$prefix);
my $defaulttablespace = 'CGS-SDB1';
my $chipchipschemafile = '/oracle/scripts/chipchip.oracle';
my $annotsschemafile = '/oracle/scripts/annotations.oracle';
my $expressionschemafile = '/oracle/scripts/expression.oracle';
unless (-e $chipchipschemafile) {
  die "Can't find $chipchipschemafile";
}
unless (-e $annotsschemafile) {
  die "Can't find $annotsschemafile";
}
unless (-e $expressionschemafile) {
  die "Can't find $expressionschemafile";
}
my $grouppwfile = '/oracle/scripts/oracle_passwords';
my $superuserrole = 'cgs';
GetOptions("prefix=s"=>\$prefix);
unless ($prefix) {
  die "Must supply --prefix   (eg fink, young, odom)";
}

print STDERR "If this doesn't ask you for the sys password in a few seconds, kill it and try\n";
print STDERR "again.  The password generation sometimes hangs.\n";

my $chipchippasswd = Crypt::RandPasswd->word(9,12);
my $annotspasswd = $chipchippasswd;
my $expressionpasswd = $chipchippasswd;
my $publicpasswd = Crypt::RandPasswd->word(9,12);
my $pwfilesdir = "/oracle/scripts/${prefix}";
my $chipchipuser = $prefix . 'chipchip';
my $annotsuser = $prefix . 'annotations';
my $expressionuser = $prefix . 'expression';
my $publicuser = $prefix . 'public';

if (-d $pwfilesdir) {
  die "$pwfilesdir already exists.  I won't overwrite anything, so I'm quitting";
}
mkdir($pwfilesdir);

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
open(PWFILE,">${pwfilesdir}/chipchip_passwd") or die "Can't open ${pwfilesdir}/chipchip_passwd : $!";
print PWFILE "dbiconnectstring=dbi:Oracle:host=olig2.csail.mit.edu;sid=cgs\n";
print PWFILE "jdbcconnectstring=jdbc:oracle:thin:\@olig2.csail.mit.edu:1521:cgs\n";
print PWFILE "user=${publicuser}\nschema=${chipchipuser}\npasswd=${publicpasswd}\n";
close PWFILE;
open(PWFILE,">${pwfilesdir}/annotations_passwd") or die "Can't open ${pwfilesdir}/annotations_passwd : $!";
print PWFILE "dbiconnectstring=dbi:Oracle:host=olig2.csail.mit.edu;sid=cgs\n";
print PWFILE "jdbcconnectstring=jdbc:oracle:thin:\@olig2.csail.mit.edu:1521:cgs\n";
print PWFILE "user=${publicuser}\nschema=${annotsuser}\npasswd=${publicpasswd}\n";
close PWFILE;
open(PWFILE,">${pwfilesdir}/expression_passwd") or die "Can't open ${pwfilesdir}/expression_passwd : $!";
print PWFILE "dbiconnectstring=dbi:Oracle:host=olig2.csail.mit.edu;sid=cgs\n";
print PWFILE "jdbcconnectstring=jdbc:oracle:thin:\@olig2.csail.mit.edu:1521:cgs\n";
print PWFILE "user=${publicuser}\nschema=${expressionuser}\npasswd=${publicpasswd}\n";
close PWFILE;



print STDERR "Appending to ${grouppwfile}\n";
open(PWFILE,">>${grouppwfile}") or die "Can't open $grouppwfile : $!";
print PWFILE "${chipchipuser}\t${chipchippasswd}\n";
print PWFILE "${annotsuser}\t${annotspasswd}\n";
print PWFILE "${expressionuser}\t${expressionpasswd}\n";
print PWFILE "${publicuser}\t${publicpasswd}\n";
close PWFILE;

print STDERR "Creating users\n";
$dbh->do("create user ${chipchipuser} identified by ${chipchippasswd} " .
	 "default tablespace \"${defaulttablespace}\" quota unlimited on \"${defaulttablespace}\" ");
foreach my $role (("ALTER SESSION","CREATE PROCEDURE","CREATE SEQUENCE","CREATE SESSION","CREATE TABLE","CREATE TRIGGER","CREATE TYPE","CREATE VIEW")) {
  $dbh->do("grant ${role} to ${chipchipuser}");
}

$dbh->do("create user ${annotsuser} identified by ${annotspasswd} " .
	 "default tablespace \"${defaulttablespace}\" quota unlimited on \"${defaulttablespace}\"");
foreach my $role (("CONNECT","ALTER SESSION","CREATE PROCEDURE","CREATE SEQUENCE","CREATE SESSION","CREATE TABLE","CREATE TRIGGER","CREATE TYPE","CREATE VIEW")) {
  $dbh->do("grant ${role} to ${annotsuser}");
}

$dbh->do("create user ${expressionuser} identified by ${expressionpasswd} " .
	 "default tablespace \"${defaulttablespace}\" quota unlimited on \"${defaulttablespace}\"");
foreach my $role (("CONNECT","ALTER SESSION","CREATE PROCEDURE","CREATE SEQUENCE","CREATE SESSION","CREATE TABLE","CREATE TRIGGER","CREATE TYPE","CREATE VIEW")) {
  $dbh->do("grant ${role} to ${expressionuser}");
}


$dbh->do("create user ${publicuser} identified by ${publicpasswd}");
$dbh->do("grant CONNECT to ${publicuser}");

print STDERR "Populating schemas\n";
system("sqlplus ${chipchipuser}/${chipchippasswd}\@cgs < ${chipchipschemafile}");
system("sqlplus ${annotsuser}/${annotspasswd}\@cgs < ${annotsschemafile}");
system("sqlplus ${expressionuser}/${expressionpasswd}\@cgs < ${expressionschemafile}");

print STDERR "Granting permissions\n";
my $stmt = $dbh->prepare("select owner, table_name from dba_tables where lower(owner) in('${chipchipuser}','${annotsuser}','${expressionuser}','core')");
$stmt->execute();
while (my @r = $stmt->fetchrow_array()) {
  $dbh->do("grant all on $r[0].$r[1] to cgs");
  $dbh->do("grant select on $r[0].$r[1] to ${publicuser}");
}

$stmt = $dbh->prepare("select sequence_owner, sequence_name from dba_sequences where lower(sequence_owner) in('${chipchipuser}','${annotsuser}','${expressionuser}')");
$stmt->execute();
while (my @r = $stmt->fetchrow_array()) {
  $dbh->do("grant all on $r[0].$r[1] to cgs");
}
