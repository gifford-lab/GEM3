#!/usr/bin/perl

use strict;
use warnings;

my $user=$ARGV[0];
print STDERR "Please enter your password: ";

system("stty -echo");
my $pass = <stdin>;
system("stty echo");

print "expdp ${user}/${pass} CONTENT=DATA_ONLY SCHEMAS=${user} DIRECTORY=AROLFE_BACKUP DUMPFILE=${user}.dmp\n";
