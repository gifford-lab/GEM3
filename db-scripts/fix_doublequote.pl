#!/usr/bin/perl

use strict;
use warnings;
use DBI;

my ($host,$sid,$user) = ('opteron.csail.mit.edu','psrg',$ENV{USER});
my $passwd = `cat ~/.oracle_passwd`;
chomp($passwd);
my $dbh = DBI->connect("dbi:Oracle:host=$host;sid=$sid", $user, $passwd);
$dbh->do("alter session set current_schema=arolfe");
my @names = qw(YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_2248_2568_238
YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_1606_1926_229
	       YIL154C_IMP2'_SGDID:S0001416,_Chr_IX_from_55021-53981,_reverse_complement,_Verified_ORF_781_1040_186
YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_1285_1605_3
	       YIL154C_IMP2'_SGDID:S0001416,_Chr_IX_from_55021-53981,_reverse_complement,_Verified_ORF_1_260_142
YIL154C_IMP2'_SGDID:S0001416,_Chr_IX_from_55021-53981,_reverse_complement,_Verified_ORF_521_780_103
	       YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_964_1284_58
YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_1927_2247_77
	       YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_322_642_240
YIL154C_IMP2'_SGDID:S0001416,_Chr_IX_from_55021-53981,_reverse_complement,_Verified_ORF_261_520_110
	       YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_643_963_244
YHR047C_AAP1'_SGDID:S0001089,_Chr_VIII_from_201302-198732,_reverse_complement,_Verified_ORF_1_321_32);

print "There are $#names names\n";
my $selecth = $dbh->prepare("select probename, probeid from probedesign where probeid = ?");
my $updateh = $dbh->prepare("update probedesign set probeid = ?, probename = ? where probeid = ?");
foreach my $name (@names) {
  my $dqpi = $name;
  $dqpi =~ s/\'/\'\'/g;
  $selecth->execute($dqpi);
  my @results = $selecth->fetchrow_array();
  my $dqpn = $results[0];  
  my $sqpn = $dqpn;
  $sqpn =~ s/\'\'/\'/g;
  unless ($results[1] eq $dqpi) {
    warn "Couldn't row for $name : @results";
  }
  $updateh->execute($name,$sqpn,$dqpi);
}
