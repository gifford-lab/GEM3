#!/usr/bin/perl

use strict;
use warnings;
use DBI;

use PSRG;
use PSRG::Database;

my $old = PSRG::Database::handleForRole('psrg');
my $new = PSRG::Database::handleForRole('cgs');

my %tables = ('arolfe.orth_mapping'=>'ylannotations.orth_mapping',
'arolfe.orth_pair'=>'ylannotations.orth_pair',
'arolfe.weightmatrix'=>'ylannotations.weightmatrix',
'arolfe.weightmatrixcols'=>'ylannotations.weightmatrixcols',
'arolfe.weightmatrixscan'=>'ylannotations.weightmatrixscan',
'arolfe.wms_properties'=>'ylannotations.wms_properties',
'arolfe.wms_scanned_regions'=>'ylannotations.wms_scanned_regions',
'arolfe.wms_hits'=>'ylannotations.wms_hits',
'arolfe.syntenyscan'=>'ylannotations.syntenyscan',
'arolfe.syntenygenome'=>'ylannotations.syntenygenome',
'arolfe.syntenyblock'=>'ylannotations.syntenyblock',
'arolfe.arraydesign'=>'ylchipchip.arraydesign',
'arolfe.galfiles'=>'ylchipchip.galfiles',
'arolfe.probedesign'=>'ylchipchip.probedesign',
'arolfe.probet'=>'ylchipchip.probet',
'arolfe.probelocation'=>'ylchipchip.probelocation',
'arolfe.fragdist'=>'ylchipchip.fragdist',
'arolfe.fragdistentry'=>'ylchipchip.fragdistentry',
'arolfe.experiment'=>'ylchipchip.experiment',
'arolfe.exptToGenome'=>'ylchipchip.exptToGenome',
'arolfe.data'=>'ylchipchip.data',
'arolfe.datatemp'=>'ylchipchip.datatemp',
'arolfe.datatemp2'=>'ylchipchip.datatemp2',
'arolfe.mleanalysis'=>'ylchipchip.mleanalysis',
'arolfe.mleparameters'=>'ylchipchip.mleparameters',
'arolfe.mleanalysisinputs'=>'ylchipchip.mleanalysisinputs',
'arolfe.mleToGenome'=>'ylchipchip.mleToGenome',
'arolfe.mleresults'=>'ylchipchip.mleresults',
'arolfe.bayesanalysis'=>'ylchipchip.bayesanalysis',
'arolfe.bayesparameters'=>'ylchipchip.bayesparameters',
'arolfe.bayesanalysisinputs'=>'ylchipchip.bayesanalysisinputs',
'arolfe.bayesToGenome'=>'ylchipchip.bayesToGenome',
'arolfe.bayesresults'=>'ylchipchip.bayesresults',
'arolfe.rosettaanalysis'=>'ylchipchip.rosettaanalysis',
'arolfe.rosettaparameters'=>'ylchipchip.rosettaparameters',
'arolfe.rosettaanalysisinputs'=>'ylchipchip.rosettaanalysisinputs',
'arolfe.rosettaToGenome'=>'ylchipchip.rosettaToGenome',
'arolfe.rosettaresults'=>'ylchipchip.rosettaresults',
'arolfe.ipmeta'=>'ylchipchip.ipmeta',
'arolfe.hybmeta'=>'ylchipchip.hybmeta',
'arolfe.scanmeta'=>'ylchipchip.scanmeta',
'arolfe.bindingscan'=>'ylchipchip.bindingscan',
'arolfe.bindingscanToExpt'=>'ylchipchip.bindingscanToExpt',
'arolfe.bindingscanToGenome'=>'ylchipchip.bindingscanToGenome',
'arolfe.bindingscanregion'=>'ylchipchip.bindingscanregion',
'arolfe.bindingscanparam'=>'ylchipchip.bindingscanparam',
'arolfe.bindingevent'=>'ylchipchip.bindingevent',
'arolfe.species'=>'core.species',
'arolfe.genome'=>'core.genome',
'arolfe.chromosome'=>'core.chromosome',
'arolfe.chromsequence'=>'core.chromsequence',
'arolfe.condition'=>'core.condition',
'arolfe.cells'=>'core.cells',
'arolfe.factors'=>'core.factors');

foreach my $oldname (keys %tables) {
  my $newname = $tables{$oldname};
  my ($oldcount) = $old->selectrow_array("select count(*) from $oldname");
  my ($newcount) = $new->selectrow_array("select count(*) from $newname");
  unless ($oldcount == $newcount) {
    print "$newname : $oldcount, $newcount\n";
  }
}
