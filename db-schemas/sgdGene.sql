-- MySQL dump 9.10
--
-- Host: localhost    Database: sacCer1
-- ------------------------------------------------------
-- Server version	4.0.17-standard

--
-- Table structure for table `sgdGene`
--

CREATE TABLE sgdGene (
  name varchar(255) NOT NULL default '',
  chrom varchar(255) NOT NULL default '',
  strand char(1) NOT NULL default '',
  txStart int(10) unsigned NOT NULL default '0',
  txEnd int(10) unsigned NOT NULL default '0',
  cdsStart int(10) unsigned NOT NULL default '0',
  cdsEnd int(10) unsigned NOT NULL default '0',
  exonCount int(10) unsigned NOT NULL default '0',
  exonStarts longblob NOT NULL,
  exonEnds longblob NOT NULL,
  proteinID varchar(40) NOT NULL default '',
  KEY name (name(16)),
  KEY chrom (chrom(8),txStart),
  KEY chrom_2 (chrom(8),txEnd),
  KEY proteinID (proteinID(10))
) TYPE=MyISAM;
