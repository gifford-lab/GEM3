-- MySQL dump 9.10
--
-- Host: localhost    Database: sacCer1
-- ------------------------------------------------------
-- Server version	4.0.17-standard

--
-- Table structure for table `sgdOther`
--

CREATE TABLE sgdOther (
  bin smallint(6) NOT NULL default '0',
  chrom varchar(255) NOT NULL default '',
  chromStart int(11) NOT NULL default '0',
  chromEnd int(11) NOT NULL default '0',
  name varchar(255) NOT NULL default '',
  score int(11) NOT NULL default '0',
  strand char(1) NOT NULL default '',
  type varchar(255) NOT NULL default '',
  KEY chrom (chrom(8),bin)
) TYPE=MyISAM;
