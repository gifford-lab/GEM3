--  -*- mode:sql  -*-

create table arraydesign (
	id integer,
	name varchar,
	genome integer
	);

create table galfiles (
	name varchar,
	id integer
	);

create table probedesign (
	id integer,
	arraydesign integer,
	blockno integer,
	colno integer,
	rowno integer,
	galfile integer,
	probename varchar,
	probeid varchar,
	type varchar,
	sequence varchar
	);

create table probetm(
	id integer,
	tm integer
	);

create table probelocation (
	id integer,
	chromosome integer,
	startpos integer,
	stoppos integer,
	strand varchar,
	loccount integer,
	bitscore integer
	);

create table fragdist (
	id integer,
	name varchar,
	version varchar,
	description varchar
	);

create table fragdistentry (
	distribution integer,
	distance integer,
	value integer
	);	

create table experiment (
	id integer,
	name varchar,
	version varchar,
	replicate varchar,
	fragdist integer,
	species integer,
	cellsone integer,
	conditionone integer,
	factorone varchar,
	cellstwo integer,
	conditiontwo integer,
	factortwo varchar,
	normalization varchar,
	active integer
	);
	
create table exptMetadata (
	experiment integer,
	key varchar,
	value clob not null
	);

create table exptToGenome (
	experiment integer,
	genome integer
	);

create table data (
	experiment integer,
	probe integer,
	channelone float,
	channeltwo float,
	mor float,
	channelratio float,
	ratio float
	);

create table bayesanalysis (
	id integer,
	species integer,
	name varchar,
	version varchar,
	active smallint
	);

create table bayesparameters (
	analysis integer,
	name varchar,
	value varchar
	);

create table bayesanalysisinputs (
	analysis integer,
	experiment integer
	);

create table bayesToGenome (
	analysis integer,
	genome integer
	);

create table bayesresults (
	analysis integer,
	chromosome integer,
	position integer,
	posterior float,
	posteriorstd float,
	strength float,
	strengthstd float,
	);

create table bindingscan (
	id integer,
	version varchar,
	type varchar
	);

create table bindingscanToExpt (
	scan integer,
	scanexpt integer,
	scantype integer
	);

create table bindingscanToGenome (
	scan integer,
	genome integer
	);

create table bindingscanregion (
	scan integer,
	chromosome integer,
	startpos integer,
	stoppos integer
	);

create table bindingscanparam (
	scan integer,
	key varchar,
	value varchar
	);

create table bindingevent (
	scan integer,
	chromosome integer,
	startpos integer,
	stoppos integer,
	eventsize float,
	eventconf float
	);

