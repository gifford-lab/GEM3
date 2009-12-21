--  -*- mode:sql  -*-

create table species (
	id integer,
	name varchar 
	);

create table genome (
	id integer,
	species integer,
	version varchar,
	description varchar
	);

create table chromosome (
	id integer,
	name varchar,
	genome integer
	);

create table chromsequence (
	id integer,
	sequence varchar
	);

create table condition (
	id integer,
	name varchar
	);

create table cells (
	id integer,
	name varchar
	);

create table factors (
	id integer,
	name varchar
	);

create table timeseries (
	id integer,
	name varchar
	);

create table timepoint (
	id integer,
	time_series integer,
	series_order integer,
	name varchar
	);


