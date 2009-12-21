create table chiplist(
	filename varchar2(300),
	designname varchar2(100),
	species varchar2(100),
	cellsone varchar2(100),
	cellstwo varchar2(100),
	conditionone varchar2(100),
	conditiontwo varchar2(100),
	factorone varchar2(100),
	factortwo varchar2(100),
	active number(1),
	loaded number(1));

create sequence species_id;
create table species (
	id number(10) constraint species_id unique not null,
	name varchar2(80) constraint species_pk primary key);

create sequence genome_id;
create table genome (
	id number(10) constraint genome_id unique not null,
	species constraint fk_genome_species references species(id) not null,
	version varchar2(100) constraint genome_version not null,
	constraint genome_pk primary key(species,version));

create sequence chromosome_id;
create table chromosome (
	id number(10) constraint chromosome_id unique not null,
	name varchar2(100) constraint chromosome_name not null,
	genome constraint fk_chromosome_genome references genome(id) not null,
	constraint chromosome_pk primary key(name,genome));
create index ix_chromosome_name on chromosome(genome,name);

create table chromsequence (
	id number(10) constraint fk_chromseq_id references chromosome(id) not null,
	sequence clob);
create index ix_chromsequence_id on chromsequence(id);

create sequence condition_id;
create table condition (
	id number(10) constraint condition_u unique not null,
	name varchar2(100) constraint condition_pk primary key);

create sequence cells_id;
create table cells (
	id number(10) constraint cells_id unique not null,
	name varchar2(100) constraint cells_pk primary key);

create sequence arraydesign_id;
create table arraydesign (
	id number(10) constraint arraydesign_id unique not null,
	name varchar2(100) constraint arraydesign_u unique not null,
	genome constraint fk_arraydesign_genome references genome(id) not null,
	constraint arraydesign_id_pk primary key(name,genome));

create sequence galfile_id;
create table galfiles (
	name varchar2(300) constraint galfiles_name primary key,
	id number(10) constraint galfiles_id unique not null);
	

create sequence probedesign_id;
create table oldprobedesign (
	id number(10) constraint probedesign_idunique unique not null,
	arraydesign constraint fk_probedesign_arraydesign references arraydesign(id),
	blockno number(10) constraint probedesign_block not null,
	colno number(10) constraint probedesign_col not null,
	rowno number(10) constraint probedesign_row not null,
	galfile number(10) constraint fk_probedesign_galfile references galfiles(id),
	probename varchar2(200),
	probeid varchar2(200),
	type varchar2(100),
	sequence varchar2(2000),
	chromosome constraint fk_probedesign_chromosome references chromosome(id) not null,
	startpos number(10) not null,
	stoppos number(10) not null,
	constraint probedesign_unique_pos unique (galfile,probeid,blockno,colno,rowno),
	constraint probedesign_pk primary key (chromosome, startpos, stoppos, id))
organization index 
compress 1
pctthreshold 2
overflow including arraydesign;

create table oldprobelocation (
	id number(10) constraint fk_probeloc_probe references probedesign(id) not null,
	chromosome constraint fk_probeloc_chromosome references chromosome(id) not null,
	startpos number(10) not null,
	stoppos number(10) not null,
	constraint probeloc_pk primary key (chromosome, startpos, stoppos, id))
organization index 
compress 1
pctthreshold 2;
create index ix_probelocation_id on probelocation(id,chromosome);

create table probedesign (
	id number(10) constraint newprobedesign_idunique not null,
	arraydesign constraint fk_newprobedesign_arraydesign references arraydesign(id),
	blockno number(10) constraint newprobedesign_block not null,
	colno number(10) constraint newprobedesign_col not null,
	rowno number(10) constraint newprobedesign_row not null,
	galfile number(10) constraint fk_newprobedesign_galfile references galfiles(id),
	probename varchar2(200),
	probeid varchar2(200),
	type varchar2(100),
	sequence varchar2(2000),
	constraint newprobedesign_unique_pos unique (galfile,probeid,blockno,colno,rowno),
	constraint newprobedesign_pk primary key (id))
organization index 
pctthreshold 2
overflow including arraydesign;

create index ix_newprobedesign_design on probedesign(arraydesign);
create index ix_newprobedesign_probename on probedesign(probename);
create index ix_newprobedesign_probeid on probedesign(probeid);
create index ix_newprobedesign_galfile on probedesign(galfile);
create index ix_newprobedesign_type on probedesign(type);

create table probetm(
	id number(10) constraint fk_probetm_probe references probedesign(id) not null,	
	tm number(10,4));
create index ix_probetm_id on probetm(id);

create table probelocation (
	id number(10) constraint fk_newprobeloc_probe references probedesign(id) not null,
	chromosome constraint fk_newprobeloc_chromosome references chromosome(id) not null,
	startpos number(10) not null,
	stoppos number(10) not null,
	loccount number(10),
	bitscore number(10),
	constraint newprobeloc_pk primary key (chromosome, startpos, stoppos, id))
organization index 
compress 1;
create index ix_newprobelocation_id on newprobelocation(id);


create sequence fragdist_id;
create table fragdist (
	id number(10) constraint fragdist_id unique not null,
	name varchar2(100),
	version varchar2(100),
	description varchar2(1000),
	constraint fragdist_id_pk primary key(name, version));
create table fragdistentry (
	distribution constraint fk_fragdistentry_dist references fragdist(id) not null,
	distance number(6),
	value number);	

create sequence experiment_id;
create table experiment (
	id number(10) constraint experiment_id unique not null,
	name varchar2(100) constraint experiment_name not null,
	version varchar2(100),
	replicate varchar2(40),
	fragdist constraint fk_experiment_fragdist references fragdist(id) not null,
	species constraint fk_experiment_species references species(id) not null,
	cellsone constraint fk_experiment_cellsone references cells(id) not null,
	conditionone constraint fk_experiment_conditionone references condition(id) not null,
	factorone varchar2(40) constraint experiment_factorone not null,
	cellstwo constraint fk_experiment_cellstwo references cells(id) not null,
	conditiontwo constraint fk_experiment_conditiontwo references condition(id) not null,
	factortwo varchar2(40) constraint experiment_factortwo not null,
	normalization varchar2(1000),
	active smallint,
	constraint experiment_pk primary key(name, species, version, replicate));

create table exptToGenome (
	experiment constraint fk_expt2genomeexpt references experiment(id) not null,
	genome constraint fk_expt2genomegenome references genome(id) not null);
create index ix_expt2genomee on exptToGenome(experiment);
create index ix_expt2genomeg on exptToGenome(genome);

create table olddata (
	experiment constraint fk_data_experiment references experiment(id) not null,
	probe constraint fk_data_probe references probedesign(id) not null,
	channelone number(10,2),
	channeltwo number(10,2),
	ratio number(10,4));
create index ix_data_probeexpt on data(experiment,probe);
create index ix_data_expt on data(experiment);
create index ix_data_probe on data(probe);

create table data (
	experiment constraint fk_newdata_experiment references experiment(id) not null,
	probe constraint fk_newdata_probe references probedesign(id) not null,
	channelone binary_float,
	channeltwo binary_float,
	mor binary_float,
-- ROM from channelone / channeltwo
	channelratio binary_float,
-- final output ratio
	ratio binary_float);
create index ix_data_exptprobe on data(experiment,probe);
create index ix_data_probe on data(probe);

-- used to some of the normalization procedures
create table datatemp (
	experiment constraint fk_datatemp_experiment references experiment(id) not null,
	probe constraint fk_datatemp_probe references probedesign(id) not null,
	channelone binary_float,
	channeltwo binary_float,
	mor binary_float,
-- ROM from channelone / channeltwo
	channelratio binary_float,
-- final output ratio
	ratio binary_float,
	controlratio binary_float);

create table datatemp2 (
	experiment constraint fk_datatemp2_experiment references experiment(id) not null,
	probe constraint fk_datatemp2_probe references probedesign(id) not null,
	channelone binary_float,
	channeltwo binary_float,
	mor binary_float,
-- ROM from channelone / channeltwo
	channelratio binary_float,
-- final output ratio
	ratio binary_float,
	controlratio binary_float);

create sequence analysis_id;
create table mleanalysis (
	id number(10) constraint analysis_id unique not null,
	species constraint fk_analysis_species references species(id),
	name varchar2(100),
	version varchar2(200),
	active smallint,
	 constraint analysis_pk primary key (name, species, version));
create table mleparameters (
	analysis constraint fk_parameters_analysis references mleanalysis(id) not null,
	name varchar2(100),
	value varchar2(1000));
create table mleanalysisinputs (
	analysis constraint fk_inputs_analysis references mleanalysis(id) not null,
	experiment constraint fk_inputs_expt references experiment(id) not null);

create table mleToGenome (
	analysis constraint fk_mle2genomeexpt references mleanalysis(id) not null,
	genome constraint fk_ml2genomegenome references genome(id) not null);
create index ix_mle2genomea on mleToGenome(analysis);
create index ix_mle2genomeg on mleToGenome(genome);


-- create table oldmleresults (
-- 	analysis constraint fk_mleresults_analysis references mleanalysis(id) not null,
-- 	chromosome constraint fk_mleresults_chromosome references chromosome(id) not null,
-- 	position number(10),
-- 	b_i number(10,4),
-- 	bindll number(20,2),
-- 	nullll number(20,2),
-- 	lograt number(20,2),
-- 	conf number(10,4));
-- create index ix_mleresults_acp on mleresults(analysis,chromosome,position);
-- create index ix_mleresults_location on mleresults(chromosome,position);

create table mleresults (
	analysis constraint fk_mleresults_analysis references mleanalysis(id) not null,
	chromosome constraint fk_mleresults_chromosome references chromosome(id) not null,
	position number(10),
	b_i number(10,4),
	bindll number(20,2),
	nullll number(20,2),
	lograt number(20,2),
	conf number(10,4),
	constraint pk_mleresults primary key (analysis, chromosome, position))
organization index compress 2;


create table bayesanalysis (
	id number(10) constraint bayesanalysis_id unique not null,
	species constraint fk_bayesanalysis_species references species(id),
	name varchar2(100),
	version varchar2(200),
	active smallint,
	constraint bayesianalysis_pk primary key (name, species, version));
create table bayesparameters (
	analysis constraint fk_bayesparameters_analysis references bayesanalysis(id) not null,
	name varchar2(100),
	value varchar2(1000));
create table bayesanalysisinputs (
	analysis constraint fk_bayesinputs_analysis references bayesanalysis(id) not null,
	experiment constraint fk_bayesinputs_expt references experiment(id) not null);

create table bayesToGenome (
	analysis constraint fk_bayes2genomeexpt references bayesanalysis(id) not null,
	genome constraint fk_bayes2genomegenome references genome(id) not null);
create index ix_bayes2genomea on bayesToGenome(analysis);
create index ix_bayes2genomeg on bayesToGenome(genome);


-- create table oldbayesresults (
-- 	analysis constraint fk_bayesresults_analysis references bayesanalysis(id) not null,
-- 	chromosome constraint fk_bayesresults_chromosome references chromosome(id) not null,	
-- 	step number(4),
-- 	startpos number(10) not null,
-- 	stoppos number(10) not null,
-- 	maxposterior number(5,3) not null,
-- 	maxstrength number(5,3) not null,
-- 	posterior blob,
-- 	posteriorstd blob,
-- 	strength blob,
-- 	strengthstd blob);


create table bayesresults (
	analysis constraint fk_bayesresults_analysis references bayesanalysis(id) not null,
	chromosome constraint fk_bayesresults_chromosome references chromosome(id) not null,	
	position number(10) not null,
	posterior binary_float,
	posteriorstd binary_float,
	strength binary_float,
	strengthstd binary_float,
	constraint pk_newbayesresults primary key(analysis,chromosome,position))
organization index compress 2;

create table rosettaparameters (
	analysis constraint fk_rosettaparameters_analysis references rosettaanalysis(id) not null,
	name varchar2(100),
	value varchar2(1000));
create table rosettaanalysisinputs (
	analysis constraint fk_rosettainputs_analysis references rosettaanalysis(id) not null,
	experiment constraint fk_rosettainputs_expt references experiment(id) not null);
create sequence rosettaanalysis_id;

create table rosettaanalysis (
	id number(10) constraint rosetta_id unique not null,
	species constraint fk_rosetta_species references species(id),
	name varchar2(100),
	version varchar2(100),
	active smallint,
	constraint rosetta_pk primary key (name, species, version));

create table rosettaToGenome (
	analysis constraint fk_rosetta2genomeexpt references rosettaanalysis(id) not null,
	genome constraint fk_rosetta2genomegenome references genome(id) not null);
create index ix_rosetta2genomea on rosettaToGenome(analysis);
create index ix_rosetta2genomeg on rosettaToGenome(genome);

create table rosettaresults (
	analysis constraint fk_rosetta_analysis references rosettaanalysis(id) not null,
	chromosome constraint fk_rosetta_chromosome references chromosome(id) not null,
	position number(10),
	ratio binary_float,
	X binary_float,
	pval binary_float,
	pval3 binary_float,
	red binary_float,
	green binary_float,
	medianofratios binary_float,
	constraint pk_rosettaresults primary key (analysis,chromosome,position))
organization index compress 2;
	

create sequence gff_id;
create table gff (
	startdate date,
	stopdate date,
	id number(10) constraint fk_gff_id unique not null,
	chromosome constraint fk_gff_chromosome references chromosome(id),
	source varchar2(100),
	type varchar2(100),
	startpos number(10),
	stoppos number(10),
	score number(20,6),
	strand char(1),
	frame char(1),
	attributes varchar2(2000))
partition by hash(chromosome) partitions 16;

create index ix_gff_date on gff(startdate, stopdate);
create index ix_gff_typesource on gff(type, source);
create index ix_gff_chrom on gff(chromosome, startpos, stoppos);
create index ix_gff_tcstartstop on gff(type, chromosome, startpos, stoppos);
create index ix_gff_tcstop on gff(type, chromosome, stoppos);
create index ix_gff_id on gff(id);
create index ix_gff_chromtypesrcpos on gff(chromosome, type, source, startpos);

create table gff_gene_name (
	id constraint fk_gff_gene_name references gff(id) unique not null,
	name varchar2(100)); 
create index ix_gff_genename on gff_gene_name(name);
create index ix_gff_genenameid on gff_gene_name(id);

create table gff_motif_name (
	id constraint fk_gff_motif_name references gff(id) unique not null,
	name varchar2(100)); 
create index ix_motif_name on gff_motif_name(name);

create table expr_experiment (
	id number(10) constraint expr_expt_id unique,
	name varchar2(1000),
	species number(10) constraint fk_expr_experiment_species references species(id)
);
create index ix_expr_experiment_id on expr_experiment(id);

create table expr (
	expt_index number(10) constraint fk_expr_index references expr_experiment(id) not null,
	probe_name varchar2(100),
	value number(20,6)
);
create index ix_expr_expt on expr(expt);
create index ix_expr_probe_name on expr(probe_name);



create table expr_class (
	id number(10) constraint fk_expr_class_id unique not null,
	name varchar2(1000)
);
create index ix_expr_class_id on expr_class(id);

create table expr_classification (
	expt number(10) constraint fk_expr_classification_expt references expr_experiment(id) not null,
	class number(10) constraint fk_expr_classification_class references expr_class(id) not null
);
create index ix_expr_classification_expt on expr_classification(expt);
create index ix_expr_classification_class on expr_classification(class);

create table expr_series (
	id number(10) constraint fk_expr_series_id unique not null,
	name varchar2(1000)
);
create index ix_expr_series_id on expr_series(id);

create table expr_series_index (
	expt number(10) constraint fk_expr_series_index_expt references expr_experiment(id) not null,
	series number(10) constraint fk_exper_series_index_series references expr_series(id) not null,
	series_order number(10) constraint fk_expr_series_order not null
);
create index ix_expr_series_index_expt on expr_series_index(expt);
create index ix_expr_series_index_series on expr_series_index(series);

create table name_relationship ( 
	parent varchar2(200) constraint fk_name_relationship_parent not null,
	child varchar2(200) constraint fk_name_relationship_child not null,
	rel_type number(5),
	version varchar2(100)
);
create index ix_name_relationship_parent on name_relationship(parent);
create index ix_name_relationship_child on name_relationship(child);

create sequence weightmatrix_id;
create table weightmatrix (
	id number(10) constraint weightmatrix_id unique not null,
	species constraint fk_weightmatrix_species references species(id) not null,
	name varchar2(200),
	version varchar2(200),
	type varchar2(100));
create index ix_wm_name on weightmatrix(species,name,version);
create table weightmatrixcols (
	weightmatrix constraint fk_wmc_id references weightmatrix(id) not null,
	position number(10) not null,
	letter char(1),
	weight binary_double);
create index ix_wmc_id on weightmatrixcols(weightmatrix,position);
create sequence weightmatrixscan_id;
create table weightmatrixscan (
	id number(10) constraint weightmatrixscan_id unique not null,	
	weightmatrix constraint fk_wms_wm references weightmatrix(id) not null,
	name varchar2(200),
	cutoff binary_float);
create table wms_properties (
	scan constraint fk_wmsprops_scan references weightmatrixscan(id) not null,
	name varchar2(100),
	value varchar2(1000));
create index ix_wms_properties_id on wmsproperties(scan);
create table wms_targets (
	scan constraint fk_wmstargets_scan references weightmatrixscan(id) not null,
	chromosome constraint fk_wmstargets_chromosome references chromosome(id) not null,
	startpos number(10),
	stoppos number(10));
create index ix_wms_targets_scanchrom on wms_targets(scan,chromosome);
create table wms_hits (
	scan constraint fk_wms_hits_scan references weightmatrixscan(id) not null,
	chromosome constraint fk_wms_hits_chromosome references chromosome(id) not null,
	startpos number(10),
	stoppos number(10),
	strand char(1),
	score binary_float,
	constraint pk_wms_hits primary key (scan,chromosome,startpos,stoppos,strand)) organization index compress 2;	


create table alignment_version (
	id number(10) constraint fk_alignment_version_id unique not null,
	name varchar2(1000)
);

create table alignment (
	id number(10),
	params varchar2(1000),
	version number(10) constraint fk_alignment_version references alignment_version(id) not null,
	score number,
	constraint aligment_pk primary key (id)
) organization index;
create index ix_alignment_version on alignment(version);

create table align_block (
	id number(10) constraint align_block_idunique unique not null,
	alignment constraint fk_align_block_alignment references alignment(id) not null,
	chromosome constraint fk_align_block_chromosome references chromosome(id) not null,
	start_pos number(10),
	stop_pos number(10),
	strand number(5),
	bit_string clob,
	gapped_length number(10),
	constraint align_block_pk primary key (chromosome, start_pos, stop_pos, id)
) organization index compress 1;
create index ix_align_block_alignment on align_block(alignment);
create index ix_align_block_id on align_block(id);

create index ix_align_block_chrpos on align_block(chromosome,start_pos,stop_pos);
create index ix_align_block_chrstop on align_block(chromosome,stop_pos);

create table ipmeta (
	id number(10) constraint fk_ipmeta_id unique not null,
	who varchar2(100), 	
	when date,
	antibody varchar2(500),
	xlink_condition varchar2(500),
	description varchar2(1000)
);

create table hybmeta (
	id constraint fk_hybmeta_id unique not null,
	ip1 constraint fk_hybmeta_ip1 references ipmeta(id) not null,
	ip2 constraint fk_hybmeta_ip2 references ipmeta(id) not null,
	who varchar2(100),
	when date,
	arraydesign constraint fk_hybmeta_galfile references galfiles(id),
	chip_id number(10),
	use_count number(10)
);

create table scanmeta (
	id number(10) constraint fk_scanmeta_id unique not null,
	hyb constraint fk_scanmeta_hyb references hybmeta(id) not null,
	experiment constraint fk_scanmeta_experiment references experiment(id) not null,
	who varchar2(100),
	when date,
	gprfilename varchar2(2000),
	tifffilename varchar2(2000)
);

create table orth_mapping (
	id number(10) constraint fk_orth_mapping_id unique not null,
	name varchar2(1000) not null,
	version varchar2(1000) not null
);

create table orth_pair (
	id number(10) constraint fk_orth_pair_id unique not null,
	mapping number(10) constraint fk_orth_pair_mapping references orth_mapping(id) not null,
	name1 varchar2(100) constraint fk_orth_pair_name1 not null,
	genome1 number(10) constraint fk_orth_pair_genome1 references genome(id) not null,
	name2 varchar2(100) constraint fk_orth_pair_name2 not null,
	genome2 number(10) constraint fk_orth_pair_genome2 references genome(id) not null
);
create index ix_orth_pair_mapping on orth_pair(mapping);
create index ix_orth_pair_name1 on orth_pair(genome1, name1);
create index ix_orth_pair_name2 on orth_pair(genome2, name2);

create table func_version (
	id number(10) constraint fk_func_version_id unique not null,
	name varchar2(500) constraint fk_func_version_name not null
);
create table func_category (
	id number(10) constraint fk_func_category_id unique not null,	
	version number(10) constraint fk_func_category_version references func_version(id) not null,
	name varchar2(100) constraint fk_func_category_name not null,
	description varchar2(1000) 
);
create index ix_func_category_version on func_category(version);

create table func_assignment (
	id number(10) constraint fk_func_assignment_id unique not null,
	version number(10) constraint fk_func_assignment_version references func_version(id) not null,
	object varchar2(100) constraint fk_func_assignment_object not null,
	category number(10) constraint fk_func_assignment_category references func_category(id) not null
);
create index ix_func_assignment_object on func_assignment(object);
create index ix_func_assignment_category on func_assignment(category);

create table func_subcategory (
	child_id number(10) constraint fk_func_sub_child references func_category(id) not null,
	parent_id number(10) constraint fk_func_sub_parent references func_category(id) not null,
	version number(10) constraint fk_func_sub_version references func_version(id) not null
);
create index ix_func_sub_child on func_subcategory(child_id);
create index ix_func_sub_parent on func_subcategory(parent_id);
create index ix_func_sub_version on func_subcategory(version);

CREATE TABLE ucsc_kgXref (
  kgID varchar(40) constraint kgXref_pk primary key,
  mRNA varchar(40),
  spID varchar(40) ,
  spDisplayID varchar(40) ,
  geneSymbol varchar(40) ,
  refseq varchar(40) ,
  protAcc varchar(40) ,
  description varchar(255))
organization index pctthreshold 20;
create index ix_kgXref_mrna on ucsc_kgXref(mRNA);
create index ix_kgXref_spid on ucsc_kgXref(spID);
create index ix_kgXref_dispid on ucsc_kgXref(spDisplayID);
create index ix_kgXref_symbol on ucsc_kgXref(geneSymbol);
create index ix_kgXref_refseq on ucsc_kgXref(refseq);
create index ix_kgXref_protacc on ucsc_kgXref(protAcc);


create table liver_conservation_motifs (
	chromosome constraint fk_liverconsmotifs_chromosome references chromosome(id) not null,
	factor varchar2(200),
	position number(10),
	score number(10,6),
	constraint liver_cons_motifs_pk primary key(factor,chromosome,position,score))
organization index compress 2;
