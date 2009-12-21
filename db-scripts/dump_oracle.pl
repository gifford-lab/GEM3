#!/usr/bin/perl

=head1 dunldr

  unload data from an oracle database

  use 'dunldr -help' for help on usage

  jared still
  10/24/2001

=cut

use warnings;
use FileHandle;
use DBI;
use strict;
use File::Path;
use IO::File;
use Data::Dumper;

use Getopt::Long;

our %optctl = ();
our %bincol = ();
our %hexcols = ();


unless (Getopt::Long::GetOptions( \%optctl,				  
				  "database=s",
				  "username=s",
				  "password=s",
				  "owner=s",
				  "directory=s",
				  "dateformat=s",
				  "header!",
				  "schemadump!",
				  "longlen=i",
				  "rowlimit=i",
				  "table=s@",
				  "bincol=s" => \%bincol,
				  "sysdba!",
				  "sysoper!",
				  "z","h","help")) {
  Usage(1); 
}

for my $table ( keys %bincol ) {
  my @bincols = split(/\,/,$bincol{$table});
  $hexcols{uc($table)} = \@bincols;
  
}

#print Dumper(\%optctl);
#print Dumper(\%hexcols);
#for my $hexdumpcol ( @{$hexcols{XML_DATA}} ) {
#print "hexdumpcol: $hexdumpcol\n";
#}
#exit;


our($db, $username, $password, $connectionMode);

$connectionMode = 0;

if ( $optctl{sysoper} ) {
  $connectionMode = 4;
} if ( $optctl{sysdba} ) {
  $connectionMode = 2;
}

Usage(1) unless $optctl{database};
Usage(1) unless $optctl{username};
Usage(1) unless $optctl{password};
Usage(1) unless $optctl{owner};


$optctl{longlen} = 65535 unless $optctl{longlen};

if ( $optctl{h} || $optctl{z} || $optctl{help} ) {
  Usage(0);
}

if ( $optctl{schemadump} ) {
  $optctl{table} = ['SCHEMADUMP'];

} else {
  Usage(1) unless $optctl{table};
}


# default hdr to off
$optctl{header} ||= 0;

#if ( $optctl{bincol} ) {
#}

$username=$optctl{username};
$password = $optctl{password};
$db = $optctl{database};

# create the working directory
unless ( $optctl{directory} ) {
  $optctl{directory} = qq{$optctl{owner}.dump};
}

# create directory path if it doesn't exist 
-d $optctl{directory} || File::Path::mkpath([$optctl{directory}]); 

our $dbh = DBI->connect(
			'dbi:Oracle:' . $db,
			$username, $password,

			{
			 RaiseError => 1,
			 AutoCommit => 0,
			 ora_session_mode => $connectionMode
		       }
		       );

die "Connect to $db failed \n" unless $dbh;

$dbh->{LongReadLen} = $optctl{longlen};

# set Oracle NLS date format
if ( $optctl{dateformat} ) {
  $dbh->do(qq{alter session set nls_date_format
	      = '$optctl{dateformat}'} );
}

my $tableHash = new Tables($dbh, \%optctl);

#print "tables: ", join(':', keys %{$tableHash}), "\n";
#for my $table ( keys %{$tableHash} ){
#print "TABLE: $table FILE: $tableHash->{$table}\n";
#}


# print console info immediately
autoflush STDOUT 1;

my $sth;

# take a dump
for my $table ( keys %{$tableHash} ) {

  print "Table: $table\n";

  my $sql = qq{select * from $optctl{owner}\.$table};

  if ( $optctl{rowlimit}) {
    $sql .= qq{ where rownum <= $optctl{rowlimit}};
  }

  $sth = $dbh->prepare($sql);

  my @columns = @{$sth->{NAME_uc}};
  my %colOrder = ();
  for my $el ( 0 ..$#columns ) {
    $colOrder{$columns[$el]} = $el;

  }

  my $dumpFile = $optctl{directory} . '/' . $tableHash->{$table}; open(DUMP, "+> $dumpFile") || die "could not create file $dumpFile - $!\n";

  if ( $optctl{header} ) {
    print DUMP join(',',@columns),"\n";
  }

  $sth->execute;

  # create the ctl and par files
  Tables->createCtl(
		    TABLE => $table,
		    COLUMNS => \@columns,
		    DUMPFILE => $tableHash->{$table},
		    DIRECTORY => $optctl{directory},
		    SCHEMA => $optctl{owner},
		    HEXCOLS => \@{$hexcols{$table}},
		    COLORDER => \%colOrder
		   );

  # turn warnings off here so that warnings are not
  # reported for null columns when printed
  # comment it out to see what I mean


  no warnings;
  while ( my $ary = $sth->fetchrow_arrayref ) {	# change column to hex if specified as binary via -bincol arg
    if ( exists $hexcols{$table} ) {
      for my $hexdumpcol ( @{$hexcols{$table}} ) {
	$ary->[$colOrder{uc($hexdumpcol)}] =

	  uc(unpack("H*",$ary->[$colOrder{uc($hexdumpcol)}]));
      }
    }
    print DUMP q{"} . join(q{","},@{$ary}) . qq{"\n}; #print "ROW: " . q{'} . join(q{','},@{$ary}) . qq{'\n}; 
  }
  use warnings;
  close DUMP;
}
$sth->finish;
$dbh->disconnect;

sub Usage {

  my ($exitCode) = @_;

  print q{

	  dunldr - data unloader for Oracle

	usage:

	  dunldr -database <database> -username <userid> -password <password> \ 
	  -directory <data unload directory> \
	  -header|noheader \
	  -owner <schema owner> \
	  -table <table1,table2,table3,...)



    -database database name

      -username user to login as

	-password password for login user

	  -owner owner of tables to dump

	    -directory directory to unload data into will default to <owner>.dump

	      -dateformat Oracle NLS date format - optional -header|noheader should first line include column names?

		-table table to dump. may be repeated as many times as necessary.

		  -schemadump dump entire schema of <owner> will ignore -table settings

		    -rowlimit limit number of rows returned

		      -longlen if longs are in the table, set this to the maximum length you want.
			defaults to 65535

			  -bincol use to specify columns that should be dumped in hex format. columns with binary data tend to cause problems in text dumps.
			    e.g. -bincol <table_name>=<column_name,column_name,...>

			      dunldr -database orcl -username system -password manager \

				-owner scott -directory scott.tables \
				  -header \
				    -table emp \
				      -table dept \
					-table sales


					  dunldr -database orcl -username system -password manager \

					    -owner scott \
					      -dateformat 'mm/dd/yyyy' \
						-header \
						  -schemadump \
						    -bincol xml_data=payload,header,authorization \
						      -bincol app_notes=text



						    }
;

exit $exitCode ? $exitCode : 0;
}

package Tables;

sub new {

  my $pkg = shift;
  my $class = ref($pkg) || $pkg;

  my ( $dbh, $optionHash ) = @_;

  my $tableHash;
  if ( grep(/^SCHEMADUMP$/, @{$optionHash->{table}} ) ) { # get all tables of owner
    my $sql = q{
		select table_name
		from all_tables
		where owner = ?
	      };
    my $sth = $dbh->prepare($sql);
    $sth->execute(uc($optionHash->{owner}));
    my @tableArray;
    while ( my $ary = $sth->fetchrow_arrayref ) {
      push(@tableArray, $ary->[0]);
    }
    $tableHash = setTables(\@tableArray);

  } else {
    $tableHash = setTables(\@{$optionHash->{table}});
  }

  bless $tableHash, $class;
  return $tableHash;

}

=head1 setTables

  make a neat hash of the form TABLE_NAME => 'table_name.dump' all table names upper case, all file names lower case for dump file names - Perl is awesome

=cut

  sub setTables {
    my ($tableArray) = shift;

    my %tables = map(
		     split(/:/, $_),
		     map(
			 $_.':'.lc($_).'.txt',

			 split(
			       /:/,
			       uc(join(':',@{$tableArray}))

			      )
			)
		    );

    # uncomment these lines to see it
    #use Data::Dumper;


    #print Dumper(\%tables);
    #exit;

    my $hashRef = \%tables; 
    return $hashRef;
  }

sub createCtl {
  my($self,%args) = @_;

  my @columns = @{$args{COLUMNS}};
  my %colOrder = %{$args{COLORDER}};

  if ( $args{HEXCOLS} ) {
    for my $hexdumpcol ( @{$args{HEXCOLS}} ) {
      $columns[$colOrder{uc($hexdumpcol)}] =
	$columns[$colOrder{uc($hexdumpcol)}] .

	  qq{ "hex_to_raw(:$columns[$colOrder{uc($hexdumpcol)}])"};
    }
  }

  my $ctlFile = $args{DIRECTORY}. '/' . lc($args{TABLE}) . '.ctl'; my $ctlFh = new IO::File();
  $ctlFh->open("> $ctlFile") || die "cannot create file $ctlFile - $!\n";
  $ctlFh->print("load data\n");
  $ctlFh->print("infile '$args{DUMPFILE}'\n");
  $ctlFh->print("into table $args{TABLE}\n");
  $ctlFh->print(q{fields terminated by ',' optionally enclosed by
		  '"'}. "\n");
  $ctlFh->print("(\n");
  $ctlFh->print( "\t" . join(",\n\t",@columns) . "\n");
  $ctlFh->print(")\n");
  $ctlFh->close;

  my $parFile = $args{DIRECTORY}. '/' . lc($args{TABLE}) . '.par'; my $parFh = new IO::File();
  $parFh->open("> $parFile") || die "cannot create file $parFile - $!\n";
  $parFh->print("userid = $args{SCHEMA}\n");
  $parFh->print("control = " . lc($args{TABLE}) . ".ctl\n");
  $parFh->print("log = " . lc($args{TABLE}) . ".log\n");
  $parFh->print("bad = " . lc($args{TABLE}) . ".bad\n");
  $parFh->close;

} 
