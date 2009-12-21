#!/usr/bin/perl

# Dumps an experiment from the database, runs JBD, saves the results back to the database, and runs the default binding scans
# Does all this on the cluster.
# 
# Input is on STDIN and is an experiment name.  Optionally, it's "experiment name;experiment version;replicate".
# This uses (and creates if it doesn't exist) the directory ~/jbd_files as a temporary working area
#
# If you're using a job-queuing system other than Torque, then you'll need to replace the
# calls to PBS::Client with something else.  One nice feature of PBS::Client are job dependencies- the
# jobs that run JBD are submitted first, a job that loads the JBD results depends on them, and finally a job 
# that runs the binding scan depends on the job that loads the results.
#
# Usage:
# cat exptlist.txt | ./run_jbd.pl --species "$SC" --genome "SGDv1" \
#         [--queue batch] [--version "JBD results version"] [--epargs "--ep_cmdline_arg1 --ep_cmdline_arg2 val"] \
#         [--generateargs "--generate_datafiles.pl_arg1"]
#
# --epargs and --generateargs get passed through to ep_cmdline.pl and generate_datafiles.pl


use strict;
use warnings;
use PBS::Client;
use Getopt::Long;

my $codedir = $0;
$codedir =~ s/run_jbd.pl//;

my $client = new PBS::Client();
my $jbdfiles = ($ENV{HOME}  || '/afs/csail.mit.edu/group/psrg/tmp'). '/jbd_files';
unless (-e $jbdfiles and -d $jbdfiles) {
  mkdir($jbdfiles);
}
chdir($jbdfiles) or die "Can't change to $jbdfiles : $!";
# queue to which jobs are submitted
my $queue = 'filler';  
my $baseversion = '5/27/08, default params';
my ($epargs,$generateargs);
my ($species,$genome);
GetOptions("queue=s"=>\$queue,
	   "version=s"=>\$baseversion,
	   "epargs=s"=>\$epargs,
	   "generateargs=s"=>\$generateargs,
	   "species=s"=>\$species,
	   "genome=s"=>\$genome);

my $roleprefix = '';
if ($ENV{CHIPCHIPROLE}) {
  $roleprefix = "export CHIPCHIPROLE=$ENV{CHIPCHIPROLE} && ";
}
if ($species =~ /;/) {
    ($species,$genome) = split(';',$species);
}
if ($genome =~ /;/) {
    ($species,$genome) = split(';',$genome);
}
print STDERR "$species ::: $genome\n";

while (<STDIN>) {
  chomp;
  my ($exptname,$exptversion,$replicate) = split(/;/,$_);  
  $exptversion ||= 'median linefit';
  my $underscore = $exptname . '__' . $exptversion;
  my $version = $baseversion;
  if ($replicate) {
    $underscore .= "__${replicate}";
    $version = "${baseversion}: ${replicate}";
  }
  $underscore =~ s/\W/_/g;
  if (-e $underscore) {
    my $count = 1;
    while (-e "${underscore}__${count}") {
      $count++;
    }
    $underscore = "${underscore}__${count}";
  }
  unless ($genome) {
      ($species,$genome) = abbrev_to_species($_);
  }
  my $spec = "${species};${genome}";
  
  my $repstring = $replicate ? ";${replicate}" : '';
  print STDERR "Underscore is $underscore    $jbdfiles\n";  
  mkdir("${jbdfiles}/${underscore}");
  next unless (-d "${jbdfiles}/${underscore}");
  my $cmd = "cd ${jbdfiles}/${underscore} && java edu.mit.csail.cgs.tools.chipchip.GenerateJBDInput --expt \"${exptname};${exptversion}${repstring}\" --species \"$spec\" ${generateargs} --outname $underscore";
  print STDERR "cmd : ${cmd}\n";
  system("${cmd}");
  $cmd = "cd ${jbdfiles}/${underscore} && /afs/csail.mit.edu/group/psrg/software/jbd-1.2/ep_cmdline.pl --dir . --dontsub  ${epargs} @ARGV";
  print STDERR "cmd : ${cmd}\n";
  system("${cmd}");
  opendir(DIR,"${jbdfiles}/${underscore}") or die "Can't open ${jbdfiles}/${underscore} : $!";
  my @files = readdir(DIR);
  closedir(DIR);
  my @jobs = ();
  chdir("${jbdfiles}/${underscore}");
  foreach (@files) {
    if (/^jbd\d+.sh/) {
      chmod(0755,"${jbdfiles}/${underscore}/${_}");
      my $jbdjob = new PBS::Client::Job(queue=>$queue,
					cmd=>"${roleprefix} ${jbdfiles}/${underscore}/${_}");
      push(@jobs,$jbdjob);
    }
  }  
  my $loadjob = new PBS::Client::Job(queue=>$queue,
				     cmd=>"${roleprefix} cd ${jbdfiles}/${underscore} && cat args.txt | xargs java edu.mit.csail.cgs.tools.chipchip.AddJBDResults --analysis '$exptname;${version}' --paramfile '${jbdfiles}/${underscore}/bayes.params' -- *.jbd");
  $loadjob->prev({ok=>\@jobs});
  my $scanjob = new PBS::Client::Job(queue=>$queue,
				     cmd=>"${roleprefix} java -classpath /afs/csail.mit.edu/group/psrg/jar/psrg/gse_tools-1.0.0.jar edu.mit.csail.cgs.tools.binding.ScanTool --species \"${species}\" --genome \"${genome}\" --type jbd --expt \"${exptname}\" --version \"${version}\"");
  $scanjob->prev({ok=>$loadjob});
  # only need to submit the last job in the chain of dependencies and the rest
  # get pushed in with it.
  $client->qsub($scanjob);
  chdir($jbdfiles);
}


sub abbrev_to_species {
  my ($abbrev) = @_;
  if ($abbrev =~ /^Mm/) {
    return ("Mus musculus","mm8");
  } elsif ($abbrev =~ /^Hs/) {
    return ("Homo sapiens","hg18");
  } elsif ($abbrev =~ /^Dm/) {
    return ("Drosophila melanogaster","dmel2");
  } elsif ($abbrev =~ /^Sc/) {
    return ("Saccharomyces cerevisiae","sacCer1");
  } elsif ($abbrev =~ /^Ce/) {
    return ("Caenorhabditis elegans","ce2");
  } else {
    die "Unknown abbreviation $abbrev";
  }
}
