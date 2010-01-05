package ReadDBClient;

use warnings;
use strict;
use IO::Socket;
use Config;
#use Authen::SASL qw(Perl);
use Authen::SASL;


sub new {
  my ($class, $hostname, $portnum, $username, $passwd) = @_;
  my $self = {debug=>0};
  bless $self,$class;
  return $self->init($hostname, $portnum, $username, $passwd);
}

sub init {
  my ($self,$hostname,$portnum,$username,$passwd) = @_;
  unless ($hostname && $portnum && $username && $passwd) {
    my %props = $self->read_config_file();
    ($hostname,$portnum,$username,$passwd) = @props{qw(hostname port username passwd)};
  }
  
  $self->{socket} = IO::Socket::INET->new(PeerAddr=>$hostname,
					  PeerPort=>$portnum,
					  Proto=>'tcp') or die "can't open socket to ${hostname}:${portnum} : $!";
  
  die "Couldn't authenticate to ${hostname}:${portnum}" unless($self->authenticate($hostname, $username, $passwd));

  $self->sendRequest("byteorder", ($Config{byteorder} =~ /^1/) ? 'little' : 'big');
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Couldn't set byte order : $response";
  }
  return $self;
}

# returns a hash of the key-value pairs from the config file, by default ~/.readdb_passwd
sub read_config_file {
  my ($self) = @_;
  my $homedir = $ENV{'HOME'};
  my $basename = 'readdb_passwd';
  if ($ENV{'READDBROLE'}) {
    $basename = $ENV{'READDBROLE'} . $basename;
  }
  open(PROPS,"${homedir}/.${basename}") or die "can't open config file .$basename in $homedir : $!";
  my %props = ();
  while(<PROPS>) {
    chomp;
    s/^\s*//;
    s/\s*$//;
    my ($k,$v) = ($_ =~ /([^=]*)\s*=\s*(.*)/);
    next unless ($k);
    $props{$k} = $v;
  }
  return %props;
}

sub debug {
  my ($self,$opt) = @_;
  if (defined($opt)) {
    $self->{debug} = $opt;
  }
  return $self->{debug};
}

sub authenticate {
  my ($self, $hostname, $username, $passwd) = @_;
  $self->sendString("${username}\n");
  my $sasl = Authen::SASL->new(
			       mechanism => 'CRAM-MD5',
			       callback => {
					    pass => $passwd,
					    user => $username,
					   }
			      );
  my $client = $sasl->client_new($hostname, "readdb");
  $client->property("POLICY_NOPLAINTEXT"=>1,
		    "POLICY_NOANONYMOUS"=>1);
  my $response = $client->client_start();
  my $challenge = '';
  my $cont = 1;
  while ($cont) {
    $self->sendString(length($response) . "\n" . $response);
    my $length = $self->readLine();
    $self->{socket}->recv($challenge, $length);
    if (length($challenge) != $length) {
      print STDERR "expected $length bytes as challenge but only read " . length($challenge) . "\n";
    } 
    $self->{socket}->recv($cont, 1);
    $response = $client->client_step($challenge);
    $cont = unpack("c",$cont);
  }
  my $status = $self->readLine();
  if ($status ne "authenticated as ${username}") {
    print STDERR "SASL error " . $sasl->error() . "\n";
    return 0;
  } else {
    return 1;
  }
}

sub sendString {
  my ($self, $string) = @_;
  if ($self->{debug}) {
    print STDERR "SEND $string";
  }
  $self->{socket}->send($string);
  $self->{socket}->flush();
}
sub sendRequest {
  my ($self,@args) = @_;
  $self->sendString(join("\n",@args) . "\nENDREQUEST\n");
}

sub readLine {
  my ($self, $string) = @_;
  my $line = '';
  my $char = ' ';
  while (1) {
    $self->{socket}->recv($char, 1);
    if (length($char) > 0) {
      if ($char eq "\n") {
	if ($self->{debug}) {
	  print STDERR "READ $line\n";
	}
	return $line;
      } else {
#	print STDERR "read char $char\n";
	$line .= $char;
      }
    } else {
      die "read error after $line";
    }
  }
}

sub shutdown {
  my ($self) = @_;
  $self->sendRequest("shutdown");
}

sub store {
  my ($self, $alignid, $chromid, $hits, $weights) = @_;
  $self->sendRequest("store", $alignid, $chromid, ($#$hits+1));
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't setup store hits : $response";
  }
  $self->sendInts($hits);
  $self->sendFloats($weights);
  $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't store hits : $response";
  }
}

sub alignExists {
  my ($self, $alignid) = @_;
  $self->sendRequest("exists",$alignid);
  my $response = $self->readLine();
  return ($response eq 'exists');
}

sub deleteAlignment {
  my ($self, $alignid) = @_;
  $self->sendRequest("deletealign",$alignid);
  my $response = $self->readLine();
  die "Can't delete alignment $alignid : $response" unless ($response eq 'OK');
}

sub getChroms {
  my ($self, $alignid) = @_;
  $self->sendRequest("getchroms",$alignid);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get chroms for $alignid : $response";
  }
  my $numchroms = $self->readLine();
  my @chroms = ();
  while ($numchroms-- > 0) {
    push(@chroms, $self->readLine());
  }
  return @chroms;
}

sub getCount {
  my ($self, $alignid, @chroms) = @_;
  if (@chroms == 0) {
    @chroms = $self->getChroms($alignid);
  }
  my $count = 0;
  for my $chrom (@chroms) {
    $self->sendRequest("count",$alignid, $chrom);
    my $response = $self->readLine();
    if ($response ne "OK") {
      die "Can't get count for $alignid $chrom : $response";
    }
    $count += $self->readLine();
  }
  return $count;
}

sub getWeight {
  my ($self, $alignid, @chroms) = @_;
  if (@chroms == 0) {
    @chroms = $self->getChroms($alignid);
  }
  my $weight = 0;
  for my $chrom (@chroms) {
    $self->sendRequest("weight",$alignid, $chrom);
    my $response = $self->readLine();
    if ($response ne "OK") {
      die "Can't get weight for $alignid $chrom : $response";
    }
    $weight += $self->readLine();
  }
  return $weight;
}

sub getCountRange {
  my ($self,$alignid,$chrom,$start,$stop,$minweight) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $self->sendRequest("countrange",$alignid, $chrom, $start, $stop, $minweight);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get count for range for $alignid $chrom $start $stop $minweight: $response";
  }
  return $self->readLine();
}

sub getWeightRange {
  my ($self,$alignid,$chrom,$start,$stop,$minweight) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $self->sendRequest("weightrange",$alignid, $chrom, $start, $stop, $minweight);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get weight for range for $alignid $chrom $start $stop $minweight: $response";
  }
  return $self->readLine();
}

sub getHits {
  my ($self,$alignid,$chrom,$minweight) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $self->sendRequest("gethits",$alignid, $chrom, $minweight);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get hits for $alignid $chrom $minweight: $response";
  }
  my $numhits = $self->readLine();
  return $self->readInts($numhits);
}

sub getWeights {
  my ($self,$alignid,$chrom,$minweight) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $self->sendRequest("getweights",$alignid, $chrom, $minweight);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get weights for $alignid $chrom $minweight: $response";
  }
  my $numhits = $self->readLine();
  return $self->readFloats($numhits);
}

sub getHitsRange {
  my ($self,$alignid,$chrom,$start,$stop,$minweight) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $self->sendRequest("gethitsrange",$alignid, $chrom, $start, $stop, $minweight);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get hits for range for $alignid $chrom $start $stop $minweight: $response";
  }
  my $numhits = $self->readLine();
  return $self->readInts($numhits);
}

sub getWeightsRange {
  my ($self,$alignid,$chrom,$start,$stop,$minweight) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $self->sendRequest("getweightsrange",$alignid, $chrom, $start, $stop, $minweight);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get weights for range for $alignid $chrom $start $stop $minweight: $response";
  }
  my $numhits = $self->readLine();
  return $self->readFloats($numhits);
}

sub getHistogram {
  my ($self, $alignid, $chrom, $start, $stop, $binsize, $minweight, $readextension) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $readextension ||= 0;
  $self->sendRequest("histogram",$alignid, $chrom, $start, $stop, $binsize, $minweight, $readextension);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get histogram for range for $alignid $chrom $start $stop $binsize $minweight $readextension : $response";
  }
  my $ints = $self->readInts($self->readLine());
  my $out = {};
  for (my $i = 0; $i < $#$ints; $i += 2) {
    $out->{$ints->[$i]} = $ints->[$i+1];
  }
  return $out;
}

sub getWeightHistogram {
  my ($self, $alignid, $chrom, $start, $stop, $binsize, $minweight, $readextension) = @_;
  $minweight = 'NaN' unless (defined($minweight));
  $readextension ||= 0;
  $self->sendRequest("weighthistogram",$alignid, $chrom, $start, $stop, $binsize, $minweight, $readextension);
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get weight histogram for range for $alignid $chrom $start $stop $binsize $minweight $readextension : $response";
  }
  my $len = $self->readLine();
  my $positions = $self->readInts($len);
  my $values = $self->readFloats($len);
  my $out = {};
  for (my $i = 0; $i < $#$positions; $i++) {
    $out->{$positions->[$i]} = $values->[$i];
  }
  return $out;
}

sub getACL {
  my ($self, $alignid) = @_;
  $self->sendRequest("getacl",$alignid);
  my $response = $self->readLine();
  if ($response ne 'OK') {
    die "Can't get acl for $alignid : $response";
  }
  my $acl = {};
  for (my $i = 0; $i < 3; $i++) {
    my $type = $self->readLine();
    my $num = $self->readLine();
    while ($num-- > 0) {
      push(@{$acl->{$type}}, $self->readLine());
    }
  }
  return $acl;
}

sub close {
  my ($self) = @_;
  $self->sendRequest("bye");
  $self->{socket}->shutdown();
}

sub sendInts {
  my ($self, $ints) = @_;
  return $self->send($ints,"i");
}
sub sendFloats {
  my ($self, $floats) = @_;
  return $self->send($floats,"f");
}
sub send {
  my ($self, $array, $format) = @_;
  my $chunksize = 8192;
  my $i = 0;
  while ($i < $#$array) {
    my $buffer = '';
    for (my $j = 0; $j < $chunksize && ($i + $j < $#$array); $j++) {
      $buffer .= pack($format, $array->[$i+$j]);
    }
    $self->{socket}->send($buffer);
    $i += $chunksize;
  }
}

# returns an arrayref 
sub readInts {
  my ($self, $numints) = @_;
  return $self->read($numints, "i");
}
sub readFloats {
  my ($self, $numints) = @_;
  return $self->read($numints, "f");
}
sub read {
  my ($self, $numitems, $format) = @_;
  my @out = ();
  while ($numitems > 0) {
    my $chunksize = (8192 < $numitems ? 8192 : $numitems);
    my $buffer;
    my $data = $self->{socket}->recv($buffer, $chunksize * 4);
    unless (length($buffer) == $chunksize * 4) {
      die "read error";
    }
    while ($buffer) {
      push(@out, unpack($format, substr($buffer, 0, 4, "")));
    }
    $numitems -= $chunksize;
  }
  return \@out;
}

1;
