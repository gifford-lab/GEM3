package Bio::DB::ReadDBClient;

use warnings;
use strict;
use IO::Socket;
use Socket qw(IPPROTO_TCP TCP_NODELAY SO_LINGER);
use Config;
use Authen::SASL qw(Perl);
#use Authen::SASL;


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
  $self->{socket}->sockopt(SO_LINGER, 0);
  
  die "Couldn't authenticate to ${hostname}:${portnum}" unless($self->authenticate($hostname, $username, $passwd));

  $self->{request} = Bio::DB::ReadDBClient::Request->new();
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
#  print STDERR "PACKAGE IS " . __PACKAGE__ . "\n";
  my $sasl = Authen::SASL->new(mechanism => 'CRAM-MD5',
			       debug=>13,
			       user=>$username,
			       pass=>$passwd,
			       callback => {pass => $passwd,
					    user => $username}
			      );
#  print STDERR "Have user as " . $sasl->user() . "\n";
  my $client = $sasl->client_new($hostname, "readdb");
#  print STDERR "Have user as " . $client->_call('user') . "\n";
#  $client->property("POLICY_NOPLAINTEXT"=>1,
#		    "POLICY_NOANONYMOUS"=>1);
  my $response = $client->client_start();
  my $challenge = '';
  my $cont = 1;
  while ($cont) {
#    print STDERR "RESPONSE is $response\n";
    $self->sendString(length($response) . "\n" . $response);
    my $length = $self->readLine();
    $self->{socket}->recv($challenge, $length);
    if (length($challenge) != $length) {
      print STDERR "expected $length bytes as challenge but only read " . length($challenge) . "\n";
    }     
#    print STDERR "CHALLENGE is $challenge\n";
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
  my ($self) = @_;
  $self->sendString($self->{request}->toString());
  $self->{request}->clear();
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
  $self->{request}{requesttype} = 'shutdown';
  $self->sendRequest();
}

sub alignExists {
  my ($self, $alignid) = @_;
  $self->{request}{requesttype} = 'exists';
  $self->{request}{alignid} = $alignid;
  $self->sendRequest();
  my $response = $self->readLine();
  return ($response eq 'exists');
}

sub getChroms {
  my ($self, $alignid, $isPaired, $isLeft) = @_;
  $self->{request}{requesttype} = 'getchroms';
  $self->{request}{alignid} = $alignid;
  $self->{request}{ispaired} = $isPaired;
  $self->{request}{isleft} = $isLeft;
  $self->sendRequest();
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
  my ($self, $alignid, $isPaired, $isLeft, $isPlusStrand, @chroms) = @_;
  if (@chroms == 0) {
    @chroms = $self->getChroms($alignid);
  }
  my $count = 0;
  for my $chrom (@chroms) {
    $self->{request}{requesttype} = 'count';
    $self->{request}{alignid} = $alignid;
    $self->{request}{chromid} = $chrom;
    $self->{request}{ispaired} = $isPaired;
    $self->{request}{isleft} = $isLeft;
    $self->{request}{isplusstrand} = $isPlusStrand;
    $self->sendRequest();
    
    my $response = $self->readLine();
    if ($response ne "OK") {
      die "Can't get count for $alignid $chrom : $response";
    }
    $count += $self->readLine();
  }
  return $count;
}

sub getWeight {
  my ($self, $alignid, $isPaired, $isLeft, $isPlusStrand, @chroms) = @_;
  if (@chroms == 0) {
    @chroms = $self->getChroms($alignid);
  }
  my $weight = 0;
  for my $chrom (@chroms) {
    $self->{request}{requesttype} = 'weight';
    $self->{request}{alignid} = $alignid;
    $self->{request}{chromid} = $chrom;
    $self->{request}{ispaired} = $isPaired;
    $self->{request}{isleft} = $isLeft;
    $self->{request}{isplusstrand} = $isPlusStrand;
    $self->sendRequest();

    my $response = $self->readLine();
    if ($response ne "OK") {
      die "Can't get weight for $alignid $chrom : $response";
    }
    $weight += $self->readLine();
  }
  return $weight;
}

sub getCountRange {
  my ($self,$alignid,$chrom,$isPaired,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  $self->{request}{requesttype} = 'count';
  $self->{request}{alignid} = $alignid;
  $self->{request}{chromid} = $chrom;
  $self->{request}{ispaired} = $isPaired;
  $self->{request}{isleft} = $isLeft;
  $self->{request}{isplusstrand} = $isPlusStrand;
  $self->{request}{minweight} = $minweight;
  $self->{request}{start} = $start;
  $self->{request}{end} = $stop;
  $self->sendRequest();

  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get count for range for $alignid $chrom $start $stop $minweight: $response";
  }
  return $self->readLine();
}

sub getWeightRange {
  my ($self,$alignid,$chrom,$isPaired,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  $self->{request}{requesttype} = 'weight';
  $self->{request}{alignid} = $alignid;
  $self->{request}{chromid} = $chrom;
  $self->{request}{ispaired} = $isPaired;
  $self->{request}{isleft} = $isLeft;
  $self->{request}{isplusstrand} = $isPlusStrand;
  $self->{request}{minweight} = $minweight;
  $self->{request}{start} = $start;
  $self->{request}{end} = $stop;
  $self->sendRequest();

  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get weight for range for $alignid $chrom $start $stop $minweight: $response";
  }
  return $self->readLine();
}

sub getHits {
  my ($self,$alignid,$chrom,$isPaired,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  $self->{request}{requesttype} = 'gethits';
  $self->{request}{alignid} = $alignid;
  $self->{request}{chromid} = $chrom;
  $self->{request}{ispaired} = $isPaired;
  $self->{request}{isleft} = $isLeft;
  $self->{request}{isplusstrand} = $isPlusStrand;
  $self->{request}{minweight} = $minweight;
  $self->{request}{start} = $start;
  $self->{request}{end} = $stop;
  $self->{request}{map} = {wantweights=>1,
			   wantpositions=>1,
			   wantlengthsandstrands=>1};
  $self->sendRequest();
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get hits for $alignid $chrom $minweight: $response";
  }
  my $numhits = $self->readLine();
  my @output =  map {{position=>$_}} @{$self->readInts($numhits)};
  my $floats = $self->readFloats($numhits);
  for (my $i = 0; $i < @output; $i++) {
    $output[$i]{weight} = $floats->[$i];
  }
  my $ints = $self->readInts($numhits);
  for (my $i = 0; $i < @output; $i++) {
    my ($len,$strand) = decodeLAS($ints->[$i]);
    $output[$i]{length} = $len;
    $output[$i]{strand} = $strand;
  }
  for (my $i = 0; $i < @output; $i++) {
    my $hit = $output[$i];
  }

  return @output;
}

sub getHistogram {
  my ($self,$alignid,$chrom,$isPaired,$extension,$binsize, $start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  $self->{request}{requesttype} = 'histogram';
  $self->{request}{alignid} = $alignid;
  $self->{request}{chromid} = $chrom;
  $self->{request}{ispaired} = $isPaired;
  $self->{request}{isleft} = $isLeft;
  $self->{request}{isplusstrand} = $isPlusStrand;
  $self->{request}{minweight} = $minweight;
  $self->{request}{start} = $start;
  $self->{request}{end} = $stop;
  $self->{request}{map} = {binsize=>$binsize};
  if ($extension) {
    $self->{request}{map}{extension}=1;
  }
  $self->sendRequest();
  my $response = $self->readLine();

  if ($response ne "OK") {
    die "Can't get histogram for range for $alignid $chrom $start $stop $minweight : $response";
  }
  my $ints = $self->readInts($self->readLine());
#  print STDERR "Read histogram of length $#$ints\n";
  my $out = {};
  for (my $i = 0; $i < $#$ints; $i += 2) {
#    print STDERR "$ints->[$i] = $ints->[$i+1]\n";
    $out->{$ints->[$i]} = $ints->[$i+1];
  }
  return $out;
}

sub getWeightHistogram {
  my ($self,$alignid,$chrom,$isPaired,$extension,$binsize,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  $self->{request}{requesttype} = 'weighthistogram';
  $self->{request}{alignid} = $alignid;
  $self->{request}{chromid} = $chrom;
  $self->{request}{ispaired} = $isPaired;
  $self->{request}{isleft} = $isLeft;
  $self->{request}{isplusstrand} = $isPlusStrand;
  $self->{request}{minweight} = $minweight;
  $self->{request}{start} = $start;
  $self->{request}{end} = $stop;
  $self->{request}{map} = {binsize=>$binsize};
  if ($extension) {
    $self->{request}{map}{extension}=1;
  }
  $self->sendRequest();
  my $response = $self->readLine();
  if ($response ne "OK") {
    die "Can't get weight histogram for range for $alignid $chrom $start $stop $binsize $minweight : $response";
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

sub close {
  my ($self) = @_;
  $self->{request}{requesttype} = 'bye';
  $self->sendRequest();
  $self->{socket}->shutdown(2);
}

sub decodeLAS {
  my ($las) = @_;
  return ($las & 0x00007fff,
	  ($las & 0x00008000) ? 1 : -1);
}

sub sendInts {
  my ($self, $ints) = @_;
  return $self->send($ints,"i>");
}
sub sendFloats {
  my ($self, $floats) = @_;
  return $self->send($floats,"f>");
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
  return $self->read($numints, "i>");
}
sub readFloats {
  my ($self, $numints) = @_;
  return $self->read($numints, "f>");
}
sub read {
  my ($self, $numitems, $format) = @_;
  my @out = ();
  my $oldbuffer = '';
  while ($numitems > 0) {
    my $chunksize = (8192 < $numitems ? 8192 : $numitems);
    my $buffer;
    my $retval = $self->{socket}->recv($buffer, $chunksize * 4);
    if ($oldbuffer) {
      $buffer = $oldbuffer . $buffer;
    }
    unless (defined($retval)) {
      die "read error $retval : expected " . ($chunksize * 4) . " bytes but got " . length($buffer);
    }
    while (length($buffer) >= 4) {
      push(@out, unpack($format, substr($buffer, 0, 4, "")));
      $numitems--;
    }
    $oldbuffer = $buffer;
  }
  return \@out;
}

package Bio::DB::ReadDBClient::Request;

use strict;
use warnings;

my @boolfields = qw(ispaired isleft isplusstrand);
my @fields = qw(requesttype alignid chromid start end minweight);

sub new {
  my ($class) = @_;
  my $self = {map=>{}, list=>[]};
  bless $self,$class;
  return $self;
}

sub clear {
  my ($self) = @_;
  foreach (@fields) {
    delete $self->{$_};
  }
  $self->{map} = {};
  $self->{list} = [];
}

sub toString {
  my ($self) = @_;
  my $s = '';
  foreach (@fields) {
    $s .= exists $self->{$_} ? "$_=$self->{$_}\n" : "$_=\n";;
  }
  foreach (@boolfields) {
    next unless (exists $self->{$_} and defined $self->{$_});
    $s .= "$_=" . ($self->{$_} ? "true" : "false") . "\n";
  }
  foreach (keys %{$self->{map}}) {
    $s .= "$_=$self->{map}{$_}\n";
  }
  foreach (@{$self->{list}}) {
    $s .= "$_\n";
  }
  $s .= "ENDREQUEST\n";
  return $s;
}

1;

