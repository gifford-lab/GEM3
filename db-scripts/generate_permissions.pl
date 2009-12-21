#!/usr/bin/perl

while (<STDIN>) {
  if (/create (sequence|table)\s(\S+[^\;\s])/) {
    print "grant all on $2 to psrgwrite;\n";
    print "grant select on $2 to psrgread;\n";
  }

}
