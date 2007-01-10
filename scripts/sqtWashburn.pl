#!/usr/bin/perl -w

use strict;

sub isTrypN {
    my $seq = shift;
    return (($seq =~ m/^[A-Z]*[KR][.][^P]/g) || ($seq =~ m/^[*][.]/g))?1:0;
}

sub isTrypC {
    my $seq = shift;
    return (($seq =~ m/[KR][.][^P][A-Z]*$/g) || ($seq =~ m/[.][*]$/g))?1:0;
}

my $fullRecords = shift @ARGV || 0;
my $found=0;

sub reportLoci {
  if ($fullRecords) {
    print $_[1];
  } else {
    my @loci = @{$_[0]};
    foreach my $l (@loci) { print "$l\n";}
  }
  $found++;
}

my @hits = ();
my @locus = ();
my $charge = 0;
my $mass = 0;
my $id = "";
my $m =0;
my $lines="";
while(my $line=<>) {
  my @fields = split(/\s+/,$line);
  if ($fields[0] eq 'S') {
    if (scalar(@hits) && $hits[0]>=0.1) {
      if ($charge==1 && $hits[1]>=1.9 && $hits[2]==1 && $hits[3]==1)
      { reportLoci(\@locus,$lines)}
      if ($charge==2 && ($hits[1]>=3.0 || ($hits[1]>=2.2 && ($hits[2]==1 || $hits[3]==1))))
      { reportLoci(\@locus,$lines)}
      if ($charge==3 && $hits[1]>=3.75 && ($hits[2]==1 || $hits[3]==1))
      { reportLoci(\@locus,$lines)}
    }
    @hits = ();
    @locus = ();
    $m=0;
    $charge=$fields[3];
    $lines="";
  }
  $lines .= $line;
  if ($fields[0] eq 'M') {
    if ($m==0) { @hits = (-1, $fields[5], isTrypN($fields[9]), isTrypC($fields[9]))};
    if ($m==1) { $hits[0] = $fields[4] };
    $m++;
  }
  if ($fields[0] eq 'L' && $m==1) {
    push @locus,$fields[1];
  }
}
print STDERR "Found $found records\n";



