#!/usr/bin/perl -w


use strict;

sub readSQT {
  my $fh = $_[0];
  my $label = $_[1];
  my $targetChar = $_[2];
  my @all_hits;
  my $charge = 0;
  while(<$fh>) {
    my @fields = split(/\s+/);
    if ($fields[0] eq 'S') {
      $charge=$fields[3];
    }
    if ($fields[0] eq 'M') {
      if ($targetChar==0 || $targetChar==$charge) {
        push (@all_hits,[$fields[5],$label]);
      }
    }
  }
  return @all_hits;
}


my $fFN = shift @ARGV;
my $rFN = shift @ARGV;
my $ch;
if ($#ARGV>-1) {
  $ch = shift @ARGV;
} else {
  $ch = 0;
}
print STDERR "Looking only at charge $ch\n" if $ch;
print STDERR "Looking at all charges\n" unless $ch;
open(*FH,"< $fFN");
my @normal = readSQT(\*FH,1,$ch);
close(*FH);

open(*FH,"< $rFN");
my @decoy = readSQT(\*FH,-1,$ch);
close(*FH);

my @both = @normal;
push @both,@decoy;

my @sorted =  sort {$b->[0] <=> $a->[0]} @both;
my @labels = map {$_->[1]} @sorted;

foreach my $label (@labels) {
    print "$label\n"
}
