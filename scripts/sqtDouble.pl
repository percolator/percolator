#!/usr/bin/perl -w


use strict;

sub similair {
  my $a = $_[0];
  my $b = $_[1];
  $a =~ s/I/L/g;
  $b =~ s/I/L/g;
  my (%aa,%bb) = ((),());
  my ($c,$diff1,$diff) = ("",0,0);
  foreach $c (split //,$a) {
    $aa{$c}++;
  }
  foreach $c (split //,$b) {
    $bb{$c}++;
  }
  foreach $c (keys %aa) {
    $diff1 += abs($aa{$c}-$bb{$c}) if defined($bb{$c});
    $diff1 += $aa{$c} unless defined($bb{$c});
  }
  foreach $c (keys %bb) {
    $diff += abs($aa{$c}-$bb{$c}) if defined($aa{$c});
    $diff += $bb{$c} unless defined($aa{$c});
  }
  $diff=$diff1>$diff?$diff1:$diff;
#  print "$diff ";
  return ($diff>3?0:1);
}

sub countDoubles {
  my @hits = @{$_[0]};
#  return (0,0) if @hits<2;
  print STDERR "Warning: file contains spectra with " . scalar(@hits) . "\n" if @hits > 3;
  my ($diffDoubles,$simDoubles)=(0,0);
  for(my $ix1=0;$ix1+1<@hits;$ix1++) {
    for(my $ix2=$ix1+1;$ix2<@hits;$ix2++) {
      if(similair($hits[$ix1],$hits[$ix2])) {
        $simDoubles++;
#        print "    similair " . $hits[$ix1] . " " . $hits[$ix2] . "\n";
      } else {
        $diffDoubles++;
#        print "not similair " . $hits[$ix1] . " " . $hits[$ix2] . "\n";
      }
    }
  }
  return ($diffDoubles,$simDoubles);
}

sub readSQT {
  my $fh = $_[0];
  my $qCut = $_[1];
  my $targetChar = $_[2];
  my @hits;
  my $charge = 0;
  my ($diffDoubles,$simDoubles)=(0,0);
  while(<$fh>) {
    my @fields = split(/\s+/);
    if ($fields[0] eq 'S') {
      my ($r1,$r2)=countDoubles(\@hits);
      if ($r1+$r2>1) { # This is ugly, we can not handle spectra with more than three hits in thus scheme
        if ($r1) {
          $r1--;
        } else {
          $r2--;
        }
      }
      $diffDoubles += $r1;
      $simDoubles += $r2;
      $charge=$fields[3];
      @hits = ();
    }
    if ($fields[0] eq 'M') {
      if (((not $targetChar) || $targetChar==$charge) && -$fields[6] < $qCut) {
#        print -$fields[6] . " " . $qCut . " " . $fields[9] . "\n";
        my $pep = $fields[9];
        $pep =~ s/^.[.]//g;
        $pep =~ s/[.].$//g;
        push (@hits,$pep);
      }
    }
  }
  print "Found " . $diffDoubles . " different doubles and " . $simDoubles . " similair doubles\n";
}


my $fFN = shift @ARGV;
my $qCut = shift @ARGV || 0.01;
my $ch = shift @ARGV || 0;
print STDERR "Looking only at charge $ch\n" if $ch;
print STDERR "Looking at all charges\n" unless $ch;
open(*FH,"< $fFN");
readSQT(\*FH,$qCut,$ch);
close(*FH);
