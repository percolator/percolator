#!/usr/bin/perl -w


use strict;

sub readSQT {
  my $fh = $_[0];
  my %all_hits;
  my @hits = ();
  my $charge = 0;
  my $mass = 0;
  my $id = "";
  while(<$fh>) {
    my @fields = split(/\s+/);
    if ($fields[0] eq 'S') {
      if (scalar(@hits)) {
        # move the deltaSC values upward one step...(weird and correct on the same time)
        foreach my $i (0 .. ($#hits-1)) {
          $hits[$i]->[2] = $hits[$i+1]->[2];
        }
	$all_hits{$id} = [@{$hits[0]}];
      }
      @hits = ();
      $id = $fields[3] . "_" . $fields[2];
      $charge=$fields[3];
      $mass=$fields[6];
    }
    if ($fields[0] eq 'M') { 
      push (@hits,[$fields[2],$mass-$fields[3],$fields[4],$fields[5],$fields[6],$fields[7]/$fields[8],$fields[9]]);
    }
  }
  return %all_hits;
}


sub isFullyTryptic {
  my $str = $_[0];
  return (($str =~ m/^[KR][.][^P]/g || $str =~ m/^-./g) && 
          ($str =~ m/[KR][.][^P]$/g || $str =~ m/[.]-$/g));
}

my $fFN = shift @ARGV;
my $rFN = shift @ARGV;

open(*FH,"< $fFN");
my %fHits = readSQT(\*FH);
close(*FH);

open(*FH,"< $rFN");
my %rHits = readSQT(\*FH);
close(*FH);

my ($label,@values,@LandS);

foreach my $psmId (keys %fHits) {
#  print $psmId . "\n";
#  print @{$fHits{$psmId}} . "\n";
  if (!defined($rHits{$psmId})){
      print $psmId;
      next;
  }
  if ($fHits{$psmId}->[3]<$rHits{$psmId}->[3]) {
    $label = -1;
    @values = @{$rHits{$psmId}};    
  } else {
    $label = 1;
    @values = @{$fHits{$psmId}};    
  }
#  print $values[1] . "\n";
  if (isFullyTryptic($values[6]) && abs($values[1])<=3) {
    push @LandS,[$values[3],$label];
    print $values[3] . "\n";
  }
}

my @sorted =  sort {$b->[0] <=> $a->[0]} @LandS;
my @labels = map {$_->[1]} @sorted;

open(*FH,"> beausoleil.res");
foreach $label (@labels) {
    print FH "$label\n";
}
close(*FH);
my $pos = scalar(grep {$_->[0]>1.4 && $_->[1]==1} @sorted); 
my $neg = scalar(grep {$_->[0]>1.4 && $_->[1]==-1} @sorted); 

my $fdr = $neg/$pos;

open(*FH,"> beausoleil.pnt");
print FH "$pos $fdr\n";
close(*FH);
