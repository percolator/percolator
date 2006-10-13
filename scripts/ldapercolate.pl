#!/usr/bin/perl -w

use Math::Matrix;

use strict;

my $templ = new Math::Matrix([0]);

my $ONLY_CHARGE = 2;
my $NITER = 1;
my $forwardFN = '/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/forward/050606-pps-2-01-forward-no-norm.sqt';
my $randomFN = '/var/noble/data/tandem-ms/maccoss/rphase/2006-05-06/random/050606-pps-2-01-random-no-norm.sqt';

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
      if ($charge eq $ONLY_CHARGE && scalar(@hits)) {
        # move the deltaSC values upward one step...(weird and correct on the same time)
        foreach my $i (0 .. ($#hits-1)) {
          $hits[$i]->[2] = $hits[$i+1]->[2];
        }
#        $#hits--; # Last element cant be shifted so we remove that one
#        $all_hits{$id} = [@hits];
# New way: only look at the best...
	$all_hits{$id} = [$hits[0]];
      }
      @hits = ();
      $id = $fields[3] . "_" . $fields[2];
      $charge=$fields[3];
      $mass=$fields[6];
    }
    if ($fields[0] eq 'M') { # dont do that && $fields[4]>0) { # only look for second best deltasc>0
#      print;
      push (@hits,[$fields[2],$mass-$fields[3],$fields[4],$fields[5],$fields[6],$fields[7]/$fields[8]]);
    }
  }
#  print "Read $#all_hits spectras\n";
  return %all_hits;
}

sub getXcorr {
  return $_[0]->[3];
}

my @w_values = (0,0,0,1,0,0);
sub getDisc {
  my @vec=@{$_[0]};
  my $sum=0;
  foreach my $n (0 .. $#vec) {
    $sum += $vec[$n]*$w_values[$n];
  }
  return $sum;
}

#my $calcScore = \&getXcorr;
my $calcScore = \&getDisc;

sub allTheBest {
  my @allHits = @{$_[0]};
  my @bestHits = ();
  my @bestScores = ();
  foreach my $refHits (@allHits) {
    my $highScore = 0;
#    my $highScoreV = [];
    my $highScoreV;
#    foreach my $refValues (@$refHits) {
    my $refValues = $refHits->[0];
    my $score=&$calcScore($refValues);
#      print "$score\n";
#      ($highScore,$highScoreV) = ($score,$refValues) if ($score>$highScore);
#    }
#    push(@bestHits,$highScoreV) if ($highScore>0);
#    push(@bestScores,$highScore) if ($highScore>0);
    push(@bestHits,$refValues);
    push(@bestScores,$score);
  }
  my @sortedScores = sort { $b <=> $a } @bestScores;
  return (\@bestHits, \@sortedScores);
}

sub printScores {
  my @bestScoreV = @{$_[0]};
#  my $pre = $_[1];
  my $fname = $_[1];
  open(FH,$fname);
  foreach my $rScoreV (@bestScoreV) {
#    print "$pre " . $rScoreV->[1] . " " . $rScoreV->[2] . "\n" if (ref($rScoreV) eq 'ARRAY');
    print FH $rScoreV->[2] . " " . $rScoreV->[3] . "\n";
  }
  close(FH);
}

sub printPyML2 {
  my @fScoreV = @{$_[0]};
  my @rScoreV = @{$_[1]};
#  my $pre = $_[1];
  my $fname = $_[2];
  open(FH,$fname);
  foreach my $aScoreV (@fScoreV) {
    print FH $aScoreV->[2] . " " . $aScoreV->[3] . " 1\n";
  }
  foreach my $aScoreV (@rScoreV) {
    print FH $aScoreV->[2] . " " . $aScoreV->[3] . " 0\n";
  }
  close(FH);
}
sub printPyML6 {
  my @fScoreV = @{$_[0]};
  my @rScoreV = @{$_[1]};
#  my $pre = $_[1];
  my $fname = $_[2];
  open(FH,$fname);
  foreach my $aScoreV (@fScoreV) {
    print FH $aScoreV->[0] . " " . $aScoreV->[1] . " " . $aScoreV->[2] . " " .
    $aScoreV->[3] . " " . $aScoreV->[4] . " " . $aScoreV->[5] . " 1\n";
  }
  foreach my $aScoreV (@rScoreV) {
    if (!defined($aScoreV->[5])) { print $aScoreV;}
    print FH $aScoreV->[0] . " " . $aScoreV->[1] . " " . $aScoreV->[2] . " " .
    $aScoreV->[3] . " " . $aScoreV->[4] . " " . $aScoreV->[5] . " 0\n";
  }
  close(FH);
}


sub normalizeDelta {
  my %A = %{$_[0]};
  my %B = %{$_[1]};
  my ($lower,$max_peak);
  foreach my $id (keys %A) {
    ($lower,$max_peak) = ($A{$id}->[0]->[2] > $B{$id}->[0]->[3] ?
                         ($B{$id},$A{$id}->[0]->[3]):
                         ($A{$id},$B{$id}->[0]->[3]));
    # recalculate DeltaXcorrs...
    foreach my $scoreV (@{$lower}) {
      $scoreV->[1] = ($max_peak-$scoreV->[3])/$max_peak;
    }
  }
}

sub falsePositives {
  my @true=@{$_[0]};
  my @false=@{$_[1]};
  my $fp=0;
  my @plot=();
  my $search = 1;
  foreach my $tp (0 .. $#true) {
    my $score = $true[$tp];
    while($false[$fp]>$score) { $fp++;}
    my $fdr = $fp/($fp+($tp+1));
    $search = 0 & print "FDR is 1% at $score\n" if ($search && $fdr>0.01);
    push @plot,$fdr;
  }
  return @plot;
}

sub printRes {
  my @true=@{$_[0]};
  my @false=@{$_[1]};
  my $fname = $_[2];
  my $tix=0;
  my $fix=0;
  open(FH,$fname);
  while(($tix<$#true) || ($fix<$#false)) {
#    print "$tix $#true $fix $#false \n";
    if(($tix>$#true)||($false[$fix]>$true[$tix])) {
      print FH "0\n";
      $fix++;
    } else {
      print FH "1\n";
      $tix++;
    }

  }
  close(FH)
}

open(*FH,"< $forwardFN");
my %fHits = readSQT(\*FH);
close(*FH);

open(*FH,"< $randomFN");
my %rHits = readSQT(\*FH);
close(*FH);

# normalizeDelta(\%fHits,\%rHits);

my (@bestFHits,@bestRHits);
my ($fscores,$rscores);

my @fH = values %fHits;
my @rH = values %rHits;
my (@fp,$fh,$rh);

($fh,$fscores) = allTheBest(\@fH);
($rh,$rscores) = allTheBest(\@rH);
@bestFHits = @$fh;
@bestRHits = @$rh;

foreach my $n (1 .. $NITER) {
  @bestFHits = @$fh;
  @bestRHits = @$rh;
  @fp = falsePositives($fscores,$rscores);
  my @fY = map {[1,$_->[2],$_->[3]]} @bestFHits;
  my @rY = map {[-1,-$_->[2],-$_->[3]]} @bestRHits;
#  my @fY = map {[1,$_->[0],$_->[1],$_->[2],$_->[3],$_->[4],$_->[5]]} @bestFHits;
#  my @rY = map {[-1, -$_->[0], -$_->[1], -$_->[2], -$_->[3], -$_->[4], -$_->[5]]} @bestRHits;

  push @fY,@rY;

  my @bb = map {[1]} @fY;

  my $yM = \@fY;
  my $bM = \@bb;
  bless $yM, ref($templ);
  bless $bM, ref($templ);

  my $aM=$yM->pinvert()->multiply($bM);

  $w_values[2] = $aM->[1]->[0];
  $w_values[3] = $aM->[2]->[0];
#   $w_values[0] = $aM->[1]->[0];
#   $w_values[1] = $aM->[2]->[0];
#   $w_values[2] = $aM->[3]->[0];
#   $w_values[3] = $aM->[4]->[0];
#   $w_values[4] = $aM->[5]->[0];
#   $w_values[5] = $aM->[6]->[0];

#  print $w_values[2] . " " . $w_values[3] . "\n";

  ($fh,$fscores) = allTheBest(\@fH);
  ($rh,$rscores) = allTheBest(\@rH);

#  if ($n==1) {
#    printScores(\@bestFHits,"> forward.out");
#    printScores(\@bestRHits,"> random.out");
#    printFP(\@fp,"> fplot.out");
#  }
}
# printScores(\@bestFHits,"> forwardItr.out");
# printScores(\@bestRHits,"> randomItr.out");
#printFP(\@fp,"> fplotItr.out");
printPyML2(\@bestFHits,\@bestRHits,"> 2d.pyml");
printPyML6(\@bestFHits,\@bestRHits,"> 6d.pyml");
#printPyML6(\@fH,\@rH,"> 6d.pyml");
printRes($fscores,$rscores, "> lda-2d.res");




