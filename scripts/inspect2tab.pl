#!/usr/bin/perl -w

# Reading file with Inspect output and reformat into a GIST file

use strict;
local (*IN,*DATA,*LABEL);

sub isTrypN {
    my $seq = shift;
    return (($seq =~ m/^[KR][.][^P]/g) || ($seq =~ m/^[*][.]/g))?1:0;
}

sub isTrypC {
    my $seq = shift;
    return (($seq =~ m/[KR][.][^P]$/g) || ($seq =~ m/[.][*]$/g))?1:0;
}

sub missedTryp {
    my $seq = shift;
    my $count = 0;
    for (; $seq =~ /[KR][^P.]/g; $count++) { }
    return $count;
}


sub chargeVec {
    my $charge= shift;
    return ($charge==1?1:0,$charge==2?1:0,$charge==3?1:0);
}

sub convert {
  my $label = shift;
  my $line = <IN>; # skipp annotation line
  my $oldId = "xx";
  my (%protCnt,%pepCnt,%uniqPep);
  while ($line=<IN>) {
      chomp $line;
      my @f = split(/\t/,$line);
      my $prot = (split / /,$f[3])[0];
      my $pep = $f[2];
      my $id =  $pep . '-' . $prot . '-' . $f[1] . '-' . $f[4];
      if ($id ne $oldId) {
        $protCnt{$prot}++;
        $pepCnt{$pep}++;
        push @{$uniqPep{$prot}},$pep unless grep {$_ eq $pep} @{$uniqPep{$prot}};
      }
      $oldId=$id;
  }
  seek(IN,0,0);
  $line = <IN>; # skipp annotation line
  $oldId = "xx";
  while ($line=<IN>) {
      chomp $line;
      next unless length($line);
      my @f = split(/\t/,$line);
      my @prots = split / /,$f[3];
      my $prot= $prots[0];
      my $pep=$f[2];
      my $id =  $pep . '-' . $prot . '-' . $f[1] . '-' . $f[4];
      if ($id ne $oldId) {
        my @of = ($id, $label, $f[5], $f[7],$f[8],$f[9],$f[10], $f[11], $f[14],
           $f[15], $f[16], length($pep)-4 ,chargeVec($f[4]),
          isTrypN($pep), isTrypC($pep), missedTryp($pep)
#          isTrypN($pep), isTrypC($pep), missedTryp($pep), $protCnt{$prot}, $pepCnt{$pep}, scalar(@{$uniqPep{$prot}}));
          ,$pep, join("\t",@prots));
        print join("\t",@of) . "\n";
        $oldId=$id;
      }
  }
}
#0            1     2          3       4      5       6        7             8              9         10
#NewSpectrumF Scan# Annotation Protein Charge MQScore Length   TotalPRMScore MedianPRMScore FractionY FractionB

#11        12  13      14      15         16              17              18        19
#Intensity NTT p-value F-Score DeltaScore DeltaScoreOther NewRecordNumber DBFilePos SpecFilePos


my $n = shift @ARGV;
my $s = shift @ARGV;

my $s2 = shift @ARGV;

#print "$n $s $s2\n";
print "Id\tlabel\tMQScore\tTotalPRMScore\tMedianPRMScore\tFractionY\tFractionB\tIntensity\tF-Score\tDeltaScore\tDeltaScoreOther\tpepLen\tz1\tz2\tz3\ttrypN\ttrypC\tmissedTryp\tnumProt\tnumPep\tpepSite\tPeptide\tProtein(s)\n";
open(*IN,"< $n");
convert("1");
close(*IN);
open(*IN,"< $s");
convert("-1");
close(*IN);
if ($s2) {
  open(*IN,"< $s2");
  convert("-2");
  close(*IN);
}


