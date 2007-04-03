#!/usr/bin/perl -w

# Reading file with Inspect output and reformat into a GIST file

use strict;
local *IN,*DATA,*LABEL;

sub isTrypN {
    my $seq = shift;
    return (($seq =~ m/^[A-Z]*[KR][.][^P]/g) || ($seq =~ m/^[*][.]/g))?1:0;
}

sub isTrypC {
    my $seq = shift;
    return (($seq =~ m/[KR][.][^P][A-Z]*$/g) || ($seq =~ m/[.][*]$/g))?1:0;
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
      split /\t/,$line;
      my $id =  $_[0] . '_' . $_[1] . '_' . $_[4];
      my $prot= (split / /,$_[3])[0];
      my $pep=$_[2];
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
      split /\t/,$line;
      my $id =  $_[0] . '_' . $_[1] . '_' . $_[4];
      my $prot= (split / /,$_[3])[0];
      my $pep=$_[2];
      if ($id ne $oldId) {
        my @of = ($id, $_[6], $_[7],$_[8],$_[5],$_[10], $_[11], $_[12],length($pep)-4 ,chargeVec($_[4]),
          isTrypN($pep), isTrypC($pep), missedTryp($pep), $protCnt{$prot}, $pepCnt{$pep}, scalar(@{$uniqPep{$prot}}));
        print DATA join("\t",@of) . "\n";
        print LABEL "$id\t$label\n";
        $oldId=$id;
      }
  }
}

#0            1     2          3       4      5       6        7         8         9       10      11      12           13           14        15
#SpectrumFile Scan# Annotation Protein Charge MQScore CutScore IntenseBY BYPresent Unused  p-value DeltaCN DeltaCNOther RecordNumber DBFilePos SpecFilePos

my $data = shift @ARGV;
my $label = shift @ARGV;
my $n = shift @ARGV;
my $s = shift @ARGV;
my $s2 = shift @ARGV;

print "$data $label $n $s $s2\n";
open(*DATA,"> $data");
print DATA "Id\tCutScore\tIntenseBY\tBYPresent\tMQScore\tp-value\tDeltaCn\tDeltaCnOther\tpepLen\tz1\tz2\tz3\ttrypN\ttrypC\tmissedTryp\tnumProt\tnumPep\tpepSite\n";
open(*LABEL,"> $label");
print LABEL "Id\tLabel\n";
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
close(*LABEL);
close(*DATA);


