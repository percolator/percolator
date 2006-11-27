#!/usr/bin/env perl 
# Reading file with Inspect output and reformat into a GIST file

#SpectrumFile   Scan#   Annotation      Protein Charge  MQScore CutScore        IntenseBY       BYPresent       Unused  p-value DeltaCN DeltaCNOther    RecordNumber    DBFilePos       SpecFilePos

sub isTrypN {
    my $seq = shift;
    return (($seq =~ m/^[A-Z]*[KR][.][^P]/g) || ($seq =~ m/^[*][.]/g));
}
 
sub isTrypC {
    my $seq = shift;
    return (($seq =~ m/[KR][.][^P][A-Z]*$/g) || ($seq =~ m/[.][*]$/g));
}

sub chargeVec {
    my $charge= shift;
    return ($charge==1,$charge==2,$charge==3);
}

open(*FH,"< " . shift ARGV);
my $line = <FH>; # skipp annotation line
print "Id\tCutScore\tp-value\tDeltaCN\tMQScore\tDeltaCnOther\tpepLen\tz1\tz2\tz3\ttrypN\ttrypC\n";
while ($line=<FH>) {
    chomp $line;
    split /\t/,$line;
    my @of = ($_[0] . '_' . $_[1] . '_' . $_[4], $_[6], $_[10], $_[11], $_[5], $_[12], chargeVec($_[4]),isTrypN($_[2]), isTrypC($_[2]));
    print join('\t',@of) . "\n";
}
close(*FH);
