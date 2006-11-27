#!/usr/bin/env perl 
# Reading file with Inspect output and reformat into a GIST file

#0            1     2          3       4      5       6        7         8         9       10      11      12           13           14        15
#SpectrumFile Scan# Annotation Protein Charge MQScore CutScore IntenseBY BYPresent Unused  p-value DeltaCN DeltaCNOther RecordNumber DBFilePos SpecFilePos

sub isTrypN {
    my $seq = shift;
    return (($seq =~ m/^[A-Z]*[KR][.][^P]/g) || ($seq =~ m/^[*][.]/g))?1:0;
}
 
sub isTrypC {
    my $seq = shift;
    return (($seq =~ m/[KR][.][^P][A-Z]*$/g) || ($seq =~ m/[.][*]$/g))?1:0;
}

sub chargeVec {
    my $charge= shift;
    return ($charge==1?1:0,$charge==2?1:0,$charge==3?1:0);
}

open(*FH,"< " . shift ARGV);
my $line = <FH>; # skipp annotation line
print "Id\tCutScore\tp-value\tDeltaCn\tMQScore\tDeltaCnOther\tpepLen\tz1\tz2\tz3\ttrypN\ttrypC\n";
my $oldId = "xx";
while ($line=<FH>) {
    chomp $line;
    split /\t/,$line;
    my $id =  $_[0] . '_' . $_[1] . '_' . $_[4];
    my @of = ($id, $_[6], $_[10], $_[11], $_[5], $_[12],length($_[2])-4 ,chargeVec($_[4]),isTrypN($_[2]), isTrypC($_[2]));
    print join("\t",@of) . "\n" if ($id ne $oldId);
    $oldId=$id;
}
close(*FH);
