my $origin_seq=<STDIN>;
my $revcomp = reverse $origin_seq;
$revcomp =~ tr/ATGCatgc/TACGtacg/;

print "$revcomp\n";




