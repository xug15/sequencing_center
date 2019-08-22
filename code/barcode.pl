print "Usage barcode.pl barcode.txt fq1.fq fq2.fq \nbarcode.txt: like below:\nbarcode1\tbarcode2\tfile name\nGAGCGCTA\tGATAGACA\tH0_4.7\nTAGCGAGT\tATGCCTAA\tH8_6.22\n\n======================\n";
open I1, "<$ARGV[0]"; # barcode file
#open I2, "<$ARGV[0]"; # fq1
#open I3, "<$ARGV[0]"; # fq2

my(%barcode,@barcode1,@barcode2,@name);
while(<I1>) #read barcode file.
{
    chomp;
    @data=split /\t/, $_;
    push @barcode1, $data[0];
    push @barcode2, $data[1];
    push @name, $data[2];
}


for my $i (0 .. $#barcode1)
{
print "$i\t$barcode1[$i]\t$barcode1[$i]\t$name[$i]\n";
}
my(@seq);
while(<I3>)
{
    chomp;
    $seq=<I3>;
    $qn =<I3>;
    $quality =<I3>;
    $seq_sub=substr($seq,0,20);
    push @seq,$seq_sub;
}




