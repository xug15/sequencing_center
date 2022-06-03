open I1, "</home/xugang/data/rnastructure_database/all_chr_metrics_forward.out_.gff3";
open I2, "</home/xugang/data/rnastructure_database/all_chr_metrics_reverse.out_.gff3";
open O, ">genome.mef.gc.tsv";
print O "chr\tstart\tend\tstrand\tmfe\tzscore\tpvalue\tdnsDiv\tFreqmef\tgc\n";
sub calcgc {
 my $seq = $_[0];
 my @seqarray = split('',$seq);
 my $count = 0;
 foreach my $base (@seqarray) {
   $count++ if $base =~ /[G|C|g|c]/i;
 }
 my $len = $#seqarray+1;
 my $num=$count / $len;
 #my ($dec)=$num =~ /(\S{6})/;
 return $num;
}




while(<I1>){
chomp;
# use split array
@data=split/\t/,$_;
$data[0]=~s/^chr//g;
$data[8]=~/MFE=(.*?);Z-Score=(.*?);P-value=(.*?);EnsDiv=(.*?);freqMFE=(.*?);Sequence=(.*?);/;
$mfe=$1;
$zscore=$2;
$pv=$3;
$dnsdiv=$4;
$freqmef=$5;
$seq=$6;
$gc=calcgc($seq);
#output
#print
print O "$data[0]\t$data[3]\t$data[4]\t$data[6]\t$mfe\t$zscore\t$pv\t$dnsdiv\t$freqmef\t$gc\n";

}


while(<I2>){
chomp;
# use split array
@data=split/\t/,$_;
$data[0]=~s/^chr//g;
$data[8]=~/MFE=(.*?);Z-Score=(.*?);P-value=(.*?);EnsDiv=(.*?);freqMFE=(.*?);Sequence=(.*?);/;
$mfe=$1;
$zscore=$2;
$pv=$3;
$dnsdiv=$4;
$freqmef=$5;
$seq=$6;
$gc=calcgc($seq);
#output
#print
print O "$data[0]\t$data[3]\t$data[4]\t$data[6]\t$mfe\t$zscore\t$pv\t$dnsdiv\t$freqmef\t$gc\n";

}



close I1;
close I2;

