open I, "<$ARGV[0]";
open O, ">$ARGV[0].txt";
open O3, ">$ARGV[0].sta.tsv";
print O "chr\tposition\tcurrent\tsd\ttime\tref\tkmer\n";
while(<I>){
chomp;
@data=split /\t/,$_;
#print "$data[0]\t$data[1]\t$data[6]\t$data[7]\t$data[8]\n";
if($data[1]=~/\d/){
$hash_current{$data[0]."\t".$data[2]}+=$data[4];
$hash_st{$data[0]."\t".$data[2]}+=$data[3];
$hash_time{$data[0]."\t".$data[2]}+=$data[5];
$hash_num{$data[0]."\t".$data[2]}++;
$kmer{$data[0]."\t".$data[2]}=$data[6];
#print  "$data[0]\t$data[1]\t$data[6]\t$data[7]\t$data[8]\n";
$data[6]=~/^((\w)\w+)/;
$ref=$2;
$kmer=$1;
print O "$data[0]\t$data[2]\t$data[4]\t$data[3]\t$data[5]\t$ref\t$kmer\n";
}
}
print O3 "chr\tposition\tcurrent\tsd\ttime\tkmer\n";
foreach(keys(%hash_num)){
	$current=$hash_current{$_}/$hash_num{$_};
	$st=$hash_st{$_}/$hash_num{$_};
	$time=$hash_time{$_}/$hash_num{$_};
print O3 "$_\t$current\t$st\t$time\t$kmer{$_}\n";

}
close O3;
close O1;

system("sort -k1,1 -k2,2n $ARGV[0].sta.tsv > $ARGV[0].sta.tmp && mv $ARGV[0].sta.tmp $ARGV[0].sta.tsv")

