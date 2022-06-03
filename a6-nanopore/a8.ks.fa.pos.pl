open I,  "<$ARGV[0]";#ks file
open O, ">$ARGV[0].transform.tsv";
#read ks file

$head=<I>;
chomp($head);
print O "$head\ttran\n";
while(<I>){
chomp;

#print "$_\n";
@data=split /\t/,$_;
$dfpvsm{$data[1]}=$data[6];
$dfpv{$data[1]}=$data[5];
$lin{$data[1]}=$data[7];

$df{$data[1]}=$data[2];
$pv{$data[1]}=$data[3];

if($data[3]=~/NA/){
$pv2pro="NA";

$pro='NA';
}else{
$pv2pro=-log($data[3]+0.0000000000000001);
$pv2pro=(2**$pv2pro)/(1+2**$pv2pro);
$pro=$pv2pro*$data[2];
}
$pro{$data[1]}=$pro;
print O "$_\t$pro\n";
}
close I;




