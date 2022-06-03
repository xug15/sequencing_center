open I,  "<$ARGV[0]";#ks file
#read ks file
$limit_value=$ARGV[2];
$na_value=$limit_value;

<I>;

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
#print "$pv2pro\t$pro\n";
$pro{$data[1]}=$pro;
}
close I;

for ($i=0;$i<30;$i++){

	#print "$i\n";

open I2, "<$ARGV[1]";#fasta file
#read file
my $seqp=0;
my $shift=$i;
open O, ">$ARGV[0].$shift.mv.tsv";

print O "RNA\tIndex\tProbability\tBase\n";
while(<I2>){
chomp;
if($_=~/>(.*)/)
{
$name=$1;
$name=~s/ /_/g;
}else{
@seq=split(undef,$_);
foreach(@seq){
$_=~s/T/U/g;
$posp=$seqp+$shift;
if(exists($pro{$posp})){
$value=$pro{$posp};

#$value=(1-$value);
}else{
$value="NA";
}
$seqp++;
#print "$name\t$seqp\t$value\t$_\n";
 if($value=~/NA/){
 $nv="NA";
 }
 else{
 $nv=$value;
 }
 print O "$name\t$seqp\t$nv\t$_\n";
}
}  
}
close I2;
close O;

}



