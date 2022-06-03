open I, "<$ARGV[0]";
open I2, "<genome.mef.tsv";
open O, ">$ARGV[0].merge.tsv";
while(<I2>){
chomp;
@data=split/\t/,$_;
$name=$data[0]."\t".$data[3]."\t".$data[1];
$hash{$name}=$data[4]."\t".$data[5]."\t".$data[6]."\t".$data[7]."\t".$data[8];
#print "$name\n";
}
close I2;

while(<I>){
chomp;
@data=split /\t/,$_;
$reminder=$data[6]%40;
$transform=$data[6]-$reminder+1;
#print "$transform\n";
#print "$_\t$transform\n";
$name=$data[5]."\t".$data[7]."\t".$transform;
#print "2:$name\n";
if(exists($hash{$name})){
print O "$_\t$transform\t$hash{$name}\n";
}
}

close I;
close O;

