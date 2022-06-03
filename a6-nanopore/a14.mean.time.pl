open I, "<$ARGV[0]";
open O, ">$ARGV[0].neam.tsv";
<I>;
while(<I>){
chomp;
@data=split/\t/,$_;
if($#data==4){
#print "$_\n";
$time{$data[0]}+=$data[4];
$current{$data[0]}+=$data[2];
$count{$data[0]}++;
}
}
print O "name\ttime\tcurrent\n";
foreach(keys(%time)){
$timem=$time{$_}/$count{$_};
$currentm=$current{$_}/$count{$_};
	print O "$_\t$timem\t$currentm\n";

}



