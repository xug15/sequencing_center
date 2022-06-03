open I, "<$ARGV[0]";
open O, ">$ARGV[0].neam.pos.tsv";
<I>;
while(<I>){
chomp;
@data=split/\t/,$_;
#print "$_\n";
$time{$data[0]."\t".$data[1]}+=$data[4];
$current{$data[0]."\t".$data[1]}+=$data[2];
$count{$data[0]."\t".$data[1]}++;
}
#print O "name\tpos\ttime\tcurrent\tcount\n";
foreach(keys(%time)){
$timem=$time{$_}/$count{$_};
$currentm=$current{$_}/$count{$_};
	print O "$_\t$timem\t$currentm\t$count{$_}\n";

}
close I;
close O;
system("sort -k1,1 -k2,2n ".$ARGV[0].".neam.pos.tsv > ".$ARGV[0].".neam.pos2.tsv");
#system('echo -e "chrosome\tpos\ttime\tcurrent\tcount\n"> tmp.head ');
open I, "<$ARGV[0].neam.pos2.tsv";
open O, ">$ARGV[0].neam.pos.tsv";
print O "name\tpos\ttime\tcurrent\tcount\n";
while(<I>){
print O "$_";
}
close O;
close I;
system("rm  ".$ARGV[0].".neam.pos2.tsv")

