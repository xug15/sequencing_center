open I, "<$ARGV[0]";
open O, ">$ARGV[0].neam.95.pos.tsv";
<I>;
while(<I>){
chomp;
@data=split/\t/,$_;
#print "$_\n";
$time{$data[0]."\t".$data[1]}.=$data[4]."\t";
$current{$data[0]."\t".$data[1]}.=$data[2]."\t";
$count{$data[0]."\t".$data[1]}++;
}
#print O "name\tpos\ttime\tcurrent\tcount\n";
foreach(keys(%time)){
	#print "$time{$_}\n";
@time=split /\t/,$time{$_};
@timesort=sort { $a <=> $b } @time;
#print "@timesort\n";
$sei=int(0.95*$#timesort);
#print "$sei\n";
$timem=$timesort[$sei];
@current=split /\t/,$current{$_};
@currentsort=sort {$a <=> $b } @current;
$sei=int(0.95*$#currentsort);
$currentm=$currentsort[$sei];
	print O "$_\t$timem\t$currentm\t$count{$_}\n";

}
close I;
close O;
system("sort -k1,1 -k2,2n ".$ARGV[0].".neam.95.pos.tsv > ".$ARGV[0].".neam.95.pos2.tsv");
system('echo -e "chrosome\tpos\ttime\tcurrent\tcount"> tmp.head ');
system("cat tmp.head ".$ARGV[0].".neam.95.pos2.tsv > " .$ARGV[0].".neam.95.pos.tsv");
system("rm ".$ARGV[0].".neam.95.pos2.tsv")
#open I, "<$ARGV[0].neam.pos2.tsv";
#open O, ">$ARGV[0].neam.pos.tsv";
#print O "name\tpos\ttime\tcurrent\tcount\n";
#while(<I>){
#print O "$_";
#}
#close O;
#close I;
#system("rm  ".$ARGV[0].".neam.pos2.tsv")


