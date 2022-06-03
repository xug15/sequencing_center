open I, "</home/xugang/data/reference/hg38/Homo_sapiens.GRCh38.100.chr.gtf.exon";
open I2, "<$ARGV[0]";
open O, ">$ARGV[0].chrome.tsv";
$old='oldna';
while(<I>){
chomp;
@data=split/\t/,$_;
$_=~/transcript_id "(.*?)";/;
$trans=$1;
$_=~/exon_number "(.*?)"/;
$exonnumer=$1;
$length=$data[4]-$data[3]+1;
if($old=~/^$trans$/){
$totallength+=$length;
}else{
$old=$trans;
$totallength=$length;
}

#print "$data[0]\t$data[3]\t$data[4]\t$length\t$totallength\t$data[6]\t$exonnumer\t$trans\n";
#$info{$trans."\t".$totallength}="$data[0]\t$data[3]\t$data[4]\t$length\t$totallength\t$data[6]\t$exonnumer\t$trans\n";

$name{$trans}.=$totallength."\t";

$chrom{$trans}=$data[0];
$start{$trans}.=$data[3]."\t";
$end{$trans}.=$data[4]."\t";
$length{$trans}.=$length."\t";
$totallength{$trans}.=$totallength."\t";
$strand{$trans}=$data[6];
$exonnum{$trans}.=$exonnumer."\t";


}
close I;

foreach(keys(%name)){
#print "$_\t$name{$_}\n";
@name=split/\t/,$name{$_};
#print "$#name\t$name[0]\t$name[1]\n";

}

$name='ENST00000666741';
$name='ENST00000341290';
$p1=100;
$p2=200;
$p3=700;
$p=$p3;
#print "$name\t$p\n";

while(<I2>){
chomp;
@d2=split /\t/,$_;
$name=$d2[0];
$p=$d2[1]+1;
$infoa=$_;
if(exists($totallength{$name})){
$totallength=$totallength{$name};
#print "$totallength\n";
@totallength=split/\t/,$totallength; # get total length.
$chrom=$chrom{$name};# get chromsome.
@start=split/\t/,$start{$name};
@end=split/\t/,$end{$name};
@length=split/\t/,$length{$name};
$strand=$strand{$name};
@exonnum=split /\t/,$exonnum{$name};
#set indexes.
$n=0;
#set distance minier 0 and next length...
$distance_minier=0;
if($strand=~/\+/){
foreach(@totallength){
#print "$_\n";
$totall=$_;
# set the last position, and the next will subtrace.
#print "$chrom\t$start[$n]\t$end[$n]\t$length[$n]\t$strand\t$exonnum[$n]\t$totall\n";
if($p<=$totall){ # ditermine the position is smaller than totall length. Then it will be use for next 
$distance=$p-$distance_minier-1;
$start_f=$start[$n]+$distance;
#print "p:$p\t distance_miner: $distance_minier\t distance: $distance\n";
print O "$infoa\t$chrom\t$start_f\t$strand\n";

last;
}else{
$distance_minier=$totall;
}
$n++; ## $n as index, and it will be useful to get strand, chrome, start, end,
}


}
if($strand=~/\-/){
foreach(@totallength){
#print "$_\n";
$totall=$_;
# set the last position, and the next will subtrace.
#print "$chrom\t$start[$n]\t$end[$n]\t$length[$n]\t$strand\t$exonnum[$n]\t$totall\n";
if($p<=$totall){ # ditermine the position is smaller than totall length. Then it will be use for next
$distance=$p-$distance_minier-1;
$start_f=$end[$n]-$distance;
#print "p:$p\t distance_miner: $distance_minier\t distance: $distance\n";
print O "$infoa\t$chrom\t$start_f\t$strand\n";

last;
}else{
$distance_minier=$totall;
}
$n++; ## $n as index, and it will be useful to get strand, chrome, start, end,
}
}
}


}
close O;


