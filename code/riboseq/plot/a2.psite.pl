use List::Util qw(min max);

%dis=(27=>12, 
28=>12, 
29=>12, 
30=>12, 
31=>13, 
32=>13, 
33=>13, 
34=>13, 
35=>14);
print $dis{28};

open I, "<$ARGV[0]";
open O, ">$ARGV[0].psite.tsv";
open O2, ">$ARGV[0].clean.tsv";

while(<I>){
chomp;
@data=split /\t/,$_;
$data[5]=~/(\d+)M/g;
$len=$1;
$pos=$data[3];
$dis=$len;
if($data[4]>0 || $data[4]<0){
if(exists($dis{$len})){
$dis=$dis{$len};
}elsif($len>34){
$dis=35;
}
}
$psite=$dis+$pos;
$hash{$psite}++;
}




print O "position\tnumber\n";
$maxb=max(keys(%hash));

for (my $i=1;$i<$maxb+1;$i++){
if(exists($hash{$i})){
print O "$i\t$hash{$i}\n";
print O2 "$i\t$hash{$i}\n";
}else{
print O "$i\t0\n";
}

}

close O;
close I;
system("sort -k1,1n ".$ARGV[0].".psite.tsv > ".$ARGV[0].".sort.tsv");
system("rm ".$ARGV[0].".psite.tsv");



