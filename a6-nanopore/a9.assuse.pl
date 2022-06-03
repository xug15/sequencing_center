open I, "<$ARGV[0]"; # ct 1
open I2, "<$ARGV[1]"; # ct2 

open N1, "<$ARGV[2]"; #a ct 1 begion
open N2, "<$ARGV[3]"; #b ct 1 end
open N3, "<$ARGV[4]"; #c ct 2 begion
open N4, "<$ARGV[5]"; #d ct 2 end.
# 
$a=10;
$b=50;
$c=20;
$d=60;
$a=$ARGV[2];
$b=$ARGV[3];
$c=$ARGV[4];
$d=$ARGV[5];

#print "$a\t$b\t$c\t$d\n";
#read ct 1 file rna name
# read ct 2 file rna name.
$ct1_name=<I>;
$ct2_name=<I2>;
#set count n1 and n2.
$n1=0;
$n2=0;

while(<I>){
chomp;
$n1++;
if($n1 >= $a && $n1 <=$b){
#print "$_\n";
$_=~s/\s+/\t/g;
@data=split /\t/,$_;
shift @data;
$dis1=$a-1;


$data[0]=$data[0]-$dis1;
$data[2]=$data[2]-$dis1;
$data[3]=$data[3]-$dis1;
$data[5]=$data[5]-$dis1;

 if($data[4] == 0){
	 #$data[4]=0;
 }else{
	 $data[4]=$data[4]-$dis1;
 }
$info=join "\t", @data;
#print "$info\n";
$ct1{$data[0]}=$data[4];
$ct1q{$data[0]}=$data[1];
}

}
close I;

# read the ct2 file
while(<I2>){
chomp;
$n2++;
if($n2 >= $c && $n2 <=$d){
#print "$_\n";
$_=~s/\s+/\t/g;
@data=split /\t/,$_;
shift @data;
$dis1=$c-1;


$data[0]=$data[0]-$dis1;
$data[2]=$data[2]-$dis1;
$data[3]=$data[3]-$dis1;
$data[5]=$data[5]-$dis1;

 if($data[4] == 0){
         #$data[4]=0;
 }else{
         $data[4]=$data[4]-$dis1;
 }
$info=join "\t", @data;
#print "$info\n";
$ct2{$data[0]}=$data[4];
$ct2q{$data[0]}=$data[1];
}

}
close I2;

$total=0;
$pair=0;
$unpair=0;


foreach(keys(%ct1)){
	#print "$_\t$ct1{$_}\t$ct1q{$_}\n";
$total++;
if(exists($ct2{$_}))
 {
 if($ct1{$_} > 0 &&  $ct2{$_} > 0 )
  {
  $pair++;
  }
  if($ct1{$_} < 0 &&  $ct2{$_} > 0 )
  {
  $pair++;
  }if($ct1{$_} > 0 &&  $ct2{$_} < 0 )
  {
  $pair++;
  }
  if($ct1{$_} < 0 &&  $ct2{$_} < 0 )
  {
  $pair++;
  }
 if($ct1{$_} == 0 &&  $ct2{$_} == 0 ){
   $unpair++;
  }

 }
}
$acc=($pair+$unpair)/$total;
print "$ARGV[1]\t$total\t$pair\t$unpair\t$acc\n";

