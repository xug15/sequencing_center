
open I, "<$ARGV[0]";

open O, ">$ARGV[0].tab.tsv";
open O2, ">$ARGV[0].fa.tsv";
$length=1000;
$length=$ARGV[1];
<I>;
while(<I>){
@data=split /\t/,$_;
push @length,$data[2];
}

close I;
@array = sort { $a <=> $b } @length;
#print "$array[0] $array[-1]\n";

$length=$array[-1];

open I, "<$ARGV[0]";
<I>;
#time tsv
while(<I>){
chomp;
#split infor by table
@data=split /\t/,$_;
#set default array length.
${$data[1]}[$length+1]='NA';

#use read name as array name,positon is index, time is the value.
${$data[1]}[$data[2]+1]=$data[5];
#get the read sequence character in array.
@seq=split(undef,$data[6]);

${$data[0]}[$data[2]]=$seq[0];

# use hash to recoder the array of time name.
$array{$data[1]}='';

# use hash to recoder the array of sequence character.
$array_seq{$data[0]}='';

#set read name as 
${$data[1]}[0]=$data[1];

}

close I;

foreach(keys(%array)){
	#print "$_\n";
$arrayname=$_;
foreach(@{$arrayname}){
 if($_=~/^$/){
  $_=NA;
  }
 }
}



foreach(keys(%array_seq)){
        #print "$_\n";
$arrayname=$_;
foreach(@{$arrayname}){
 if($_=~/^$/){
  $_=N;
  }
 }
}



foreach(keys(%array)){

$arrayname=$_;
        #print "@{$_}\n";
$info=join "\t", @{$arrayname};

=pod
n=0;
foreach( @{$arrayname}){
print "$n\t$_\n";
$n++;
}
=cut

print O "$info\n";

}

foreach(keys(%array_seq)){

$arrayname=$_;	
	#print "@{$_}\n";
$info=join "\t", @{$arrayname};


print O2 "$arrayname\t$info\n";

}
close O;
close O2;


