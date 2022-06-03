open I, "<$ARGV[0]"; # read with na and time.
$posf=$ARGV[0];
$outname=$ARGV[0];
$posf=~s/\.(\w+)?mod.tsv/.combined.RData.tsv/g;
$outname=~s/tsv$/pos.tsv/g;
open I2, "<$posf"; # each position and the kmer.
print "$outname\t$outname\n";
open O, ">$outname";

#read the each position and kmer.
while(<I2>){
chomp;
@data=split/\t/,$_;
$position{$data[2]}=$data[6];

}
close I2;
foreach(keys(%position))
{
	#print "$_\t$position{$_}\n";
}

while(<I>)
{
chomp;
@data=split /\t/,$_;
$readname=shift @data;
#print "$readname\n";
$n=0;
foreach(@data){
if($_=~/NA/){
}else{
print O "$n\t$_\t$position{$n}\n";
}

	$n++;
}

}

close O;

