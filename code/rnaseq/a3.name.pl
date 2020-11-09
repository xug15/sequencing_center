open D1, "<b4-gfold/YXXKO2-1.read_cnt";


while(<D1>)
{
chomp;
@data=split /\t/,$_;

$hash{$data[0]}=$data[1];
#print "$data[0]\t$data[1]\n";
$data[0]=~/(\w+)\./g;
$hashg{$1}=$data[1];
#print "$1\t$data[1]\n";
}
close D1;
open D2, "<$ARGV[0]";
open O2, ">$ARGV[0].g.tsv";
$head=<D2>;
print O2 "gene\t$head";
while(<D2>)
{
    chomp;
@data=split/\t/,$_;
if(exists($hash{$data[0]})){
    $data[0]=$hash{$data[0]};
    $info=join ("\t",@data);
    print O2 "$info\n";
}

}








