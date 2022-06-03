open I, "<$ARGV[0]";
open O, ">$ARGV[0].ch";
$head=<I>;
if($head=~/^contig/){
print O $head;
}else{
print O "contig\tposition\treference_kmer\tread_name\tstrand\tevent_index\tevent_level_mean\tevent_stdv\tevent_length\tmodel_kmer\tmodel_mean\tmodel_stdv\tstandardized_level\n";
}
while(<I>){

print O "$_";
}
close I;
close O;

