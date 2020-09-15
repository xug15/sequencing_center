open DATA, "<$ARGV[0]";
$name=$ARGV[0];
$name=~s/.txt//g;
$name="c-out".$name;
system("mkdir $name");
open OUT1, ">$name/$ARGV[0].transcript_up.txt";
open OUT2, ">$name/$ARGV[0].transcript_do.txt";
open OUT3, ">$name/$ARGV[0].translatio_up.txt";
open OUT4, ">$name/$ARGV[0].translatio_do.txt";
$head=<DATA>;
print OUT1 "$head";
print OUT2 "$head";
print OUT3 "$head";
print OUT4 "$head";
while(<DATA>)
{
chomp;
@data=split /\t/,$_;
if($data[13]=~/transcription_only/)
{
if($data[1]>0)
{
   print OUT1 "$_\n";
}else{
   print OUT2 "$_\n";
}
}
if($data[13]=~/translation_only/)
{
if($data[2]>0)
{
   print OUT3 "$_\n";
}else{
   print OUT4 "$_\n";
}
}

}




