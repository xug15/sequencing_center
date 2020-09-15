open DATA, "<$ARGV[0]";
$name=$ARGV[0];
$name=~s/.txt//g;
$name="c-out".$name;
system("mkdir $name");


while(<DATA>)
{
   chomp;
   @data=split /\t/,$_;
   if($data[13]=~/homodirectional/){}
   elsif($data[13]=~/opposite_change/){}
   elsif($data[13]=~/^$/){}
   elsif($data[13]=~/stable/){}else
   {
      $hash{$data[13]}='';
   }
}
close DATA;
@name=keys(%hash);
print "@name\n";
open DATA, "<$ARGV[0]";
open OUT1, ">$name/$name[0].up.txt";
open OUT2, ">$name/$name[0].do.txt";
open OUT3, ">$name/$name[1].up.txt";
open OUT4, ">$name/$name[1].do.txt";
$head=<DATA>;
print OUT1 "$head";
print OUT2 "$head";
print OUT3 "$head";
print OUT4 "$head";
while(<DATA>)
{
   chomp;
   @data=split /\t/,$_;
   if($data[13]=~/$name[0]/)
   {
      if(abs($data[5])>1)
      {
         if($data[5]>0){
            print OUT1 "$_\n";
         }else{
            print OUT2 "$_\n";
         }
      }
      if(abs($data[6])>1)
      {
         if($data[6]>0)
         {
            print OUT1 "$_\n";
         }else{
            print OUT2 "$_\n";
         }
      }
   }
   if($data[13]=~/$name[1]/)
   {
            if(abs($data[5])>1)
      {
         if($data[5]>0){
            print OUT3 "$_\n";
         }else{
            print OUT4 "$_\n";
         }
      }
      if(abs($data[6])>1)
      {
         if($data[6]>0)
         {
            print OUT3 "$_\n";
         }else{
            print OUT4 "$_\n";
         }
      }
   }  

}


