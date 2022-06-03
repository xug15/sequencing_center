open I, "<$ARGV[0]"; # iput dot parent file.

while(<I>){
chomp;
if($_=~/^>/)
 {
 $name=$_;
 $file=$name;
 $file=~s/>//g;
   print "$name\n$file.ct\n";
   open O, ">$file.ct";
   print O "$name\n";
 $seq=<I>;
 $str=<I>;
 chomp($seq);
 chomp($str);
 #print "seq:\n$seq\n";
 #print "structure:\n$str\n";
 @seq=split(undef,$seq);
 @str=split(undef,$str);
  for($a=0;$a<$#seq;$a++){
   $b=$a+1;
   $c=$a+2;
   if($str[$a]=~/\./){
   $pair=0;
   }elsif($str[$a]=~/\(/){
   $pair=1;
   }elsif($str[$a]=~/\)/){
   $pair=2;
   }else{
   $pair=3;
   }
   print O "  $b  $seq[$a]  $a  $c  $pair  $b\n";
  }
  close O;
 }

}





