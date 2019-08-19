#!/usr/bin/perl
($file1,$file2,$bar,$ori,$end)=@ARGV;

if((not defined $file1)||(not defined $file2)||(not defined $bar)||(not defined $ori)||(not defined $end)) {
	die "Usage= ./grep_.sh  $file1,$file2,$bar,$ori,$end"
}
$flag=0;
my $len=$end-$ori+1;
@code=split(//,$bar);
$b1=$code[0].".".$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b2=$code[0].$code[1].".".$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b3=$code[0].$code[1].$code[2].".".$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b4=$code[0].$code[1].$code[2].$code[3].".".$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b5=$code[0].$code[1].$code[2].$code[3].$code[4].".".$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b6=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].".".$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b7=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].".".$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b8=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].".".$code[9].$code[10].$code[11].$code[12].$code[13].$code[14].;
$b9=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].".".$code[10].$code[11].$code[12].$code[13].$code[14].;
$b10=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].".".$code[11].$code[12].$code[13].$code[14].;
$b11=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].".".$code[12].$code[13].$code[14].;
$b12=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].".".$code[13].$code[14].;
$b13=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].".".$code[14].;
$b14=$code[0].$code[1].$code[2].$code[3].$code[4].$code[5].$code[6].$code[7].$code[8].$code[9].$code[10].$code[11].$code[12].$code[13].".".;


open (FILE1,$file1) or die $!;
open (FILE2,$file2) or die $!;
$file3="m_".$file1;
$file4="m_".$file2;
open (FILE3,">$file3") or die $!;
open (FILE4,">$file4") or die $!;
$file5="mis_".$file1;
$file6="mis_".$file2;
open (FILE5,">$file5") or die $!;
open (FILE6,">$file6") or die $!;
while (my $line_1=<FILE1> ) {
        my $line_2=<FILE2> ;
	chomp $line_1;
       chomp $line_2;
      if(($flag%4)==0) {
           $L1_0=$line_1;
           $L2_0=$line_2;
           $flag++;
      }elsif(($flag%4)==1) {
           $L1_1=$line_1;
           $L2_1=$line_2;
           $flag++;
      }elsif(($flag%4)==2) {
           $L1_2=$line_1;
           $L2_2=$line_2;
           $flag++;
      }elsif(($flag%4)==3) {
           $L1_3=$line_1;
           $L2_3=$line_2;
           $flag++;
           my $st=substr($L1_1,$ori,$len);
           if ($st=~m/$b1|$b2|$b3|$b4|$b5|$b6/) {
		  print FILE3 "$L1_0\n" ;
		  print FILE3 "$L1_1\n" ;
               print FILE3 "$L1_2\n" ;
               print FILE3 "$L1_3\n" ;
            
               print FILE4 "$L2_0\n" ;
		  print FILE4 "$L2_1\n" ;
               print FILE4 "$L2_2\n" ;
               print FILE4 "$L2_3\n" ;
	    }else{
               print FILE5 "$L1_0\n" ;
		  print FILE5 "$L1_1\n" ;
               print FILE5 "$L1_2\n" ;
               print FILE5 "$L1_3\n" ;
            
               print FILE6 "$L2_0\n" ;
		  print FILE6 "$L2_1\n" ;
               print FILE6 "$L2_2\n" ;
               print FILE6 "$L2_3\n" ;

          }

      }
}


close FILE1 or die $!;
close FILE2 or die $!;
close FILE3 or die $!;
close FILE4 or die $!;
