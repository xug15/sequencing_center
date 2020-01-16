open DATA, '<lbxfc_results.txt';
open OUT1, ">transcript_only_up.txt";
open OUT2, ">transcript_only_down.txt";
open OUT3, ">translation_only_up.txt";
open OUT4, ">translation_only_down.txt";
open OUT5, ">transcript_translation_up.txt";
open OUT6, ">transcript_translation_down.txt";
open OUT7, ">transcript_up_translation_down.txt";
open OUT8, ">transcript_down_translation_up.txt";

open AOUT1, ">transcript_only_up.list";
open AOUT2, ">transcript_only_down.list";
open AOUT3, ">translation_only_up.list";
open AOUT4, ">translation_only_down.list";
open AOUT5, ">transcript_translation_up.list";
open AOUT6, ">transcript_translation_down.list";
open AOUT7, ">transcript_up_translation_down.list";
open AOUT8, ">transcript_down_translation_up.list";

#homodirectional
#opposite_change
#stable
#transcription_only
#translation_only
$head=<DATA>;
$head="GeneId\t".$head;
print OUT1 "$head";
print OUT2 "$head";
print OUT3 "$head";
print OUT4 "$head";
print OUT5 "$head";
print OUT6 "$head";
print OUT7 "$head";
print OUT8 "$head";
while(<DATA>)
{
chomp;
@data=split /\t/,$_;
if($data[13]=~/transcription_only/)
    {
        if($data[1]>0)
        {
    print OUT1 "$_\n";
    print AOUT1 "$data[0]\n";
        }else{
    print OUT2 "$_\n";
    print AOUT2 "$data[0]\n";
        }
        
    }
if($data[13]=~/translation_only/)
    {
        if($data[2]>0)
        {
    print OUT3 "$_\n";
    print AOUT3 "$data[0]\n";
        }else{
    print OUT4 "$_\n";
    print AOUT4 "$data[0]\n";
        }
        
    }
if($data[13]=~/homodirectional/)
    {
        if($data[1]>0)
        {
    print OUT5 "$_\n";
    print AOUT5 "$data[0]\n";
        }else{
    print OUT6 "$_\n";
    print AOUT6 "$data[0]\n";
        }
        
    }
if($data[13]=~/opposite_change/)
    {
        if($data[1]>0)
        {
    print OUT7 "$_\n";
    print AOUT7 "$data[0]\n";
        }else{
    print OUT8 "$_\n";
    print AOUT8 "$data[0]\n";
        }
        
    }

}


