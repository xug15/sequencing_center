open DATA,  "<mouse-tRNAs-confidence-set.txt";
open OUT, ">mouse-tRNAs-confidence.txt";
<DATA>;

print OUT "Chr\ttRNA\tBegin\tEnd\tIsotype\tAnticodon\tOther1\tOther2\n";
while(<DATA>)
{
chomp;
@data=split /\t/,$_;
@info=split /:/,$data[2];
@posi=split/-/,$info[1];
@end=split / /,$posi[1];
print OUT "$info[0]\t$data[0]\t$posi[0]\t$end[0]\t$data[4]\t$data[3]\t$data[5]\t$data[6]\n";
}
