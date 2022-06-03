sub calcgc {
 my $seq = $_[0];
 my @seqarray = split('',$seq);
 my $count = 0;
 foreach my $base (@seqarray) {
   $count++ if $base =~ /[G|C|g|c]/i;
 }
 my $len = $#seqarray+1;
 my $num=$count / $len;
 #my ($dec)=$num =~ /(\S{6})/;
 print "$seq\t@seqarray\t$count\t$len\t$num\t$dec\n";
 return $num;
}


$se='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCUAACCCU';
$va=calcgc $se;
print "$se\n";
print "$va\n";


