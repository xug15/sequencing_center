open I, "<$ARGV[0]";
#open O, ">";
$head=<I>;
while(<I>){
#chomp;
@data=split/\t/,$_;
$hash{$data[0]}.=$_;

}
close I;
system("mkdir tmp");
foreach(keys(%hash)){
$file=$_;
	open O, ">>tmp/$file";
print O $head;
print O "$hash{$_}";
	close O;
}






