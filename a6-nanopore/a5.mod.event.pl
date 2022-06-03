open I, "<$ARGV[0]"; # modification table Tetra_NAI_N3_1.combined.RDatamod_mat.tsv
open I2, "<$ARGV[1]"; # read events. Tetra_NAI_N3_1.combined.RData.tsv

$outname=$ARGV[0];
$outname=~s/.combined.RDatamod_mat.tsv/.mod.tsv/g;
print "perl a5.mod.event.pl Tetra_NAI_N3_1.combined.RDatamod_mat.tsv Tetra_NAI_N3_1.combined.RData.tsv $outname\n";
open O, ">$outname"; #output.

# set rna length.

$rnalength=440;
$rnalength=$ARGV[2];


<I>;
while(<I>){
$_=~s/"//g;
@data=split/\t/,$_;
# get read name.
$readn=shift(@data);
#print "$readn\n$#data\n";
$info=join ",",@data;
#print "$info\n";
#get each base modification 0,1,NA as x. array start as 0.
$hash{$readn}=$info;
}

close I;
#set empty array.
@list=();
while(<I2>){
chomp;
@data=split /\t/,$_;
# test the read name is exists in the mod matrix or not. If exists, then next.
if(exists($hash{$data[1]})){
#The content of array that each position modification: 0,1,x.
@mod=split ',',$hash{$data[1]};
$name=$data[1];
#print "$data[1]\n$hash{$data[1]}\n";
# Using list to records the remembered data.
# If the read recoder before, that will create a array named by read's name.
if(exists($list{$name}))
{
	#	print "mod label".$mod[$data[1]]."\n";
# the $data[2] is the position in the RNA, also which is the array id number. They are same.
# get the situation of position . 
	if($mod[$data[2]]=~/0/){
#    
		#${$name}[$data[2]]=$data[5];
#print "$data[2]\t$data[5]\n${$name}[$data[2]]\n";
}
# the situation is 1 means the position is modification.
if($mod[$data[2]]=~/1/)
 {
#  

	 ${$name}[$data[2]]=$data[5];

 }
# if we do not read the read, we will records the data, and create the array named by reads name.
}else{ # 3fdef
$list{$name}='';
@{$name}=();
# we will put all NA as first value.
 for ($i=0;$i<$rnalength;$i++)
  {
  push @{$name}, 'NA';
  }

  #print "mod label".$mod[$data[1]]."\n";
# the $data[2] is the position in the RNA, also which is the array id number. They are same.
# get the situation of position .
        if($mod[$data[2]]=~/0/){
#
		#${$name}[$data[2]]=$data[5];
#print "$data[2]\t$data[5]\n${$name}[$data[2]]\n";
}
# the situation is 1 means the position is modification.
if($mod[$data[2]]=~/1/)
 {
#

         ${$name}[$data[2]]=$data[5];

 }



 }# end with else : 3fdef
}

}

# get the read array.
foreach(keys(%list)){
print O "$_\t";
$info=join "\t",@{$_};
print O $info;
print O "\n";

}


