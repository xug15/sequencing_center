open DATA, '<Mus_musculus.GRCm38.87.gtf';
open OUT, ">gene_id_transcript_id.txt";

while(<DATA>)
{
@data=split/\t/,$_;

if($data[8]=~/gene_id/ )
{
    if($data[8]=~/transcript_id/)
    {
        $data[8]=~/gene_id "(.*?)"/;
    $gene_id=$1;
    $data[8]=~/transcript_id "(.*?)"/;
    $transcript_id=$1;
    $key="$gene_id\t$transcript_id";
    $hash{$key}='';

    }
    
    }
}

@data=keys(%hash);
foreach(@data)
{
print OUT "$_\n";
@gene=split '\t',$_;

if(!exists($uniq{$gene[0]}))
{
$uniq{$gene[0]}=$gene[1];
}else{
   $uniq{$gene[0]}="\n".$gene[1]; 
}

}
foreach(keys(%uniq))
{
    print "$_\n$uniq{$_}\n";
}
