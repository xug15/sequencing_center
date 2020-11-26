outdir=output

rRNAf(){
[[ -d ${outdir}/a9-summary ]] || mkdir -p ${outdir}/a9-summary
cd ${outdir}/a5-rmrRNA
array=($(ls|grep err))
head='gene'
for i in "${array[@]}";
do
	#echo $i;
	head+=":${i}";
done
echo -e $head > merge.txt
# file with name add .count
# extract the first and second elements.
begin1=${array[0]};
begin2=${array[1]};
# array remove the first and second element.
array2=("${array[@]:2}");
# join the frist and second file.
join -t ':' ${begin1} ${begin2} > merge.tmp
# use the loop to join each file into the merge file.
for i in ${array2[@]};
do
#echo ${i};
join -t ':' merge.tmp ${i} >>merge.tmp2;
mv merge.tmp2 merge.tmp
done
# merge header and
cat merge.txt merge.tmp > merge2.tmp;
# delete the merge.tmp file.
rm merge.tmp merge.txt
# rename the file.
mv merge2.tmp rRNA.tsv
# replace plates with blanks.
sed -i 's/:/\t/g' rRNA.tsv
cd -

mv ${outdir}/a5-rmrRNA/rRNA.tsv  ${outdir}/a9-summary/

}

#rRNAf

mapf(){
[[ -d ${outdir}/a6-map/b1-summary ]] || mkdir -p ${outdir}/a6-map/b1-summary
cd ${outdir}/a6-map
for i in `ls|grep .final.out`;
do echo $i
	name=${i/Log.final.out/}
sed -n '9,11'p $i >b1-summary/$name.txt
sed -n '24,27'p $i >> b1-summary/$name.txt
sed -i 's/\t//g' b1-summary/$name.txt
done

cd b1-summary

array=($(ls|grep txt))


#echo ${array[@]}

head='gene'

for i in "${array[@]}";
do
	#echo $i;
	head+="|${i}";
done

echo -e $head > merge.txt
# file with name add .count
# extract the first and second elements.
begin1=${array[0]};
begin2=${array[1]};
# array remove the first and second element.
array2=("${array[@]:2}");
# join the frist and second file.
join -t '|' ${begin1} ${begin2} > merge.tmp
# use the loop to join each file into the merge file.
for i in ${array2[@]};
do
#echo ${i};
join -t '|' merge.tmp ${i} >>merge.tmp2;
mv merge.tmp2 merge.tmp
done
# merge header and
cat merge.txt merge.tmp > merge2.tmp;
# delete the merge.tmp file.
rm merge.tmp merge.txt
# rename the file.
mv merge2.tmp map.tsv
# replace plates with blanks.
sed -i 's/|/\t/g' map.tsv
mv map.tsv ../../a9-summary/
#rm *txt
cd -

}

mapf



