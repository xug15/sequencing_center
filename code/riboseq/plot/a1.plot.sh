outdir=outdir

[[ -d ${outdir}/b1-sam ]] || mkdir -p ${outdir}/b1-sam


b1sam(){
for i in `ls|grep Aligned.toTranscriptome.out.bam`;
do echo $i;
	samtools view $i > ${outdir}/b1-sam/$i.sam
done
}
#b1sam


b2rename()
{
cd ${outdir}/b1-sam
for i in `ls`;do 
	echo $i;
	IFS='.' read -ra ADDR <<< "$i";
	mv ${i} ${ADDR[0]}.sam;
done
cd -
}
#b2rename

b3extract(){
name=$1
[[ -d ${outdir}/b2-trans ]] || mkdir -p ${outdir}/b2-trans
for i in `ls ${outdir}/b1-sam|grep sam$`;do 
grep ${name}  ${outdir}/b1-sam/$i > ${outdir}/b2-trans/${name}.${i}.sam
done
cat ${outdir}/b2-trans/${name}*R* > ${outdir}/b2-trans/${name}.total.sam


}

#b3extract ENST00000619449

b4psite(){
name=$1
[[ -d ${outdir}/b3-psite/b1-clean ]] || mkdir -p ${outdir}/b3-psite/b1-clean
for i in `ls ${outdir}/b2-trans|grep sam$`;do
	echo $i;
perl a2.psite.pl ${outdir}/b2-trans/$i
mv ${outdir}/b2-trans/${i}.clean.tsv ${outdir}/b3-psite/b1-clean
mv ${outdir}/b2-trans/${i}.sort.tsv ${outdir}/b3-psite
done
}
#b4psite ENST00000619449

b5rplot(){
name=$1
begin=$2
end=$3
[[ -d ${outdir}/b4-pdf/$name ]] || mkdir -p ${outdir}/b4-pdf/$name
for i in `ls ${outdir}/b3-psite |grep tsv$|grep $name`;do 
	echo ${i};
	echo Rscript a3.rplot.R ${outdir}/b3-psite/${i} ${begin} ${end}
	Rscript a3.rplot.R ${outdir}/b3-psite/${i} ${begin} ${end} ${name} 
	mv ${outdir}/b3-psite/${i}*.pdf ${outdir}/b4-pdf/$name
done

}

#b5rplot ENST00000619449 1908 2030 

b6rplot(){
name=$1
begin=$2
end=$3
[[ -d ${outdir}/b5-pdf/$name ]] || mkdir -p ${outdir}/b5-pdf/$name
for i in `ls ${outdir}/b3-psite |grep tsv$|grep $name`;do
        echo ${i};
        echo Rscript a4.rplot.R ${outdir}/b3-psite/${i} ${begin} ${end}
        Rscript a4.rplot.R ${outdir}/b3-psite/${i} ${begin} ${end} ${name}
	echo -e "mv ${outdir}/b3-psite/${i}*.pdf ${outdir}/b5-pdf/$name"
        mv ${outdir}/b3-psite/${i}*.pdf ${outdir}/b5-pdf/$name
done

}

plot_onestep(){
name=$1
start=$2
end=$3
b3extract ${name}
b4psite ${name}
b5rplot ${name} ${start} ${end}
}

#plot_onestep ENST00000619449 1908 2030
#plot_onestep ENST00000499732 867 1055
#plot_onestep ENST00000418747 717 1112
#plot_onestep ENST00000501122 1086 1274 
#plot_onestep ENST00000578497 119 289 
#plot_onestep ENST00000652112 1917 2603

b6rplot ENST00000619449 1908 2030
b6rplot ENST00000499732 867 1055
b6rplot ENST00000418747 717 1112
b6rplot ENST00000501122 1086 1274
b6rplot ENST00000578497 119 289
b6rplot ENST00000652112 1917 2603


