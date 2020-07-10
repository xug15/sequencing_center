
#!/bin/bash

results=../results

mkdir -p $results/EnrichmentAnalysis/GUS $results/EnrichmentAnalysis/MES

RiboDensityAtEachPosition -c ../results/Ref/longest.transcripts.info.extended.txt -f ../data/GSE116570/GUS/attributes.txt -o $results/EnrichmentAnalysis/GUS/GUS  -U codon
RiboDensityAtEachPosition -c ../results/Ref/longest.transcripts.info.extended.txt -f ../data/GSE116570/MES/attributes.txt -o $results/EnrichmentAnalysis/MES/MES -U codon

enrichmentMeanDensity -i $results/EnrichmentAnalysis/GUS/GUS_GUS1-IP-1_cds_codon_density.txt,$results/EnrichmentAnalysis/GUS/GUS_GUS1-IP-2_cds_codon_density.txt -o $results/EnrichmentAnalysis/GUS/GUS_IP
enrichmentMeanDensity -i $results/EnrichmentAnalysis/GUS/GUS_GUS1-total-1_cds_codon_density.txt,$results/EnrichmentAnalysis/GUS/GUS_GUS1-total-2_cds_codon_density.txt -o $results/EnrichmentAnalysis/GUS/GUS_total

enrichmentMeanDensity -i $results/EnrichmentAnalysis/MES/MES_MES1-IP-1_cds_codon_density.txt,$results/EnrichmentAnalysis/MES/MES_MES1-IP-2_cds_codon_density.txt -o $results/EnrichmentAnalysis/MES/MES_IP
enrichmentMeanDensity -i $results/EnrichmentAnalysis/MES/MES_MES1-total-1_cds_codon_density.txt,$results/EnrichmentAnalysis/MES/MES_MES1-total-2_cds_codon_density.txt -o $results/EnrichmentAnalysis/MES/MES_total



EnrichmentAnalysis --ctrl $results/EnrichmentAnalysis/GUS/GUS_total_mean_density.txt --treat $results/EnrichmentAnalysis/GUS/GUS_IP_mean_density.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o $results/EnrichmentAnalysis/GUS/GUS -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500

EnrichmentAnalysis --ctrl $results/EnrichmentAnalysis/MES/MES_total_mean_density.txt --treat $results/EnrichmentAnalysis/MES/MES_IP_mean_density.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o $results/EnrichmentAnalysis/MES/MES -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500

EnrichmentAnalysisForSingleTrans -i $results/EnrichmentAnalysis/MES/MES_codon_ratio.txt -s GUS1 -o $results/EnrichmentAnalysis/MES/MES_GUS1 -c ../results/Ref/longest.transcripts.info.extended.txt  --id-type gene_name --slide-window y --axhline 1
EnrichmentAnalysisForSingleTrans -i $results/EnrichmentAnalysis/MES/MES_codon_ratio.txt -s ARC1 -o $results/EnrichmentAnalysis/MES/MES_ARC1 -c ../results/Ref/longest.transcripts.info.extended.txt  --id-type gene_name --slide-window y --axhline 1

EnrichmentAnalysisForSingleTrans -i $results/EnrichmentAnalysis/GUS/GUS_codon_ratio.txt -s MES1 -o $results/EnrichmentAnalysis/GUS/GUS_MES1 -c ../results/Ref/longest.transcripts.info.extended.txt  --id-type gene_name --slide-window y --axhline 1
EnrichmentAnalysisForSingleTrans -i $results/EnrichmentAnalysis/GUS/GUS_codon_ratio.txt -s ARC1 -o $results/EnrichmentAnalysis/GUS/GUS_ARC1 -c ../results/Ref/longest.transcripts.info.extended.txt  --id-type gene_name --slide-window y --axhline 1


