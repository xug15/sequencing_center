#!/bin/bash

results=../results

mkdir -p $results/FeatureAnalysis

RiboDensityAtEachKindAAOrCodon -f ../data/GSE89704/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/FeatureAnalysis/up -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
RiboDensityAtEachKindAAOrCodon -f ../data/GSE89704/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/FeatureAnalysis/down -M RPKM -S ../results/down_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
RiboDensityAtEachKindAAOrCodon -f ../data/GSE89704/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/FeatureAnalysis/unblocked -M RPKM -S ../results/unblocked_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/up_all_codon_density.txt -o ../results/FeatureAnalysis/up -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --level codon
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/up_all_codon_density.txt -o ../results/FeatureAnalysis/up -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --level AA
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/down_all_codon_density.txt -o ../results/FeatureAnalysis/down -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --level codon
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/down_all_codon_density.txt -o ../results/FeatureAnalysis/down -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --level AA

PausingScore -f ../data/GSE89704/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/FeatureAnalysis/up -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
ProcessPausingScore -i ../results/FeatureAnalysis/up_si-Ctrl-1_pausing_score.txt,../results/FeatureAnalysis/up_si-Ctrl-2_pausing_score.txt,../results/FeatureAnalysis/up_si-eIF5A-1_pausing_score.txt,../results/FeatureAnalysis/up_si-eIF5A-2_pausing_score.txt -o ../results/FeatureAnalysis/up -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --mode raw --ratio_filter 2 --pausing_score_filter 10

conda activate python27
./seq2logo-2.1/Seq2Logo.py -f ../results/FeatureAnalysis/up_pwm.txt  -u probability -I 5 -o ../results/FeatureAnalysis/up --format PDF
conda deactivate

RiboDensityAroundTripleteAAMotifs -f ../data/GSE89704/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/FeatureAnalysis/up_DD -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa --type2 DDD --type1 DD
RiboDensityAroundTripleteAAMotifs -f ../data/GSE89704/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/FeatureAnalysis/up_PP -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa --type2 PPP --type1 PP
RiboDensityAroundTripleteAAMotifs -f ../data/GSE89704/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/FeatureAnalysis/up_KK -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa --type2 KKK --type1 KK

PlotRiboDensityAroundTriAAMotifs -i ../results/FeatureAnalysis/up_PP_motifDensity_dataframe.txt -o ../results/FeatureAnalysis/up_PPP -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --mode mean --ymax 0.2
PlotRiboDensityAroundTriAAMotifs -i ../results/FeatureAnalysis/up_DD_motifDensity_dataframe.txt -o ../results/FeatureAnalysis/up_DDD -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --mode mean --ymax 0.2
PlotRiboDensityAroundTriAAMotifs -i ../results/FeatureAnalysis/up_KK_motifDensity_dataframe.txt -o ../results/FeatureAnalysis/up_KKK -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 --mode mean --ymax 0.2

bash GetCdsSeqForGeneSets.sh
tAI -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -N ../data/GSE89704/tRNA_GCNs_Saccharomyces_cerevisiae.txt -o ../results/FeatureAnalysis/eIF5A -u 0 -d 500 -t up-regulated-genes,down-regulated-genes,unblocked-genes
cAI -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -o  ../results/FeatureAnalysis/eIF5A -u 0 -d 500 -t up-regulated-genes,down-regulated-genes,unblocked-genes --reference ../results/Ref/longest_cds_sequences.fa

hydropathyCharge -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -t up-regulated-genes,down-regulated-genes,unblocked-genes -o ../results/FeatureAnalysis/eIF5A_hydropathy -u 0 -d 500 --index ../data/GSE89704/hydropathy.txt
hydropathyCharge -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -t up-regulated-genes,down-regulated-genes,unblocked-genes -o ../results/FeatureAnalysis/eIF5A_Charge -u 0 -d 500 --index ../data/GSE89704/AA_charge.txt

PlotHydropathyCharge -i ../results/FeatureAnalysis/eIF5A_Charge_values_dataframe.txt -o ../results/FeatureAnalysis/eIF5A_charge -u 0 -d 500 --mode all --ylab "Average Charges"
PlotHydropathyCharge -i ../results/FeatureAnalysis/eIF5A_hydropathy_values_dataframe.txt -o ../results/FeatureAnalysis/eIF5A_hydro -u 0 -d 500 --mode all --ylab "Average Hydrophobicity"
cAIPlot -i ../results/FeatureAnalysis/eIF5A_local_cAI_dataframe.txt -o ../results/FeatureAnalysis/eIF5A_cAI -u 0 -d 500 --mode all --start 5 --window 7 --step 1
tAIPlot -i ../results/FeatureAnalysis/eIF5A_tAI_dataframe.txt -o ../results/FeatureAnalysis/eIF5A_tAI -u 0 -d 500 --mode all --start 5 --window 7 --step 1