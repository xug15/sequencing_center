
docker run -v /home/xugang/data/rnastructure/a4-human/output/a0_annotation:/home/xugang/data/rnastructure/a4-human/output/a0_annotation  fjossinet/assemble2 bash -c 'RNAfold -p 30 < /home/xugang/data/rnastructure/a4-human/output/a0_annotation/transcripts_sequence.fa  > /home/xugang/data/rnastructure/a4-human/output/a0_annotation/mfe.tsv'





