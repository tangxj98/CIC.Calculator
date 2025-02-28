#!/usr/bin/bash
sample_id=$1  # pure numeric
path=/research/bsi/projects/PI/tertiary/Kannan_Nagarajan_m161624/s214639.CIC_calculator/processing/batch_Feb2023
file=$( ls $path/FastQ/BC_Seq$sample_id/*_1.fq.gz)
python $path/index_extact_and_count.pe.py $file  > $path/results/BC_Seq${sample_id}.temp 
cat $path/results/BC_Seq${sample_id}.temp | grep "," > $path/results/BC_Seq${sample_id}.counts.txt
cat $path/results/BC_Seq${sample_id}.temp | grep -v "," | gzip > $path/results/BC_Seq${sample_id}.missing.gz
rm $path/results/BC_Seq${sample_id}.temp
