# CIC.Calculator

The scripts are used to count and analyze the barcodes sequences. Scripts for generating MM plot and frequency tables are also provided.   

## Requriement
Python3 and R (v>=3.6) are required to run the pipeline. 

## Pipeline
1. Process Fastq file for index counts
2. Merge small clusters with only minimal difference. 
3. Filter to remove tiny clusters that might come from contamination and random barcode replication.
4. Generate MM plot and frequency count table. 

![Example](https://github.com/tangxj98/CIC.Calculator/blob/main/Example/Example.png)
