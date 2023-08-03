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

Please see more details in the publications:
Syed Mohammed Musheer Aalam and others, DNA barcoded competitive clone-initiating cell analysis reveals novel features of metastatic growth in a cancer xenograft model, NAR Cancer, Volume 4, Issue 3, September 2022, zcac022, https://doi.org/10.1093/narcan/zcac022
