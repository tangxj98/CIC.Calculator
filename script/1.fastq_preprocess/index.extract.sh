# a wrapper of index extract and count

if [ $# != 2 ] 
then 
	echo "Usage: ./index.extract.sh <fastq file> <output prefix>"
else 
	input=$1
	name=$2
	python index_extact_and_count.py $input |  awk -v name=$name '{if($0~/,/){print $0 >>  (name ".summary")}else{print $0 | gzip > (name ".missing.gz")} } '  
fi
