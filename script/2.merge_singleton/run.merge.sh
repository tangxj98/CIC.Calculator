# a wrapper used for HPC cluster job submission.
if [ $# != 4 ] 
then
	echo "Usage: ./run.merge.R <input file> <index1 sequence> <index2 sequence> <out prefix>"
else
	input=$1
	index1=$2
	index2=$3
	name=$4
	Rscript merge.R $input $index1 $index2 $name
fi
