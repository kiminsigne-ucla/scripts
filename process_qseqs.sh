# This script concatenates all qseq files with the same prefix into one large file and converts it to
# FASTQ format, only keeping those reads that pass filter
data=$1
prefix=$2

cat $data/qseqs/$prefix*qseq.txt | python /data/home/kinsigne/scripts/qseq2fastq.py > $data/$prefix.fastq