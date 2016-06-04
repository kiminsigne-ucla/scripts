input=$1
# ref="reference_sequences_sbfI.fasta"
ref=$2

echo "Cleaning up extra tags from file..."
python /data/home/kinsigne/lib/python2.7/sam_clean.py $input

ext=".sam"
basename=${input%$ext}

echo "Diff'ing file..."
input="${basename}_clean.sam"
output="${basename}_diff.sam"
# generate file to find mismatches and insertions
# in alignments
samtools calmd -S -e $input $ref > $output
rm $input

echo "Sorting sam file..."
# sort to make extracting relevant reads quicker

input="${basename}_diff.sam"
output="${basename}_diff_sorted.bam"
samtools view -bS $input | samtools sort -o $output
rm $input

output="${basename}_diff_sorted"
input="${output}.bam"
output="${output}.sam"
# convert bam to sam
samtools view -h $input > $output
rm $input
