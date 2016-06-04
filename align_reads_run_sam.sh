ref_name=$1
lib_file=$2
reads_file=$3
echo $ref_name
echo $lib_file
echo $reads_file

# create bowtie2 index
# make directory if it doesn't exist
mkdir -p $ref_name
cd $refname

echo "../$lib_file"
bowtie2-build "../$lib_file" $ref_name

# align reads
ext=".txt"
basename=${reads_file%$ext}
sam_name="${basename}.sam"
echo $sam_name
echo "Running bowtie2..."
bowtie2 -x $ref_name -r $reads_file -S $sam_name

# create samtools index
samtools faidx $lib_file

# process sam file
./process_sam_for_samparse.sh $sam_name $lib_file

# activate python3
source activate py3k

# run samparse
python samparse_for_multiple_ref.py $lib_file $sam_name samparse_output.txt
