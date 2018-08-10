#usage: ./pipeline.sh FILE OUTPUT_DIR KMERSIZE NB_INPUT_KMER NB_OUT_KMER
#fetching arguments
filePath=$1
workPath=$2
kmerSize=$3
kCutoff=$4
lowComplexity=$5
# Getting filename without extension
filename=$(basename -- "$filePath")
filename="${filename%.*}"

#creating workfolder
echo "GENERATING WORKSPACE"
if [ ! -d "$workPath" ]; then
    mkdir $workPath
fi

#Creating kmer list
echo "EXTRACTING EXACT KMERS FROM INPUT"
./kExtract $filePath -k $kmerSize -lc $lowComplexity | head -n $kCutoff > $workPath'/'$filename'_kmers_first'$kCutoff'.txt'

#Converting to fasta
echo "CONVERTING KMER LIST TO FASTA"
python3 kExtract.py  $workPath'/'$filename'_kmers_first'$kCutoff'.txt' >  $workPath'/'$filename'_kmers_first'$kCutoff'.fasta'


#Counting kmers at ~ 2err
echo "COUNTING KMERS"
./adaptFinder3  $filePath -kf  $workPath'/'$filename'_kmers_first'$kCutoff'.fasta' -nt 4 -o $workPath'/'$filename'_COMPTAGE.txt'

echo 'Ressearching overlap for 1st kmer'
./kOverlap $workPath'/'$filename'_COMPTAGE.txt'
echo "DONE"