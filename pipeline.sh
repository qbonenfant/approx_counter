#usage: ./pipeline.sh FILE OUTPUT_DIR KMERSIZE NB_INPUT_KMER NB_OUT_KMER
#fetching arguments
filePath=$1
workPath=$2
kmerSize=$3
kCutoff=$4
lowComplexity=$5
kExtract=/home/cube/Documents/Code/levAdapt/kExtract
adaptFinder=/home/cube/Documents/Code/levAdapt/adaptFinder
kOverlap=/home/cube/Documents/Code/levAdapt/kOverlap
directOv=/home/cube/Documents/Code/levAdapt/Scripts/directOverlap.py
convert=/home/cube/Documents/Code/levAdapt/Scripts/convert.py
# Getting filename without extension
filename=$(basename -- "$filePath")
filename="${filename%.*}"

#creating workfolder
echo "CREATING WORKSPACE"
if [ ! -d "$workPath" ]; then
    mkdir $workPath
fi

#Creating kmer list
echo "EXTRACTING EXACT KMERS FROM INPUT"&&
$kExtract $filePath -k $kmerSize -lc $lowComplexity | head -n $((kCutoff*2)) > $workPath'/'$filename'_kmers_first'$kCutoff'.fasta'&&

#Converting to k-mer list
echo "CONVERTING KMER LIST FROM FASTA to TABBED TEXT"&&
python3 $convert $workPath'/'$filename'_kmers_first'$kCutoff'.fasta' >  $workPath'/'$filename'_kmers_first'$kCutoff'.txt'&&

echo "SEARCHING OVERLAPS USING EXACT KMERS"&&
python3 $directOv $workPath'/'$filename'_kmers_first'$kCutoff'.txt' >  $workPath'/'$filename'_DIRECTOV.csv'&&

#Counting kmers at ~ 2err
echo "COUNTING APPROXIMATE KMERS"&&
$adaptFinder  $filePath -kf  $workPath'/'$filename'_kmers_first'$kCutoff'.fasta' -nt 4 -o $workPath'/'$filename'_COMPTAGE.txt'&&

echo 'Ressearching overlap for 1st kmer'&&
$kOverlap $workPath'/'$filename'_COMPTAGE.txt' > $workPath'/'$filename'_KOVERLAP.csv'&&
echo "DONE"
