#usage: ./pipeline.sh [FILE] [OUTPUT_DIR] [NB_INPUT_KMER]
#fetching arguments
filePath=$1
workPath=$2
kCutoff=$3
kmList=(12 14 16 18 20)
lcList=(100 75 60 50 45)

# Getting filename without extension
filename=$(basename -- "$filePath")
filename="${filename%.*}"

#creating workfolder
echo "GENERATING WORKSPACE"
if [ ! -d "$workPath" ]; then
    mkdir $workPath
fi

for kmerSize in ${kmList[@]};do
    for lowComplexity in ${lcList[@]};do
        currentWorkPath=$workPath'/k'$kmerSize'_lc'$lowComplexity
        currentFileName=$filename'_k'$kmerSize'_lc'$lowComplexity
        if [ ! -d "$currentWorkPath" ]; then
            mkdir $currentWorkPath
        fi
        echo "WORKING ON FILE "$currentFileName
        date +%D-%H:%M:%S > $currentWorkPath'/timelog.txt'
        #Creating kmer list
        echo "EXTRACTING EXACT KMERS FROM INPUT"
        ./kExtract $filePath -k $kmerSize -lc $lowComplexity | head -n $kCutoff > $currentWorkPath'/'$currentFileName'_kmers_first'$kCutoff'.txt'

        #Converting to fasta
        echo "CONVERTING KMER LIST TO FASTA"
        python3 kExtract.py  $currentWorkPath'/'$currentFileName'_kmers_first'$kCutoff'.txt' >  $currentWorkPath'/'$currentFileName'_kmers_first'$kCutoff'.fasta'

        #Counting kmers at ~ 2err
        echo "COUNTING KMERS"
        ./adaptFinder  $filePath -kf  $currentWorkPath'/'$currentFileName'_kmers_first'$kCutoff'.fasta' -nt 4 -o $currentWorkPath'/'$currentFileName'_COMPTAGE.txt'

        echo 'RESSEARCHING ADAPTER FROM 1ST KMER OVERLAPS'
        ./kOverlap $currentWorkPath'/'$currentFileName'_COMPTAGE.txt'
        date +%D-%H:%M:%S >> $currentWorkPath'/timelog.txt'
        echo ""
    done
done
echo "DONE"
