folder=$1
outputFolder=$2
for f in $folder/*
do
    echo 'Working on '$f
    baseNameFile=$(basename $f .fasta)
    mkdir $outputFolder/$baseNameFile
    date +%D-%H:%M:%S > $outputFolder/$baseNameFile/log.txt
    ./pipeline.sh $f $outputFolder/$baseNameFile 16 500 24 >> $outputFolder/$baseNameFile/log.txt
    date +%D-%H:%M:%S >> $outputFolder/$baseNameFile/log.txt
    cp $f $outputFolder/$baseNameFile/$baseNameFile'.fasta'
done