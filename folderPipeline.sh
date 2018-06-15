folder=$1
outputFolder=$2
mkdir done
for f in $folder/*
do
    echo 'Working on '$f
    baseFolder=$(basename $f .fasta)
    time ./pipeline.sh $f $outputFolder/$baseFolder 16 500 10000
    mv $f 'done/'$baseFolder
done
