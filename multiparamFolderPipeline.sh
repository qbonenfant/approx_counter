folder=$1
outputFolder=$2
mkdir done
for f in $folder/*
do
    echo 'Working on '$f
    baseFolder=$(basename $f .fasta)
    mkdir $outputFolder/$baseFolder
    date +%D-%H:%M:%S > $outputFolder/$baseFolder/'log.txt'
    ./multiparamPipeline.sh $f $outputFolder/$baseFolder 500 >> $outputFolder/$baseFolder/'log.txt'
    date +%D-%H:%M:%S >> $outputFolder/$baseFolder/'log.txt'
    #mv $f 'done/'$baseFolder
done
 
