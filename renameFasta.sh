folder=$1
for f in $folder/*
do
        mv  $f $f.fasta
done 
