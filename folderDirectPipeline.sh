folder=$1
for f in $folder/*
do
    echo 'Working on '$f
    baseFolder=$(basename $f)
    python3 kDirectOverlap.py  $f'/'$baseFolder'_kmers_first500.txt' > $f'/'$baseFolder'_kDIRECT.txt'
done