folder=$1
for f in $folder/*
do
    for anl in $f/*
    do
        workfile=$(basename $f)"_"$(basename $anl)"_COMPTAGE.txt"
        echo 'Working on '$workfile
        python3 Scripts/graphOverlap.py $anl/$workfile > null.txt
    done
done
 
