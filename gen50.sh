folder=$1
for i  in $(seq 1 50)
do
    python3 patternGen.py > $folder'/Generated_set'$i'.fasta'
done 
