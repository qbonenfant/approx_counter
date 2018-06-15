mkdir randOutput
filename=$1
f=$(basename $filename .fasta)
nbFiles=$2
nbseq=$3

for i in `seq 1 $nbFiles`;
        do
        	python3 Scripts/extractRandSeq.py $filename $nbseq > randOutput/$f'_'$i.fasta
        done
