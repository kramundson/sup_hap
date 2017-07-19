dirlist=(*.fastq.gz)
files=$(ls *.fastq.gz | wc -l)
count=$(($files - 2))
for i in `seq 0 2 $count`
do
    j=$(($i + 1))
    echo $i
    echo $j
    BASE=$(basename ${dirlist[$i]} _1.fastq.gz)
    echo $BASE
    echo ${dirlist[$i]}
    echo ${dirlist[$j]}
    paste <(gunzip -c ${dirlist[$i]} | paste - - - - ) <(gunzip -c  ${dirlist[$j]} | paste - - - -) | tr '\t' '\n' | gzip -c > $BASE'.fastq.gz'
done
