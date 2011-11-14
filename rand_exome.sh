#/bin/sh

cause=$1
cases=$2
count=$3
ratio=$4
annot_file=Annotate-it_1000G_annotations.txt.uniq
Nonsynonymous_COUNT=1793
Synonymous_COUNT=778
Nonsense_COUNT=30
exon_file=exon_list.txt

# placeholder numbers
p=1
id_len=1000
pos_int=1000

for iter in `seq 1 $cases`
do
awk -v rat=$ratio 'BEGIN {srand()} { p=rand(); while(p<0.00089) p=rand(); if(p<$12) { r=rand(); if(rat>r) {print $9 "\t" $12 "\t" 1 "\t" $11-$10+1 "\t" $2-$10+1} } }' > filtered_${cause}_${iter}_${cases}.txt < ${annot_file}.${cause}
for i in `seq 1 $count`
do
	fid=$(echo "($RANDOM/32767)*160887" | bc -l)
	id=`head -${fid/.*} $exon_file | tail -1`
	echo "$id	$p	1	$id_len $pos_int" >> filtered_${cause}_${iter}_${cases}.txt
done
done

cat filtered_${cause}_*_${cases}.txt | sort | awk -f run_sum_2.awk > filtered_${cause}_${cases}.txt

