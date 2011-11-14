#/bin/sh

cause=$1
cases=$2
count=$3
ratio=$4
annot_file=Annotate-it_1000G_annotations.txt.uniq
Nonsynonymous_COUNT=1793
Synonymous_COUNT=778
Nonsense_COUNT=30
gene_file=gene_list.txt

# placeholder numbers
p=1
id_len=1000
pos_int=1000

for iter in `seq 1 $cases`
do
awk -v rat=$ratio 'BEGIN {srand()} { p=rand(); while(p<0.00089) p=rand(); if(p<$12) { r=rand(); if(rat>r) {print $6 "\t" $12 "\t" 1 "\t" $8-$7+1 "\t" $2-$7+1} } }' > filtered_gene_${cause}_${iter}_${cases}.txt < ${annot_file}.${cause}
for i in `seq 1 $count`
do
	fid=$RANDOM
	if [[ $fid -gt 19921 ]]
	then
		let "fid=$fid-19921"
	fi
	if [[ $fid -eq 0 ]]
	then
		fid=1
	fi
	id=`head -${fid} $gene_file | tail -1`
	echo "$id	$p	1	$id_len $pos_int" >> filtered_gene_${cause}_${iter}_${cases}.txt
done
done

cat filtered_gene_${cause}_*_${cases}.txt | sort | awk -f run_sum_2.awk > filtered_gene_${cause}_${cases}.txt

