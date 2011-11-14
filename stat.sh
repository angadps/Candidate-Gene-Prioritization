#/bin/sh

exon_id=$1
cause=$2
n_cases=$3
var_cases=$4
cases=$5
variants=$6
tot_ctrls=1119
annot_file=Annotate-it_1000G_annotations.txt.uniq.cat
exon_file=exon_list.txt

var_ctrls=`grep -w $exon_id $annot_file | grep -c $cause`
let dorm_ctrls="$tot_ctrls-$var_ctrls"

rm -f input/${exon_id}_${cases}_${variants}.txt
for i in `seq 1 $var_ctrls`
do
	echo "0	1" >> input/${exon_id}_${cases}_${variants}.txt
done
for i in `seq 1 $dorm_ctrls`
do
	echo "0	0" >> input/${exon_id}_${cases}_${variants}.txt
done
for i in `seq 1 $n_cases`
do
	n=`grep -wc $exon_id filtered_${cause}_${i}_${cases}_${variants}.txt`
	# assuming #causative variants(var_cases) = #cases(n_cases)
	if [ $var_cases -gt 0 ]
	then
		let n="$n+1"
		let var_cases="$var_cases-1"
	fi
	echo "1	$n" >> input/${exon_id}_${cases}_${variants}.txt
done

