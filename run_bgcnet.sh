# #!/bin/bash

[ -z "$1" ] && echo "Missing three letters strain name, please input it right after script name" && exit 1

#This part of the script will create the table 1 parsing all the genes in each cluster and then annotate each gene's category according to subcluster_dictionary:

ls ../antiSMASH_input/*.gbk | cat > clusters.txt

num=`cat clusters.txt | wc -l`

rm clusters.txt

for i in {1..9}
do
    python ./table_1_gen.py ../antiSMASH_input/*cluster00${i}.gbk $1_00${i}
    python ./category_gen.py $1_00${i}
done

for i in $(seq 10 $num);
do
    python ./table_1_gen.py ../antiSMASH_input/*cluster0${i}.gbk $1_0${i}
    python ./category_gen.py $1_0${i}
done

if [ -d ../outputs/ ]; then rm -r ../outputs/; mkdir ../outputs/; else mkdir ../outputs/; fi

if [ -d ../outputs/tables/ ]; then rm -r ../outputs/tables/; mkdir ../outputs/tables/; mv $1_* ../outputs/tables/; else mkdir ../outputs/tables/; mv $1_* ../outputs/tables/; fi

#In order to run the multigeneBLAST, we need to create a database:

cd ../multigeneblast_1.1.14/

ls -d ../database_clusters/*/ | cat > subjects.txt

command=`cat subjects.txt | awk -v ORS='* '  '{ print $1 }' | sed 's/,$/\n/'`

python makedb.py $1_db $command

rm subjects.txt

#Now, we'll run the multigeneBLAST to detect the best match per gene:

if [ -d ./pre_mgb_result/ ]; then rm -r ./pre_mgb_result/; mkdir ./pre_mgb_result/; else mkdir ./pre_mgb_result/; fi

for i in {1..9}
do
	python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_00${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 90 -from 0 -to 1000000 -out pre_mgb_result/$1_00${i}
done

for i in $(seq 10 $num);
do
    python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_0${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 90 -from 0 -to 1000000 -out pre_mgb_result/$1_0${i}
done

mv pre_mgb_result ../outputs/

#This step will extend table 1:

for i in {1..9}
do
    python ../bin/table_1_extender.py $1_00${i}
done

for i in $(seq 10 $num);
do
    python ../bin/table_1_extender.py $1_0${i}
done

rm temp.txt

#Next, this step will create all possible subcluster divisions (table 2) using distance matrixes and dbscan:

cd ../outputs/tables

for i in {1..9}
do
    python ../../bin/subcluster_gen.py $1_00${i}
done

for i in $(seq 10 $num);
do
    python ../../bin/subcluster_gen.py $1_0${i}
done

#Now, we need to run again multigeneBLAST, but this time using the subcluster divisions. Then, only the subcluster division with the best summed score will be kept:

cd ../../multigeneblast_1.1.14/

if [ -d ./mgb_result/ ] ; then rm -r ./mgb_result/ ; mkdir ./mgb_result/ ; else mkdir ./mgb_result/ ; fi

for i in {1..9}
do
	ls ../outputs/tables/$1_00${i}_table2_*.csv | cat > itineration.txt
	iti=`cat itineration.txt | wc -l`
	for j in $(seq 1 $iti);
	do
		table2=../outputs/tables/$1_00${i}_table2_${j}.csv
		while IFS= read -r line; do
			case "$line" in 'BGC'*) continue ;; esac
			start=`echo $line | awk '{print $4}'`
			stop=`echo $line | awk '{print $5}'`
			sub=`echo $line | awk '{print $6}'`
			python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_00${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -from $start -to $stop -out mgb_result/$1_00${i}_${j}_$sub
		done < $table2
	done
	rm itineration.txt
done

for i in $(seq 10 $num);
do
	ls ../outputs/tables/$1_0${i}_table2_*.csv | cat > itineration.txt
	iti=`cat itineration.txt | wc -l`
	for j in $(seq 1 $iti);
	do
		table2=../outputs/tables/$1_0${i}_table2_${j}.csv
		while IFS= read -r line; do
			case "$line" in 'BGC'*) continue ;; esac
			start=`echo $line | awk '{print $4}'`
			stop=`echo $line | awk '{print $5}'`
			sub=`echo $line | awk '{print $6}'`
			python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_0${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -from $start -to $stop -out mgb_result/$1_0${i}_${j}_$sub
		done < $table2
	done
	rm itineration.txt
done

mv mgb_result ../outputs/

#Next, the script below will obtain the edges table and then filter it to keep only the best itineration for the subclustering

cd ../outputs/mgb_result/

for i in {1..9}
do
	if [ $i -eq 1 ]
	then
		python ../../bin/edges_gen.py $1_001 1
		mv $1_001_edges_1.txt $1_edges_all_itineration.txt
		ls ../tables/$1_00${i}_table2_*.csv | cat > itineration.txt
		iti=`cat itineration.txt | wc -l`
		for j in $(seq 2 $iti);
		do
	    	python ../../bin/edges_gen.py $1_00${i} ${j}
	    	cat $1_00${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
	    	rm $1_00${i}_edges_${j}.txt
		done
	else
		ls ../tables/$1_00${i}_table2_*.csv | cat > itineration.txt
		iti=`cat itineration.txt | wc -l`
		for j in $(seq 1 $iti);
		do
	    	python ../../bin/edges_gen.py $1_00${i} ${j}
	    	cat $1_00${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
	    	rm $1_00${i}_edges_${j}.txt
		done
	fi
done

for i in $(seq 10 $num);
do
	ls ../tables/$1_0${i}_table2_*.csv | cat > itineration.txt
	iti=`cat itineration.txt | wc -l`
	for j in $(seq 1 $iti);
	do
    	python ../../bin/edges_gen.py $1_0${i} ${j}
    	cat $1_0${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
    	rm $1_0${i}_edges_${j}.txt
	done
done

rm itineration.txt

python ../../bin/keep_best_itineration.py $1

#Last, the scripts below will filter the network and also generate the features table

python ../../bin/filter_edges.py $1

python ../../bin/feature_gen.py $1 $1_edges_filtered.txt

mv ./$1_table3.csv ../../outputs/tables/
mv ./$1_features.txt ../../outputs/
mv ./$1_edges_all_itineration.txt ../../outputs/
mv ./$1_edges_best_itineration.txt ../../outputs/
mv ./$1_edges_filtered.txt ../../outputs/
