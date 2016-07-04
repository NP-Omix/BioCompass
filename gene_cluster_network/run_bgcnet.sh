#!/bin/bash

# A. Inputs

start_time=`date +%s`

[ -z "$1" ] && echo "Missing three letters strain name, please input it right after script name" && exit 1

ls ../antiSMASH_input/$1/*.gbk | cat > clusters.txt

num=`cat clusters.txt | wc -l`

rm clusters.txt

#if [ -d ../outputs/ ]; then rm -r ../outputs/; mkdir ../outputs/; else mkdir ../outputs/; fi

#if [ -d ../outputs/tables/ ]; then rm -r ../outputs/tables/; mkdir ../outputs/tables/; mv $1_* ../outputs/tables/; else mkdir ../outputs/tables/; mv $1_* ../outputs/tables/; fi

# E. Creating multigeneBLAST database:

cd ../../multigeneblast_1.1.14/

ls -d ../database_clusters/*/ | cat > subjects.txt

command=`cat subjects.txt | awk -v ORS='* '  '{ print $1 }' | sed 's/,$/\n/'`

python makedb.py $1_db $command

rm subjects.txt

# F. Running multigeneBLAST:

if [ -d ./mgb_result/ ] ; then rm -r ./mgb_result/ ; mkdir ./mgb_result/ ; else mkdir ./mgb_result/ ; fi

if [ $num -gt 9 ];
        then
			for i in {1..9}
			do
				ls ../outputs/tables/$1_00${i}_table2_*.csv | cat > itineration.txt
				iti=`cat itineration.txt | wc -l`
				for j in $(seq 1 $iti);
				do
					table2=../outputs/tables/$1_00${i}_table2_${j}.csv
					while IFS= read -r line; do
						case "$line" in 'BGC'*) continue ;; esac
						genes=`echo $line | awk '{print $4}'`
						sub=`echo $line | awk '{print $5}'`
						python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_00${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -genes $genes -out mgb_result/$1_00${i}_${j}_$sub
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
						genes=`echo $line | awk '{print $4}'`
						sub=`echo $line | awk '{print $5}'`
						python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_0${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -genes $genes -out mgb_result/$1_0${i}_${j}_$sub
					done < $table2
				done
				rm itineration.txt
			done
		else
			for i in $(seq 1 $num);
			do
				ls ../outputs/tables/$1_00${i}_table2_*.csv | cat > itineration.txt
				iti=`cat itineration.txt | wc -l`
				for j in $(seq 1 $iti);
				do
					table2=../outputs/tables/$1_00${i}_table2_${j}.csv
					while IFS= read -r line; do
						case "$line" in 'BGC'*) continue ;; esac
						genes=`echo $line | awk '{print $4}'`
						sub=`echo $line | awk '{print $5}'`
						python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_00${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -genes $genes -out mgb_result/$1_00${i}_${j}_$sub
					done < $table2
				done
				rm itineration.txt
			done
fi

mv mgb_result ../outputs/


# G. Creating table of edges:

cd ../outputs/mgb_result/

if [ $num -gt 9 ];
        then
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
		else
			for i in $(seq 1 $num);
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
fi

rm itineration.txt

python ../../bin/keep_best_itineration.py $1

# H. Filtering and generating final outputs:

python ../../bin/filter_edges.py $1

python ../../bin/feature_gen.py $1 $1_edges_filtered.txt

# mv ./$1_table3.csv ../../outputs/tables/
mv ./$1_features.txt ../../outputs/
mv ./$1_edges_all_itineration.txt ../../outputs/
mv ./$1_edges_best_itineration.txt ../../outputs/
mv ./$1_edges_filtered.txt ../../outputs/

end_time=`date +%s`
echo The execution time was `expr $end_time - $start_time`s.