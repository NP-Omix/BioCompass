#!/bin/bash

# A. Inputs

start_time=`date +%s`

[ -z "$1" ] && echo "Missing three letters strain name, please input it right after script name" && exit 1

ls ../antiSMASH_input/$1/*.gbk | cat > clusters.txt

num=`cat clusters.txt | wc -l`

rm clusters.txt

#if [ -d ../outputs/ ]; then rm -r ../outputs/; mkdir ../outputs/; else mkdir ../outputs/; fi

#if [ -d ../outputs/tables/ ]; then rm -r ../outputs/tables/; mkdir ../outputs/tables/; mv $1_* ../outputs/tables/; else mkdir ../outputs/tables/; mv $1_* ../outputs/tables/; fi



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