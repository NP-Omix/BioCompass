

SHELL = /bin/bash

REFNAME = PAL

ROOT=$(shell cd ../; pwd)

SRC=$(ROOT)/gene_cluster_network
INPUTDIR=$(ROOT)/antiSMASH_input
MULTIGENEBLASTDIR=$(ROOT)/multigeneblast
OUTPUTDIR=$(ROOT)/outputs

MBGDIR=$(OUTPUTDIR)/mgb_result

TABLEDIR=$(OUTPUTDIR)/tables
MAKEDB=$(MULTIGENEBLASTDIR)/makedb.py
MULTIGENEBLAST=$(MULTIGENEBLASTDIR)/multigeneblast.py

ClustersDatabase=../database_clusters

CLUSTER_LIST=`find  $(INPUTDIR)/$(REFNAME)/clusterblast -iname 'cluster*txt'`

# ============================================================================
.SECONDARY:

help:
	@echo "To run:"
	@echo "make INPUTDIR='path/to/antiSMASH_input' MULTIGENEBLASTDIR='path/to/multigeneblast' ALL; run_bgcnet.sh"

ALL:
	for file in `ls $(INPUTDIR)/$(REFNAME)/*cluster*.gbk`; do $(MAKE) $(TABLEDIR)/$(REFNAME)_`basename $$file | sed -e s/.*cluster// | sed -e s/\.gbk//`_table2_1.csv; done
	$(MAKE) $(TABLEDIR)/$(REFNAME)_db
	$(MAKE) Run_multigeneBLAST

#$(REFNAME)_%.gbk:
$(TABLEDIR)/$(REFNAME)_%.gbk:
	mkdir -p `dirname $@`
	# B. Building table 1
	# python table_1_gen.py ../antiSMASH_input/PAL/c00001_Moorea_...cluster001.gbk PAL_001
	python table_1_gen.py $(INPUTDIR)/$(REFNAME)/c00001_Moorea_...cluster$*.gbk $(TABLEDIR)/$(REFNAME)_$*
	# python category_gen.py PAL_001
	python category_gen.py $(TABLEDIR)/$(REFNAME)_$*
	# C. Extending table_1 with best hits
	python table_1_extender.py $(REFNAME) $*
	-rm temp.txt

$(TABLEDIR)/$(REFNAME)_%_table2_1.csv: $(TABLEDIR)/$(REFNAME)_%.gbk
	# D. Creating distance matrixes and running DBScan for subcluster divisions (table 2):
	cd $(TABLEDIR); python $(ROOT)/subcluster_gen.py $(REFNAME)_$*

$(ClustersDatabase):
	mkdir -p $@
	for file in $(CLUSTER_LIST); \
		do BioCompass $$file $@; \
		done

$(TABLEDIR)/$(REFNAME)_db: $(ClustersDatabase)
	# E. Creating multigeneBLAST database:
	#
	# cd ../../multigeneblast_1.1.14/
	#
	# ls -d ../database_clusters/*/ | cat > subjects.txt
	#
	# command=`cat subjects.txt | awk -v ORS='* '  '{ print $1 }' | sed 's/,$/\n/'`
	#
	# python makedb.py $1_db $command
	#
	# rm subjects.txt

	mkdir -p `dirname $@`
	export PATH="$(PATH):$(MULTIGENEBLASTDIR)"; \
	python $(MAKEDB) $@ \
		`find ../database_clusters -name '*gbk' | awk '{print $$0}' ORS=' '`


Run_multigeneBLAST: $(TABLEDIR)/$(REFNAME)_db

	# F. Running multigeneBLAST:

	#if [ -d ./mgb_result/ ] ; then rm -r ./mgb_result/ ; mkdir ./mgb_result/ ; else mkdir ./mgb_result/ ; fi
	#
	#if [ $num -gt 9 ];
	#        then
	#                        for i in {1..9}
	#                        do
	#                                ls ../outputs/tables/$1_00${i}_table2_*.csv | cat > itineration.txt
	#                                iti=`cat itineration.txt | wc -l`
	#                                for j in $(seq 1 $iti);
	#                                do
	#                                        table2=../outputs/tables/$1_00${i}_table2_${j}.csv
	#                                        while IFS= read -r line; do
	#                                                case "$line" in 'BGC'*) continue ;; esac
	#                                                genes=`echo $line | awk '{print $4}'`
	#                                                sub=`echo $line | awk '{print $5}'`
	#                                                python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_00${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -genes $genes -out mgb_result/$1_00${i}_${j}_$sub
	#                                        done < $table2
	#                                done
	#                                rm itineration.txt
	#                        done
	#                        for i in $(seq 10 $num);
	#                        do
	#                                ls ../outputs/tables/$1_0${i}_table2_*.csv | cat > itineration.txt
	#                                iti=`cat itineration.txt | wc -l`
	#                                for j in $(seq 1 $iti);
	#                                do
	#                                        table2=../outputs/tables/$1_0${i}_table2_${j}.csv
	#                                        while IFS= read -r line; do
	#                                                case "$line" in 'BGC'*) continue ;; esac
	#                                                genes=`echo $line | awk '{print $4}'`
	#                                                sub=`echo $line | awk '{print $5}'`
	#                                                python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_0${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -genes $genes -out mgb_result/$1_0${i}_${j}_$sub
	#                                        done < $table2
	#                                done
	#                                rm itineration.txt
	#                        done
	#                else
	#                        for i in $(seq 1 $num);
	#                        do
	#                                ls ../outputs/tables/$1_00${i}_table2_*.csv | cat > itineration.txt
	#                                iti=`cat itineration.txt | wc -l`
	#                                for j in $(seq 1 $iti);
	#                                do
	#                                        table2=../outputs/tables/$1_00${i}_table2_${j}.csv
	#                                        while IFS= read -r line; do
	#                                                case "$line" in 'BGC'*) continue ;; esac
	#                                                genes=`echo $line | awk '{print $4}'`
	#                                                sub=`echo $line | awk '{print $5}'`
	#                                                python ./multigeneblast.py -db $1_db -in ../outputs/tables/$1_00${i}.gbk -hitspergene 50 -minseqcov 80 -minpercid 50 -genes $genes -out mgb_result/$1_00${i}_${j}_$sub
	#                                        done < $table2
	#                                done
	#                                rm itineration.txt
	#                        done
	#fi
	#
	#mv mgb_result ../outputs/

	mkdir -p ${OUTPUTDIR}/mgb_result

	num=`ls $(INPUTDIR)/${REFNAME}/*gbk | wc -l`; \
	for i in $$(seq 1 $$num); \
	    do n=`printf "%03d\n" $$i`; \
	    echo $$n; \
	    nclusters=`ls ${TABLEDIR}/${REFNAME}_$${n}_table2_*.csv  | wc -l`; \
	    echo "nclusters: "$${nclusters}; \
	    for j in $$(seq 1 $${nclusters}); \
	        do echo "j: "$$j; \
	        table2=${TABLEDIR}/${REFNAME}_$${n}_table2_$${j}.csv; \
		echo $$j; \
		while IFS= read -r line; \
		    do case "$$line" in 'BGC'*) continue ;; esac; \
		    genes=`echo $$line | awk '{print $$4}'`; \
		    sub=`echo $$line | awk '{print $$5}'`; \
		    python $(MULTIGENEBLAST) \
		        -db $(TABLEDIR)/$(REFNAME)_db \
		        -in ${TABLEDIR}/${REFNAME}_$${n}.gbk \
			-hitspergene 50 -minseqcov 80 -minpercid 50 \
			-genes $${genes} \
			-out ${OUTPUTDIR}/mgb_result/${REFNAME}_$${n}_$${j}_$${sub}; \
	            done < $$table2; \
	        done; \
	    done


step_G:
	## G. Creating table of edges:
	#
	#cd ../outputs/mgb_result/
	#
	#if [ $num -gt 9 ];
	#        then
	#                        for i in {1..9}
	#                        do
	#                                if [ $i -eq 1 ]
	#                                then
	#                                        python ../../bin/edges_gen.py $1_001 1
	#                                        mv $1_001_edges_1.txt $1_edges_all_itineration.txt
	#                                        ls ../tables/$1_00${i}_table2_*.csv | cat > itineration.txt
	#                                        iti=`cat itineration.txt | wc -l`
	#                                        for j in $(seq 2 $iti);
	#                                        do
	#                                        python ../../bin/edges_gen.py $1_00${i} ${j}
	#                                        cat $1_00${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
	#                                        rm $1_00${i}_edges_${j}.txt
	#                                        done
	#                                else
	#                                        ls ../tables/$1_00${i}_table2_*.csv | cat > itineration.txt
	#                                        iti=`cat itineration.txt | wc -l`
	#                                        for j in $(seq 1 $iti);
	#                                        do
	#                                        python ../../bin/edges_gen.py $1_00${i} ${j}
	#                                        cat $1_00${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
	#                                        rm $1_00${i}_edges_${j}.txt
	#                                        done
	#                                fi
	#                        done
	#                        for i in $(seq 10 $num);
	#                        do
	#                                ls ../tables/$1_0${i}_table2_*.csv | cat > itineration.txt
	#                                iti=`cat itineration.txt | wc -l`
	#                                for j in $(seq 1 $iti);
	#                                do
	#                                python ../../bin/edges_gen.py $1_0${i} ${j}
	#                                cat $1_0${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
	#                                rm $1_0${i}_edges_${j}.txt
	#                                done
	#                        done
	#                else
	#                        for i in $(seq 1 $num);
	#                        do
	#                                if [ $i -eq 1 ]
	#                                then
	#                                        python ../../bin/edges_gen.py $1_001 1
	#                                        mv $1_001_edges_1.txt $1_edges_all_itineration.txt
	#                                        ls ../tables/$1_00${i}_table2_*.csv | cat > itineration.txt
	#                                        iti=`cat itineration.txt | wc -l`
	#                                        for j in $(seq 2 $iti);
	#                                        do
	#                                        python ../../bin/edges_gen.py $1_00${i} ${j}
	#                                        cat $1_00${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
	#                                        rm $1_00${i}_edges_${j}.txt
	#                                        done
	#                                else
	#                                        ls ../tables/$1_00${i}_table2_*.csv | cat > itineration.txt
	#                                        iti=`cat itineration.txt | wc -l`
	#                                        for j in $(seq 1 $iti);
	#                                        do
	#                                        python ../../bin/edges_gen.py $1_00${i} ${j}
	#                                        cat $1_00${i}_edges_${j}.txt | sed 1d >> $1_edges_all_itineration.txt
	#                                        rm $1_00${i}_edges_${j}.txt
	#                                        done
	#                                fi
	#                        done
	#fi
	#
	#rm itineration.txt
	#
	#python ../../bin/keep_best_itineration.py $1

	mkdir -p $(MBGDIR)

	cd $(MBGDIR); \
	num=`echo $(GBK_INPUT_LIST) | wc -w`; \
	for i in $$(seq 1 $$num); \
	    do n=`printf "%03d\n" $$i`; \
	    if [ $$i -eq 1 ]; \
	        then \
		    python $(SRC)/edges_gen.py $(REFNAME)_$${n} 1; \
		    mv $(REFNAME)_001_edges_1.txt $(REFNAME)_edges_all_itineration.txt; \
		    ls $(TABLEDIR)/$(REFNAME)_$${n}_table2_*.csv | cat > itineration.txt; \
		    iti=`cat itineration.txt | wc -l`; \
		    for j in $$(seq 2 $$iti); \
		        do \
	                python $(SRC)/edges_gen.py $(REFNAME)_$${n} $${j}; \
		        cat $(REFNAME)_$${n}_edges_$${j}.txt | \
			    sed 1d >> $(REFNAME)_edges_all_itineration.txt; \
		        rm $(REFNAME)_$${n}_edges_$${j}.txt; \
			done; \
		else \
	            ls $(TABLEDIR)/$(REFNAME)_$${n}_table2_*.csv | cat > itineration.txt; \
		    iti=`cat itineration.txt | wc -l`; \
		    for j in $$(seq 1 $$iti); \
		        do \
			python $(SRC)/edges_gen.py $(REFNAME)_$${n} $${j}; \
			cat $(REFNAME)_$${n}_edges_$${j}.txt | \
			    sed 1d >> $(REFNAME)_edges_all_itineration.txt; \
			rm $(REFNAME)_$${n}_edges_$${j}.txt; \
			done; \
		fi; \
	    done

	rm $(MBGDIR)/itineration.txt

	cd $(MBGDIR); \
	python $(SRC)/keep_best_itineration.py $(REFNAME)



#legacy:
#	./run_bgcnet.sh

