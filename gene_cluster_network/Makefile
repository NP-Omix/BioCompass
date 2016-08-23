
# To run:
# $ make INPUTDIR='path/to/antiSMASH_input'; run_bgcnet.sh

SHELL = /bin/bash

REFNAME = PAL

INPUTDIR=../antiSMASH_input
OUTPUTDIR=../outputs

TABLEDIR=$(OUTPUTDIR)/tables


# ============================================================================
.SECONDARY:

ALL:
	for file in `ls $(INPUTDIR)/$(REFNAME)/*cluster*.gbk`; do $(MAKE) $(TABLEDIR)/$(REFNAME)_`basename $$file | sed -e s/.*cluster// | sed -e s/\.gbk//`_table2_1.csv; done

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
	cd $(TABLEDIR); python ../../gene_cluster_network/subcluster_gen.py $(REFNAME)_$*

	# E. Creating multigeneBLAST database:
	#cd ../multigeneblast/; \
	#	ls -d ../database_clusters/*/ > subjects.txt \
	#	python makedb.py $(REFNAME)_db \
	#	`cat subjects.txt | awk -v ORS='* '  '{ print $$1 }' | sed 's/,\$/\n/'`

	cd ../multigeneblast/; \
		find ../database_clusters -name '*gbk' | awk '{print $0}' ORS=' ' > gbk_list.txt; \
		python makedb.py PAL_db `cat gbk_list.txt` \
		#rm gbk_list.txt
#legacy:
#	./run_bgcnet.sh

