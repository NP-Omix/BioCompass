# BGCNet (24JUN16) Tiago Leao

## A. Inputs

#### Instructions for users:
1. Please, copy your input antiSMASH file (as is) to the folder named "antiSMASH_input";
2. Copy the antiSMASH files for genomes to be included into the database to the folder "database_clusters/costumized"
3. As suggestion, take a look on the file named "subcluster_dictionary.csv". Change classifications as your will;
4. Please, confirm the number of genomes and clusters per genome;
5. Now, select what type of gene cluster would like to network (e.g. only t1PKS or all clusters). PS: be aware that this selection only filters your input, not the database;
6. Please provide the three letter abbreaviation for each genome;
7. Now, sit and wait till the run is done and the filtering options are prompted. Maybe your help is required assigning category for a product absent in the dictionary, maybe not. =)

#### The script `check.py`(to be writen) will:
- [ ] Check how many genomes are being inputed and how many BGCs (.gbk) are per genome, then prompt to the user;
- [ ] Ask if user want to network a specific type of BGCs, if so, prompt all types of BGCs in the input and ask each ones (prompt options with numbers);
- [ ] Check if all reamining clusters in the first genome contaning the same "DEFINITION"? And in the second genome and so on;
- [X] Ask for 3 letters abbreaviation for each genome (e.g. STR), then replace the "ACCESSION" field for each cluster kepping for example "STR_00X". Then, do the same with the file names;
- [X] Rename the folder name to "STR";
- [X] Once everything is okay, it will start running the network with genome 1. After conclusion, it will run for genome 2, and so on.

#### Observations
```
Be aware of the antiSMASH concatenated .gbk file, either delete it or do not count it;
Also, each genome can not have more than 99 BGCs;
Remember that the user selection of BGC type only filters the input, not the database. Take some time to think if filtering the database is necessary.
```
## B. Building table 1

#### The script `table_1_gen.py` will (for each cluster):
- [X] Parse the gbk file for each CDS, collecting info for locus_tag, start, stop, strand and product;
- [X] Create table 1 with the parsed information.

#### The script `category_gen.py` will (for each cluster):
- [ ] Access the file `subcluster_dictionary.csv`, contaning the most updated version of all antiSMASH smCOGs (numbers) and their respective category annotation;
- [X] Use this dict to extend table_1, adding a column with the category for each gene product;
- [ ] If applied, prompt a smCOG gene product that is not in the dict and ask the user for a proper category annotation.
 
#### Observations
```
In the dictionary, only include (and search it by) the smCOG number, avoiding the current problem of one gene product having more than one annotation.
```

## C. Extending table_1 with best hits

#### The script `table_1_extender.py` wiil (for each cluster):
- [X] Take as input the antiSMASH multigeneBLAST file, parse it and create a tmp file for each hit (which should be a cluster from NCBI);
- [ ] Save every NCBI cluster hits (ID and gene loci) into a list named 'database_cluster_IDs' to be used in step E;
- [X] Scan that temp file to find the best hit per gene;
- [X] Extend the table 1 including best hit per gene, ID, coverage and which "cluster" it came from.

#### Observations
```
By using antiSMASH multigeneBLAST file as input, the number of hits increase significantly, since their cutoff was probably very low. This causes the running time of the script to increase, although, still much smaller than running pre-multigeneBLAST (e.g. up to this step took less than 20 min, instead more than one hour);
One problem is that private (unpublished) gene clusters won't be consider for the subcluster division;
Also pay attention that antiSMASH multigeneBLAST does not have cluster IDs, instead, has scaffold IDs. Therefore, if a gene cluster is divided between two scaffolds in the best hit, this will look like two different cluster during DBScan step;
```

#### Table_1 (extended) example:

`PAL_005_table_1`

| BGC     | locus_tag | product                                          | start | stop  | strand | category     | best_hit_%cov | best_hit_%id | best_hit_BGC | best_hit_gene  | best_hit_gene_loc |
|---------|-----------|--------------------------------------------------|-------|-------|--------|--------------|---------------|--------------|--------------|----------------|-------------------|
| PAL_005 | ctg1_888  | None                                             | 285   | 678   | -1     | hypothetical | None          | None         | None         | None           | None              |
| PAL_005 | ctg1_889  | None                                             | 870   | 1818  | -1     | hypothetical | 100           | 95           | GL890970     | LYNGBM3L_65230 | EGJ29298          |
| PAL_005 | ctg1_890  | None                                             | 2193  | 2442  | -1     | hypothetical | None          | None         | None         | None           | None              |
| PAL_005 | ctg1_891  | None                                             | 2938  | 4123  | -1     | hypothetical | 100           | 93           | GL890970     | LYNGBM3L_65260 | EGJ29300          |
| PAL_005 | ctg1_892  | None                                             | 4091  | 7535  | -1     | hypothetical | 100           | 95           | GL890970     | LYNGBM3L_65280 | EGJ29301          |
| PAL_005 | ctg1_893  | None                                             | 7555  | 10150 | -1     | hypothetical | 100           | 93           | GL890970     | LYNGBM3L_65290 | EGJ29302          |
| PAL_005 | ctg1_894  | None                                             | 10250 | 10877 | -1     | hypothetical | 95            | 58           | GL890971     | LYNGBM3L_69450 | EGJ28870          |
| PAL_005 | ctg1_895  | SMCOG1284:transcriptional_regulator,_MarR_family | 12078 | 12513 | -1     | regulatory   | None          | None         | None         | None           | None              |
| PAL_005 | ctg1_896  | None                                             | 12739 | 14209 | 1      | hypothetical | None          | None         | None         | None           | None              |
| PAL_005 | ctg1_897  | None                                             | 14712 | 14937 | -1     | hypothetical | None          | None         | None         | None           | None              |
| PAL_005 | ctg1_898  | SMCOG1004:thioesterase                           | 15086 | 15806 | -1     | biosynthetic | 95            | 39           | GL890970     | LYNGBM3L_67160 | EGJ29281          |


## D. Creating distance matrixes and running DBScan for subcluster divisions (table 2):

#### The script `subcluster_gen.py` wiil (for each cluster):
- [X] Create a distance matrix using table_1's annotations as input, gene vs gene, scoring according to "category", "best_hit_BGC" and "best_hit_gene_loc";
- [X] Use DBScan to divide the matrix into "density clusters", where the eps (the maximum distance between two samples for them to be considered as in the same neighborhood) will start vary from 1 to length of table_1, only keeping the non-repeated itineration;
- [X] Create a python dictionary for each DBScan "density clusters" output, then, convert this into table_2 (one for ech itineration), the final output of this step.

#### Observations
```
The overall category can be defined by the majority of the gene categories OR by a single biosynthetic gene;
Enumeration of categories per cluster got remove, should it be re-included?;
The fact that the matrix highly depends on the rules makes necessary to optimize those rules, which is not executed in the current version.
```

#### Table_2 examples:

`PAL_005_table_2_1`

| BGC | CDSs_count | category     | loci                                                                                                                                                           | subcluster |
|-----|------------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| PAL | 16         | hypothetical | ctg1_888; ctg1_890; ctg1_896; ctg1_897; ctg1_903; ctg1_905; ctg1_906; ctg1_910; ctg1_911; ctg1_912; ctg1_914; ctg1_916; ctg1_918; ctg1_919; ctg1_920; ctg1_922 | A          |
| PAL | 4          | hypothetical | ctg1_889; ctg1_891; ctg1_892; ctg1_893                                                                                                                         | B          |
| PAL | 15         | biosynthetic | ctg1_894; ctg1_895; ctg1_898; ctg1_899; ctg1_900; ctg1_901; ctg1_902; ctg1_904; ctg1_907; ctg1_908; ctg1_909; ctg1_913; ctg1_915; ctg1_917; ctg1_921           | C          |

`PAL_005_table_2_2`

| BGC | CDSs_count | category     | loci                                                                                                                                                           | subcluster |
|-----|------------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| PAL | 16         | hypothetical | ctg1_888; ctg1_890; ctg1_896; ctg1_897; ctg1_903; ctg1_905; ctg1_906; ctg1_910; ctg1_911; ctg1_912; ctg1_914; ctg1_916; ctg1_918; ctg1_919; ctg1_920; ctg1_922 | A          |
| PAL | 4          | hypothetical | ctg1_889; ctg1_891; ctg1_892; ctg1_893                                                                                                                         | B          |
| PAL | 2          | biosynthetic | ctg1_898; ctg1_902                                                                                                                                             | C          |
| PAL | 2          | mobilome     | ctg1_900; ctg1_908                                                                                                                                             | D          |
| PAL | 11         | hypothetical | ctg1_894; ctg1_895; ctg1_899; ctg1_901; ctg1_904; ctg1_907; ctg1_909; ctg1_913; ctg1_915; ctg1_917; ctg1_921                                                   | E          |

## E. Creating multigeneBLAST database:

#### The script `download_clusters_JGI.py`(to be writen) wiil (for each cluster):
- [ ] Use the list 'database_cluster_IDs' created at step C to download all required gene clusters from NCBI (using ID and loci) that are hits for the input;
- [X] Prepare the bash command to create the multigeneBLAST database using the downloaded cluster from JGI, the most updated MiBIG clusters and costumized antiSMASH clusters (any other antiSMASH clusters that the user want to include)

## F. Running multigeneBLAST:

#### The general script will:
- [X] Use table_2 to run run multigeneBLAST seperatly to each subcluster;
- [X] Move the output to the proper folder;

## G. Creating table of edges:

#### The script `edges_gen.py` will (for each cluster):
- [X] Go through each multigeneBLAST result (for each subcluster and itineration), collect scores and hits IDs and output in a table (STR_edges_all_itineration.txt) where each line represent a multigeneBLAST result;

#### The script `keep_best_itineration.py` will:
- [ ] Go through the output of the script above and remove self hits;
- [X] Filter for best itineration only and generate a second file with those itineration (STR_edges_best_itineration.txt);
- [X] Generate table with statistics for itinerations.

## H. Filtering and generating final outputs:

#### The script `filter_edges.py` will:
- [X] Filter the file STR_edges_best_itineration.txt and keep only the edges that respect the cutoffs specified by the user (% matching genes, min multigeneBLAST score, min cumulative BLAST score);
- [ ] Generate files for subcluster visualization (only for best itineration after filtering).
- [X] Generate table of nodes.

## Visualizing outputs at Cytoscape

#### Example of network:



## Final considerations
```
1. It's important to understand that gene similarity score is unilateral. For example, a gene A can be 80% similar to 50% of B, but in the opposite direction this may increase or even substantially decrease, because in the other direction we need to consider the other 50% of B that was neglected at the first time;
2. Given cpnsideration 1, implementing Jaccard index of pFAM domains should be very useful to validate each gene similarity score;
3. Hits within clusters from same genome input are not included yet;
4. We need to work on how to merge networks (at the end) for more than one input. PS: make sure that input 2,3,4... are in the database of input 1 and  vice-versa;
5. At the filtering step, modify the min values to the actual min values obtained. Also, filtering by MultiGeneBLAST can removed edges with good % of matching genes, good cumulative BLAST score but small number of genes;
6. At the network (final output), larger subcluster tend to have higher cumulative and multigeneBLAST scores, therefore, generating thicker edges. Should we think of a way to normalize these scores? No, because subclusters with more genes can have lower multigeneBLAST scores, which is what the network is trying to represent. Instead, would be nice to have a multigeneBLAST score distribution graph, so we could mirror cytoscape visualization with it.
```

