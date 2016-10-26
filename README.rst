====================
BioCompass
====================


.. image:: https://img.shields.io/pypi/v/gene_cluster_network.svg
        :target: https://pypi.python.org/pypi/gene_cluster_network

.. image:: https://img.shields.io/travis/castelao/gene_cluster_network.svg
        :target: https://travis-ci.org/castelao/gene_cluster_network

.. image:: https://readthedocs.org/projects/gene-cluster-network/badge/?version=latest
        :target: https://gene-cluster-network.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Python package for gene clustering


* Free software: BSD license
* Documentation: https://BioCompass.readthedocs.io.


What is BioCompass?
-------------------

An important skill from bottom-up approaches is to de-replicate NP at the genome level. This de-replication consists on grouping the biosynthetic gene clusters (BGCs) into families according to nucleotide sequence homology. This procedure is known as ‘gene cluster networking’. Since the source code from published networking approach still not publicly available, we adapted our own strategy for the discovery of gene cluster families, referred as BioCompass. 


How BioCompass Works?
---------------------

BioCompass groups BGCs into gene cluster families based on synteny and homology. These cluster need to be identified by antiSMASH[link here] first, ClusterFinder option preferentially turned off. A similarity matrix is used to divide each given BGC into subclusters based on synteny at best MultiGeneBLAST[link here] hits (obtained using antiSMASH 3.0) and the functional annotation of each gene in the queried cluster. This information is then incorporated into a query-specific database to search for the best matches for each subcluster. This newly created database includes microbial BGCs identified by antiSMASH (downloaded from NCBI database, Genbank NR [link here]) and the latest repository for the well database of known gene clusters (MIBiG)[link here]. Additional gene clusters, for example missing from NCBI and MIBiG, can also be inputted by the user. Final similarity scores are calculated via MultiGeneBLAST for each subcluster and stored as tables.

Future Implementations
----------------------

One of the problems of using network approaches at the field of natural products concerns the concept of networking itself. For accurate dereplication of families (both molecular families, like used in GnPS, and gene cluster families), it’s required to define a threshold for which once trespassed, two gene clusters are not part of the same family anymore. Analogously to species definition via 16S rRNA gene, this threshold is empirical and can be imprecise in some cases (reference here). Hence, BioCompass envision to implement a cutoff calibration feature to minimize this issue. This new feature consists on the user evaluating the network diagram for both gene homology (scored via multigene blast) and domain homology (score via Jaccard index, another feature to be implemented soon), visually deciding which cutoff would better represent those scores for his query. The user will count with an internal standard to aid its decision making process.
