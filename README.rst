====================
BioCompass
====================


.. image:: https://img.shields.io/pypi/v/BioCompass.svg
        :target: https://pypi.python.org/pypi/BioCompass

.. image:: https://img.shields.io/travis/castelao/BioCompass.svg
        :target: https://travis-ci.org/castelao/BioCompass

.. image:: https://readthedocs.org/projects/biocompass/badge/?version=latest
        :target: https://biocompass.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/travis/castelao/BioCompass.svg
        :target: https://travis-ci.org/castelao/BioCompass

.. image:: https://codecov.io/gh/castelao/BioCompass/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/castelao/BioCompass


Python package for gene clustering


* Free software: BSD license
* Documentation: https://BioCompass.readthedocs.io.


What is BioCompass?
-------------------

An emergent need in the field of natural products is to dereplicate biosynthetic pathways at the genomic level. This dereplication consists of grouping the biosynthetic gene clusters (BGCs) into families according to their nucleotide sequence homology, and a procedure known as ‘gene cluster networking’. Since the source code from published networking approaches are still not publicly available, we adapted our own strategy for the discovery of gene cluster families, named Biosynthetic Gene Cluster Comparative Synteny Software (BioCompass).Please note that this is a beta version of BioCompass which is still undergoing final testing before its official release. The website, its software and all content found on it are provided on an “as is” and “as available” basis.



How BioCompass Works?
---------------------

BioCompass groups BGCs into gene cluster families based on synteny and homology. These clusters need to be identified by antiSMASH, with the ClusterFinder option preferentially turned off. A similarity matrix is used to divide each given BGC into subclusters based on synteny with the best MultiGeneBLAST hits (obtained using antiSMASH 3.0) and the functional annotation of each gene in the queried cluster. This information is then incorporated into a query-specific database to search for the best matches for each subcluster. This newly created database includes microbial BGCs identified by antiSMASH (downloaded from NCBI database, Genbank NR) and the latest version of the well annotated MIBiG database of known gene clusters. Additional gene clusters, for example missing from NCBI and MIBiG, can also be added in by the user. Final similarity scores are calculated by  MultiGeneBLAST for each subcluster and then stored as tables. The outputs can be displayed as a network diagram using Cytoscape v3.2.1.



Future Implementations
----------------------

One of the issues of using networking approaches in the natural products research area concerns the concept of networking itself. For accurate dereplication of families (both molecular families, as used in GnPS (link), and gene cluster families), it’s required to define a threshold for which once trespassed, two gene clusters are not part of the same family anymore. Analogously to the species definition when using the 16S rRNA gene, this threshold is empirical and can be imprecise in some cases. Hence, BioCompass envisions to implement a cutoff calibration feature to minimize this issue. The new feature consists on the user evaluating the network diagram for both gene homology (scored via multigene blast) and domain homology (score via Jaccard index, another feature to be implemented soon), visually deciding which cutoff would better represent those scores for the particular query. The user will use  an internal standard to aid in the decision making process.
