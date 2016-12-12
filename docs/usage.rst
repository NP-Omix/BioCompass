=====
Basic Usage
=====

The current version of BioCompass can network the biosynthetic gene clusters from ONE strain/genome against the biosynthetic gene clusters from `NCBI/GenBank <https://www.ncbi.nlm.nih.gov/genbank/>`_ and `MIBiG databases <http://mibig.secondarymetabolites.org>`_, as well as user's input (unpublished data).

In order to see an example of using BioCompass in a real dataset and how to use it for multiple strain, check out our tutorials.

Using BioCompass in your queried genome at your machine
####

If your installation of BioCompass succeeded, follow the instructions below. If not, please `report an issue on our GitHub repo <https://github.com/NP-Omix/BioCompass/issues>`_ and include the error message you obtained during your installation.

A. Submit your query genome to `antiSMASH <http://antismash.secondarymetabolites.org>`_. Once the run is complete, download results using "Download all results" (top right corner, an icon shaped as a download arrow). Decompress and save the results into a folder named "antiSMASH_input". We advise renaming the subfolder inside antiSMASH_input to a shorter name representing your strain (hereafter generically referred to as "REFNAME");

B. Download and decompress multigeneBLAST for command line;

	For linux::

	wget 'http://downloads.sourceforge.net/project/multigeneblast/multigeneblast_1.1.13_linux64.tar.gz'
	tar -xvzf multigeneblast_1.1.13_linux64.tar.gz

	For mac:

	wget 'https://superb-sea2.dl.sourceforge.net/project/multigeneblast/multigeneblast_1.1.14_macosx_commandline.zip'

	unzip multigeneblast_1.1.14_macosx_commandline.zip -d multigeneblast

C. Now, execute BioCompass using the following comand::

    cd path/to/BioCompass/BioCompass
    make INPUTDIR='path/to/antiSMASH_input' REFNAME='NAME' MULTIGENEBLASTDIR='path/to/multigeneblast' ALL
    

D. Towards the end of the run, three questions will be prompted asking to select the cutoff you would like to use for filtering your results. We advise using a low cutoff first (near the minimum), checking the network diagram and then revisiting this step using the code below (that only reruns this step, not the whole pipeline) to find the best cutoff for your data::
    
    make INPUTDIR='path/to/antiSMASH_input' REFNAME='NAME' MULTIGENEBLASTDIR='path/to/multigeneblast' step_H


E. The table outputs (REFNAME_edges.txt and REFNAME_nodes.txt) can be visualized as a network diagram using `Cytoscape 3.2.1 <http://www.cytoscape.org/download.php>`_.


Using BioCompass in your project
####

To use BioCompass in a project::

    import BioCompass
