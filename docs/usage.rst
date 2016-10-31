=====
Usage
=====

The current version of BioCompass can network the gene clusters from ONE strain/genome against the gene clusters from `NCBI/GenBank <https://www.ncbi.nlm.nih.gov/genbank/>`_ and `MIBiG databases <http://mibig.secondarymetabolites.org>`_. Future Implementations will provide the option to analyze multiple strains.

Using BioCompass in your queried genome
####

If your installation of BioCompass succeeded. Follow the instructions below. If not, please `report an issue on our GitHub repo <https://github.com/NP-Omix/BioCompass/issues>`_ including the error message you obtained during your installation.

1. Submit your query genome to `antiSMASH <http://antismash.secondarymetabolites.org>`_ and "Download all results" (top right corner, icon shaped as a download arrow);

2. Once you downloaded the files, descompress and save them into a folder named for example "antiSMASH_input". We advise renaming the subfolder inside antiSMASH_input to a shorter name representing your strain. As an example for this tutorial, since we're using the genome of `Moorea producens PAL 15AUG08-1 <https://www.ncbi.nlm.nih.gov/assembly/GCA_001767235.1>`_, we'll rename the folder as "PAL";


3. Download `multigeneBLAST 1.1.14 <https://sourceforge.net/projects/multigeneblast/files/>`_ for command line;

4. Now, execute BioCompass using the following comand::

    cd path/to/BioCompass/BioCompass
    make INPUTDIR='path/to/antiSMASH_input' REFNAME='NAME' MULTIGENEBLASTDIR='path/to/multigeneblast' ALL
    
For the example using Moorea producens PAL 15AUG08-1, assuming that both antiSMASH_input and multigeneBLAST are in the main BioCompass folder, the command should be::

    make INPUTDIR='/home/Desktop/BioCompass/antiSMASH_input' REFNAME='PAL' MULTIGENEBLASTDIR='/home/Desktop/BioCompass/multigeneblast_1.1.14_macosx_commandline' ALL

PS: the folder name inside antiSMASH_input and the REFNAME must be the same!

5. Select the cutoff you would like to use for filtering your results. We advise using a low cutoff first, checking the network diagram and then revisiting this step using the code (only rerun this step, not the whole pipeline) below to better define a good cutoff for your data.::
    
    make ???


6. The table outputs (REFNAME_edges.txt and REFNAME_nodes.txt) can be vizualized as a network diagram using `Cytoscape 3.2.1 <http://www.cytoscape.org/download.php>`_, according to the tutorial here (tutorial under construction!)

Using BioCompass in your project
####

To use BioCompass in a project::

    import BioCompass
