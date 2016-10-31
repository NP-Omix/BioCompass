=====
Usage
=====

The current version of BioCompass can network the gene clusters from ONE strain/genome against the gene clusters from `NCBI/GenBank <https://www.ncbi.nlm.nih.gov/genbank/>`_ and `MIBiG databases <http://mibig.secondarymetabolites.org>`_. Future Implementations will provide the option to analyze multiple strains.

Using BioCompass in your queried genome at your machine
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

PS2: The running time at an average laptop (e.g. MacBook Pro 2.6 GHz Intel Core i5 8 GB 1600 MHz DDR3) for Moorea producens PAL 15AUG08-1 (which contains 44 gene clusters) was XhXmin. We advise the use of the command `screen <https://www.linode.com/docs/networking/ssh/using-gnu-screen-to-manage-persistent-terminal-sessions>`_.

Using BioCompass in your queried genome at `Amazon Web Service <https://aws.amazon.com>`_
####

We expect that tutorial above would work at your machine, but in case it doesn't, we can garantee that it will work on AWS if followed the instructions in this topic. The instructions for running BioCompass at AWS are very similar than explained in the tutorial above, although requires the following setup first.

A. Start an Amazon Web Services EC2 computer following the tutorial `here <http://2016-metagenomics-sio.readthedocs.io/en/latest/aws/boot.html>`_

B. Log into your instance following the instructions `here <http://2016-metagenomics-sio.readthedocs.io/en/latest/aws/login-shell.html>`_

C. Once logged into your AWS computer, install BioCompass running the following command::

    git clone git://github.com/NP-Omix/BioCompass
    sudo apt-get install gcc
    curl https://bootstrap.pypa.io/get-pip.py | sudo python #acquires most updated version of pip
    sudo apt-get install python-dev #to install missing python headers
    sudo pip install BioCompass

D. Upon successful instalation, follow steps 1 to 3 from the tutorial above to generate the antiSMASH and multigeneBLAST folders. Open a new terminal window and Copy the folders from your machine into your AWS instance using the command::

    scp -r -i ~/Downloads/amazon-key.pem path/to/antiSMASH_input/ ubuntu@MACHINE_NAE:/home/
    scp -r -i ~/Downloads/amazon-key.pem path/to/multigeneBLAST_1.1.14_macosx_commandline/ ubuntu@MACHINE_NAE:/home/

E. Follow steps 1 to 3 from the tutorial above to conclude your analysis.

Using BioCompass in your project
####

To use BioCompass in a project::

    import BioCompass
