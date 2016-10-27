=====
Usage
=====

Using BioCompass in your queried genome
####

The current version of BioCompass can network the gene cluster from ONE strain/genome against the gene cluster at `NCBI/GenBank <https://www.ncbi.nlm.nih.gov/genbank/>`_ and `MIBiG databases <http://mibig.secondarymetabolites.org>`_. `Future Implementations`_ will provide the option to analyze multiple strains.

1. Submit your genome (fasta file) to `antiSMASH 3.0 <http://antismash.secondarymetabolites.org>`_. Once the run is completed, "Download all files" and move the zip file to the BioCompass folder created by the `Installation`_ instructions.

2. Descompact the zip file and rename it to a shorter name. Example, if your strain is called "Streptomyces coelicolor A3", use the name "A3" to rename the antiSMASH file.

3. Create a new folder named "antiSMASH_input" and move your antiSMASH files into it (in this example, your antiSMASH files are inside the folder named "A3", move the whole folder into "antiSMASH_input").

4. Download the proper version of `MultiGeneBlast <https://sourceforge.net/projects/multigeneblast/files/>`_ for command line. Move the zip file to the BioCompass folder created by the `Installation`_ instructions and descompact it.

5. Open a terminal window and access the BioCompass directory where the Makefile is located:: 

    cd path/to/BioCompass
    
    cd BioCompass/

6. Run BioCompass using the following command::
    
    make REFNAME = 'NAME_HERE' INPUTDIR='path/to/antiSMASH_input' MULTIGENEBLASTDIR='path/to/multigeneblast_1.1.14_commandline' ALL

In the example we're using here, the command should be something like::
    
    make REFNAME = 'A3' INPUTDIR='/Users/user/Desktop/BioCompass/antiSMASH_input' MULTIGENEBLASTDIR='/Users/user/Desktop/BioCompass/multigeneblast' ALL

Using BioCompass in your project
####

To use BioCompass in a project::

    import BioCompass
