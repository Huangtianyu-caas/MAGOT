![alt tag](https://github.com/biorover/GenomePy/blob/master/GenomePyLogo75p.png) 
# GenomePy

A functional and simple python library for genomics

This project has two principle components: The genome module, which contains many useful functions and classes for genomic analyses; and the genome_tools script, which has lots of useful tools and wrappers for genomic analyses. The genome_tools_config library is actually basically just a list of paths for genome_tools to use. If all the programs that you'll be calling through genome_tools are in the system PATH, then you don't need to modify genome_tools_config and you're good to go. Otherwise, just open genome_tools_config.py in a text editor, add paths to the programs after the "=" signs, and then save in with the same name and in the same directory. Make sure all three files (genome.py, genome_tools.py, and genome_tools_config.py) are in the same directory as each other!

##Installation
Dependencies: Obviously, this program requires python (V2.7). I've tried to use only the most common/default modules- currently just sys, subprocess, and copy. Additionally many tools within genome_tools.py script require other programs and some are essentially just wrappers for other programs. A list of tools:dependancies is below (hopefully up to date, this thing changes pretty fast right now).  
--nucmer_plot: Mummer and associated executables nucmer, dnadiff, and mummerplot; ps2pdf       
--genewise_wrapper: genewise (note, this must be configured properly- at least the default ubuntu apt-get install DOES NOT WORK.
 please make sure that genewise works by itself before filing a bug here)      
--prep4apollo: Depending on how you use this tool, this requires blast (and associated executables makeblastdb and tblastn) and exonerate

Installation: Download/copy all three python modules (genome.py, genome_tools.py, and genome_tools_config.py) to the same directory. If the dependency executables aren't in your $PATH environment variable (i.e. they don't run from the terminal just by typing the executable name), open genome_tools_config.py in a text editor and set the paths to each executable accordingly. If you want to be able to use the genome.py library within python, either copy all three modules to your python library folder or add the directory containing all three modules to your $PYTHONPATH environment variable (in Linux/Unix, either type into the terminal "export PYTHONPATH=$PYTHONPATH;{path to directory containing all three modules}" at the beggining of each terminal session or add this line to your bash profile).

##Running genome_tools
This is just run as a python script, with the tool you want to use given as the first argument and the arguments for this tool given as subsequent arguments. E.g. to use the "prep4apollo" tool to map exons and full length protein sequences to a genome and output a series of Apollo readable gff3 files (warning, the exonerate step here can take a long time if you have many full length peptide sequences!):
```sh
python {path to module directory}/genome_tools.py prep4apollo genome_sequence={path to fasta containing genome sequence (nucleotide obviously)} exon_fasta={path to fasta containing exon peptide sequences} full_length_seqs={path to fasta containing full length peptide sequences}
```

