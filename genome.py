#/usr/bin/python
#GenomePy: a functional and simple library for genomic analysis
import copy


class Sequence(string):
    pass


class GenomeSequence(dic):
    """genome sequence class, currently takes input in multi-fasta format."""
    def __init__(self,genome_sequence = None):
        #tries to read "genome_sequence" into string- first as file
        #    location, then file, then string.
        try:
            sequence_string = open(genome_sequence).read()
        except:
            try:
                sequence_string = genome_sequence.read()
            except:
                sequence_string = genome_sequence
        #breaks sequence into sequence blocks (contigs, scaffolds, or chromosomes), adds sequence 
        #   from each block as dictionary entry into self with block name as key.
        for block in sequence_string.split('>')[1:]:
            block = block.split('\n')
            seqid = block[0]
            seqstring = sequence("".join(block[1:]))
            self[seqid] = seqstring
        


class Genome():
    """genome class, which contains sequence and annotations"""
    
