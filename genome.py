#/usr/bin/python
#GenomePy: a functional and simple library for genomic analysis
import copy


def read_to_string(potential_file):
    """Tries to read "potential_file" into string- first as file location,
    then file, then string."""
    try:
        output_string = open(potential_file).read()
    except IOError:
        try:
            output_string = potential_file.read()
        except AttributeError:
            output_string = potential_file
    return output_string


def read_gff3(gff3, high_level_features=['gene'], exon_features=['exon','CDS'],other_features = None, annotation_set = False):
    #reads input file location, file, or string
    gff3_list = read_to_string(gff3).split('\n')
    #removes empty lines
    while "" in gff3_list:
        gff3_list.remove("")
    #sets AnnotationSet based on annotation_set, or creates new one if annotation_set = None
    if annotation_set:
        annotations = annotation_set
    else:
        annotations = AnnotationSet()
    #populates AnnotationSet
    active_annotation = None
    for gff_line in gff3_list:
        line_split = gff_line.split('\t')
        seqid = line_split[0]
        source = line_split[1]
        feature_type = line_split[2]
        coords = [int(line_split[3]),int(line_split[4])]
        coords.sort()
        coords = tuple(coords)
        score = line_split[5]
        strand = line_split[6]
        attributes_split = linesplit[8].split(';')
        attribute_dict = {}
        for attribute in attributes_split:
            x = attributes.split('=')
            attribute_dict[x[0]] = x[1]
        ID = attribute_dict[ID]
        if feature_type in high_level_features:
            annotations[ID] = Annotation(seqid = seqid, coords = coords,
                                         annotation_class = feature_type,
                                         ID = ID, strand = strand)
            active_annotation = annotations[ID]
            active_annotation.score = score
            active_annotation.attribute_dict = copy.copy(attribute_dict)
        elif not feature_type in exon_features:
            active_annotation.feature_type = MiscFeature()
            active_annotation.feature_type.ID = ID
            active_annotation.feature_type.coords = coords
            active_annotation.feature_type.score = score
            active_annotation.feature_type.strand = strand
            active_annotation.feature_type.attribute_dict = copy.copy(attribute_dict)
        elif feature_type == exon_features[0]:
            
            
            
            


class Sequence(str):
    """DNA sequence. Will eventually have methods allowing reverse complimenting,
        translating, etc."""
    pass


class ExonList(list):
    pass


class Exon(tuple):
    pass


class CDS(tuple):
    pass

class MiscFeature():
    pass



class GenomeSequence(dict):
    """genome sequence class, currently takes input in multi-fasta format."""
    def __init__(self,genome_sequence = None):
        #reads input file location, file, or string
        sequence_string = read_to_string(genome_sequence)
        #breaks sequence into sequence blocks (contigs, scaffolds, or chromosomes), adds sequence 
        #   from each block as dictionary entry into self with block name as key.
        if sequence_string != None:
            for locus in sequence_string.split('>')[1:]:
                block = block.split('\n')
                seqid = block[0]
                seqstring = Sequence("".join(block[1:]))
                self[seqid] = seqstring


class Annotation():
    """A single annotation on a genome. For example, a gene or a contig"""
    def __init__(self, seqid = None, coords = (), annotation_class = None, ID = None, strand = None):
        self.seqid = seqid
        self.coords=coords
        self.annotation_class = annotation_class
        self.ID = ID
        self.strand = strand


class AnnotationSet(dict):
    """A set of annotations of a single genome."""
    def __init__(self):
        pass



class Genome():
    """genome class, which contains sequence and annotations"""
    def __init__(self,genome_sequence = None, annotations = None):
        self.genome_sequence = GenomeSequence(genome_sequence)
        
