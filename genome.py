#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas


import copy
import StringIO
from annotation_parse import *
from magot_smallfuncs import *
from magot_variants import *

import numpy

verbose = True


class Genotype():
    """individual genotype to populate a GenotypeDict"""
    def name(self, ):
        pass
    

class Variant():
    """individual variants from reference genome sequence. Can be SNPs or polimorphism"""
    def __init__(self, ref_allele, alt_alleles, genotypes = None, variant_set = None, vcf_header_index = None):
        self.variant_set = variant_set
        self.ref_allele = ref_allele
        self.alt_alleles = alt_alleles
        self.genotypes = genotypes
    

class VariantSet(dict):
    """set of variants from reference genome sequence. Can be SNPs or polimorphism"""
    def __init__(self, genome = None, vcf_headers = None):
        self.genome = genome
        self.vcf_headers = vcf_headers
    

class AnnotationSet():
    """A set of annotations of a single genome. Each feature type (e.g. gene, transcript, exon, etc.)
    is stored in it's own dictionary as Annotations with their ID as their key (see "Annotation" class).
    The AnnotationSet itself also functions losely as a dictionary, in that any feature can be returned
    by indexing the AnnotationSet with the ID as a key (e.g. my_annotation_set["my_feature_ID"])"""
    def __init__(self, genome = None):
        self.gene = {}
        self.transcript = {}
        self.CDS = {}
        self.UTR = {}
        self.genome = genome
    
    def __getitem__(self,item):
        all_dicts = {}
        for attribute in dir(self):
            if type(eval("self." + attribute)) == dict:
                try:
                    all_dicts[item] = eval("self." + attribute)[item]
                except:
                    pass
        return all_dicts[item]
    
    def read_gff(self, gff, *args, **kwargs):
        kwargs["annotation_set_to_modify"] = self
        read_gff(gff, *args, **kwargs)
    
    def get_seqid(self, seqid):
        seqid_annotation_set = AnnotationSet()        
        for attribute in dir(self):
            if type(eval("self." + attribute)) == dict:
                setattr(seqid_annotation_set,attribute,{})
                for feature in eval('self.'+ attribute):
                    feature_obj = eval('self.' + attribute)[feature]
                    if feature_obj.seqid == seqid:
                        eval('seqid_annotation_set.' + feature_obj.feature_type)[feature] = feature_obj
        return seqid_annotation_set
    
    def get_all_seqids(self):
        seqid_list = []
        for attribute in dir(self):
            if type(eval("self." + attribute)) == dict:
                for feature in eval('self.'+ attribute):
                    seqid_list.append(eval('self.' + attribute)[feature].seqid)
        return list(set(seqid_list))
    
    def read_exonerate(self, exonerate_output):
        read_exonerate(exonerate_output,annotation_set_to_modify = self)
    
    def read_blast_csv(self, blast_csv, hierarchy = ['match','match_part'], source = 'blast', find_truncated_locname = False):
        read_blast_csv(blast_csv, annotation_set_to_modify = self, hierarchy = hierarchy, source = source, find_truncated_locname = find_truncated_locname)
    
    def read_cegma_gff(self, cegma_gff):
        read_cegma_gff(cegma_gff, annotation_set_to_modify = self)
    
    def get_fasta(self,feature,seq_type = "nucleotide",longest=False,genomic = False):
        fasta_list = []
        for annotation in eval('self.' + feature):
                fasta_list.append(eval('self.' + feature)[annotation].get_fasta(seq_type = seq_type,longest = longest, genomic = genomic))
        return "\n".join(fasta_list)
        
    

class BaseAnnotation():
    """Bottom-most level annotation on a genome, for example CDS, UTR, Match_part, etc. Anything that should have no children"""
    def __init__(self, ID, seqid, coords, feature_type, parent = None, strand = ".", other_attributes = {}, annotation_set = None):
        #Sets up most attributes
        self.ID = ID
        self.seqid = seqid
        self.coords = coords
        self.feature_type = feature_type
        self.annotation_set = annotation_set
        for attribute in other_attributes:
            setattr(self, attribute, other_attributes[attribute])
        self.parent = parent
        self.strand = strand

    def get_coords(self):
        return self.coords
    
    def get_seq(self):
        try:
            if self.strand == '+' or self.strand == '.':
                return Sequence(self.annotation_set.genome.genome_sequence[self.seqid][self.coords[0]-1:self.coords[1]])
            elif self.strand == '-':
                return Sequence(self.annotation_set.genome.genome_sequence[self.seqid][self.coords[0]-1:self.coords[1]]).reverse_compliment()
            else:
                print self.ID + ' has invalid strand value "' + self.strand + '"'
        except:
            print "either base_annotation has not annotation_set, or annotation_set has no genome, or genome has no\
            genome sequence, or genome sequence has no matching seqid, or coords are out of range on that seqid"
            print self.seqid
    
    def get_gff(self, gff_format = "simple gff3"):
        if self.annotation_set != None:
            fields = ['seqid','source','feature_type','get_coords()[0]','get_coords()[1]','score','strand','phase']
            fields_list = []
            for field in fields:
                try:
                    fields_list.append(str(eval('self.' + field)))
                except AttributeError:
                    fields_list.append('.')
            if gff_format == "simple gff3" or "extended gff3":
                defline = 'ID=' + self.ID
                if self.parent != None:
                    defline = defline + ';Parent=' + self.parent                
            if gff_format == "extended gff3":
                for attribute in self.__dict__:
                    if type(self.__dict__[attribute]).__name__ == 'str' and not attribute in ['ID','Parent','score','strand','seqid','feature_type','phase']:
                        defline = defline + ';' + attribute + '=' + self.__dict__[attribute]
            elif gff_format[:13] == "augustus hint":
                gff_format_fields = gff_format.split()
                fields_list[2] = gff_format_fields[2]
                defline = "src=" + gff_format_fields[3]
                if len(gff_format_fields) > 4:
                    defline = defline + ";pri=" + gff_format_fields[4]
            # elif gff_format[:3] == "gtf":
            #     defline = ""
            #     for parent_type in gff_format.split()[1:]:
            #         defline = defline + parent_type + ''
            #     
            fields_list.append(defline)
            return '\t'.join(fields_list)
    


class ParentAnnotation():
    """Parent of any BaseAnnotation. Examples include genes and transcripts. Suggested hierarchy for genes is
    CDS (as BaseAnnotation) -> transcript -> gene."""
    def __init__(self, ID, seqid, feature_type, child_list = [], parent = None, strand = ".", annotation_set = None, other_attributes = {}):
        self.ID = ID
        self.seqid = seqid
        self.feature_type = feature_type
        self.child_list = copy.copy(child_list)
        self.parent = parent
        self.annotation_set = annotation_set
        self.strand = strand
        for attribute in other_attributes:
            setattr(self, attribute, other_attributes[attribute])
    
    def get_coords(self):
        if len(self.child_list) > 0 and self.annotation_set != None:
            coords_list = []
            for child in self.child_list:
                child_object = self.annotation_set[child]
                if child_object.__class__.__name__ == 'ParentAnnotation':
                    coords_list = coords_list + list(child_object.get_coords())
                elif child_object.__class__.__name__ == 'BaseAnnotation':
                    coords_list = coords_list + list(child_object.coords)
                else:
                    print "for some reason you have children in ParentAnnotation " + self.ID + " which are neither \
                    ParentAnnotation objects nor BaseAnnotation object. Get your act together"
            return (min(coords_list),max(coords_list))
    
    def get_fasta(self, seq_type = "nucleotide", longest=False, genomic = False, name_from = 'ID'):
        """Returns fasta of this annotation's sequence. If this feature has multiple subfeatures (e.g. this is a gene
        and it has multiple transcripts), the sequence of each subfeature will be an entry in the fasta string."""
        if genomic == True:
            if self.annotation_set.genome != None:
                return ">" + self.ID + '\n' + self.annotation_set.genome.genome_sequence[self.seqid][self.get_coords()[0] - 1:self.get_coords()[1]] + '\n'
        elif len(self.child_list) > 0 and self.annotation_set != None:
            if self.annotation_set.genome != None:
                fasta_list = []
                child_type = self.annotation_set[self.child_list[0]].__class__.__name__
                if child_type == 'BaseAnnotation':
                    seq_list = []
                    child_dict = {}
                    for child in self.child_list:
                        child_obj = self.annotation_set[child]
                        try:
                            child_dict[child_obj.coords] = child_obj.get_seq()
                        except AttributeError:
                            print "ParentAnnotation has both ParentAnnotation and BaseAnnotation children!"
                            print self.ID
                        strand = child_obj.strand
                    children_in_correct_order = list(child_dict)
                    children_in_correct_order.sort()
                    if strand == '-':
                        children_in_correct_order.reverse()
                    for child in children_in_correct_order:
                        seq_list.append(child_dict[child])
                    if seq_type == "nucleotide":
                        new_seq = Sequence("".join(seq_list))
                    elif seq_type == "protein":
                        new_seq = Sequence("".join(seq_list)).translate()
                    else:
                        print seq_type + ' is not valid seq_type. Please specify "protein" or "nucleotide".'
                    fasta_list.append('>' + self.__dict__[name_from] + '\n' + new_seq)
                else:
                    for child in self.child_list:
                        try:
                            child_fasta = self.annotation_set[child].get_fasta(seq_type=seq_type, name_from = name_from)
                            if child_fasta != "":
                                fasta_list.append(child_fasta)
                        except AttributeError:
                            print "ParentAnnotation has both ParentAnnotation and BaseAnnotation children!"
                            print self.ID
                if longest == True:
                    seqlens = {}
                    for seq in fasta_list:
                        seqlens[len("".join(seq.split('\n')[1:]))] = seq
                    return seqlens[max(list(seqlens))]
                else:
                    try:
                        return '\n'.join(fasta_list)
                    except:
                        return ""
            else: return ""
        else: return ""
    
    def get_gff(self, gff_format = "simple gff3"):
        """presets will eventually include "simple gff3", "simple gff2", "apollo gff3", and more by request"""
        if self.annotation_set != None:
            fields = ['seqid','source','feature_type','get_coords()[0]','get_coords()[1]','score','strand','phase']
            fields_list = []
            parent_line = True
            for field in fields:
                try:
                    fields_list.append(str(eval('self.' + field)))
                except AttributeError:
                    fields_list.append('.')
            if gff_format == "simple gff3" or "extended gff3":
                defline = 'ID=' + self.ID
                if self.parent != None:
                    defline = defline + ';Parent=' + self.parent
            if gff_format == "extended gff3":
                for attribute in self.__dict__:
                    if type(self.__dict__[attribute]).__name__ == 'str' and not attribute in ['ID','Parent','score','strand','seqid','feature_type','phase']:
                        defline = defline + ';' + attribute + '=' + self.__dict__[attribute]
            elif gff_format[:13] == "augustus hint":
                parent_line = False
            fields_list.append(defline)
            if parent_line:
                lines_list = ['\t'.join(fields_list)]
            else:
                lines_list = []
            child_dict = {}
            for child in self.child_list:
                child_object = self.annotation_set[child]
                if child_object.get_coords() not in child_dict:
                    child_dict[child_object.get_coords()] = child_object.get_gff(gff_format)
                else:
                    child_dict[(child_object.get_coords()[0],child_object.get_coords()[1] + len(child_dict))] = child_object.get_gff(gff_format)
            child_coords = list(child_dict)
            child_coords.sort()
            for child_index in child_coords:
                lines_list.append(child_dict[child_index])
            return '\n'.join(lines_list)


class Sequence(str):
    """DNA sequence. Has methods allowing reverse complimenting,
        translating, etc."""
    def reverse_compliment(self):
        """returns reverse compliment of self"""
        new_sequence_list = []
        compliment_dict = {'a':'t','t':'a','g':'c','c':'g','A':'T','T':'A','G':'C','C':'G','n':'n','N':'N','-':'-'}
        for residue in self[::-1]:
            try:
                new_sequence_list.append(compliment_dict[residue])
            except KeyError:
                new_sequence_list.append('n')
        return Sequence(''.join(new_sequence_list))
    
    def translate(self,library = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                                  'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                                  'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                                  'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                                  'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
                                  'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
                                  'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                                  'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'},
                  frame = 0, strand = '+',trimX = True):
        triplet = ""
        newseq = ""
        if strand == '+':
            seq = self
        elif strand == '-':
            seq = self.reverse_compliment()
        if len(seq) > (2 + frame):
            for residue_position in range(frame, len(self)):
                triplet = triplet + seq[residue_position].upper()
                if (residue_position + frame) % 3 == 2:
                    try:
                        newseq = newseq + library[triplet]
                    except KeyError:
                        newseq = newseq + 'X'
                    triplet = ""
            if trimX:
                if newseq[0] == 'X':
                    newseq = newseq[1:]
            return newseq
    
    def get_orfs(self, longest = False, strand = 'both', from_atg = False):
        orflist = []
        if longest:
            candidate_list = []
            longest_orf_len = 0
        for frame in [0,1,2]:
            for strand in ['-','+']:
                translated_seq = self.translate(frame=frame,strand=strand)
                if translated_seq:
                    translated_seq_list = translated_seq.split('*')
                    for orf in translated_seq_list:
                        if from_atg:
                            try:
                                output_orf = 'M' + ''.join(orf.split('M')[1:])
                            except IndexError:
                                continue
                        else:
                            output_orf = orf
                        if longest:
                            if len(output_orf) > longest_orf_len:
                                candidate_list.append(output_orf)
                                longest_orf_len = len(output_orf)
                        else:
                            orflist.append(output_orf)
        if longest:
            return candidate_list[-1]
        else:
            return orflist


class GenomeSequence(dict):
    """genome sequence class, currently takes input in multi-fasta format."""
    def __init__(self,genome_sequence = None, truncate_names = False):
        #reads input file location, file, or string
        sequence_file = ensure_file(genome_sequence)
        #breaks sequence into sequence blocks (contigs, scaffolds, or chromosomes), adds sequence 
        #   from each block as dictionary entry into self with block name as key.
        if sequence_file != None:
            seq = ""
            seqname = ""
            for line in sequence_file:
                if line[0] == ">":
                    if truncate_names == True:
                        seqid = line[1:].replace('\r','').replace('\n','').split()[0]
                    else:
                        seqid = line[1:].replace('\r','').replace('\n','')
                    if seq != "":
                        self[seqname] = seq
                        seq = ""
                    seqname = seqid
                else:
                    seq = seq + line.replace('\r','').replace('\n','')
            if seq != "":
                self[seqname] = seq


class Genome():
    """genome class, which contains sequence and annotations. Annotations can be given as annotation_set object, gff3, cegma_gff,
    blast_csv, or exonerate_output (just set annotation_format)."""
    def __init__(self,genome_sequence = None, annotations = None, varients = None, annotation_format = 'annotation_set', truncate_names = False):
        if genome_sequence.__class__.__name__ == 'GenomeSequence' or genome_sequence == None:
            self.genome_sequence = genome_sequence
        else:
            self.genome_sequence = GenomeSequence(genome_sequence, truncate_names = truncate_names)
        if annotations != None:
            if annotations.__class__.__name__ == "AnotationSet" and annotation_format == 'annotation_set':
                self.annotations = annotations
                self.annotations.genome = self
            elif annotation_format == 'gff3':
                self.annotations = read_gff(annotations)
                self.annotations.genome = self
            elif annotation_format == 'cegma_gff':
                self.annotations = read_cegma_gff(annotations)
                self.annotations.genome = self
            elif annotation_format == 'blast_csv':
                self.annotations = read_blast_csv(annotations)
                self.annotations.genome = self
            elif annotation_format == 'exonerate_output':
                self.annotations = read_exonerate(annotations)
                self.annotations.genome = self
        else:
            self.annotations = annotations
    
    def get_scaffold_fasta(self, seqid):
        return '>' + seqid + '\n' + self.genome_sequence[seqid]
    
    def get_genome_fasta(self, remove_spaces = False):
        fasta_list = []
        for seqid in self.genome_sequence:
            if remove_spaces:
                fasta_header = seqid.split()[0]
            else:
                fasta_header = seqid
            fasta_list.append('>' + fasta_header + '\n' + self.genome_sequence[seqid])
        return "\n".join(fasta_list)
    
    def write_apollo_gff(self, seqid, suppress_fasta = False):
        if self.genome_sequence != None and self.annotations != None:
            try:
                apollo_gff = write_longform_gff(self.annotations.get_seqid(seqid))
                if not suppress_fasta:
                    apollo_gff = apollo_gff + '\n' + self.get_scaffold_fasta(seqid)
                return apollo_gff
            except:
                if not suppress_fasta:
                    return self.get_scaffold_fasta(seqid)
                else:
                    return ""
        else:
            print "genome object is either missing genome_sequence or annotations"
    
    def get_seqids(self, from_annotations = False):
        seqid_list = []
        warning = False
        if self.genome_sequence != None:
            for seqid in self.genome_sequence:
                seqid_list.append(seqid)
        if self.annotations != None and from_annotations:
            for seqid in self.annotations.get_all_seqids():
                if seqid not in seqid_list:
                    seqid_list.append(seqid)
                    warning = True
        if warning:
            print "warning, some annotations possessed seqids not found in sequence dictionary"
        return seqid_list
    
    def read_exonerate(self, exonerate_output):
        if self.annotations != None:
            self.annotations.read_exonerate(exonerate_output)
        else:
            self.annotations = read_exonerate(exonerate_output)
            self.annotations.genome = self
    
    def read_blast_csv(self, blast_csv, hierarchy = ['match','match_part'], source = 'blast', find_truncated_locname = False):
        if self.annotations == None:
            self.annotations = AnnotationSet()
            self.annotations.genome = self
        self.annotations.read_blast_csv(blast_csv, hierarchy = hierarchy, source = source, find_truncated_locname = find_truncated_locname)
    
    def read_cegma_gff(self, cegma_gff):
        if self.annotations != None:
            self.annotations.read_cegma_gff(cegma_gff)
        else:
            self.annotations = read_cegma_gff(cegma_gff)
            self.annotations.genome = self

    def read_gff(self, gff, *args, **kwargs):
        if self.annotations != None:
            self.annotations.read_gff(gff, *args, **kwargs)
        else:
            self.annotations = read_gff(gff, *args, **kwargs)
            self.annotations.genome = self
    
    def read_vcf(self, vcf):
        self.variants = read_vcf(vcf)
    

class position_dic(dict):
    def __init__(self, genome_sequence, dtype=bool):
        for seqid in genome_sequence:
            self[seqid] = numpy.zeros(len(genome_sequence[seqid]),dtype=dtype)

    def fill_from_annotations(self, annotation_set, feature, fill_type = "coords",fill_with = "1"):
        """accepts "start" and "coords" for fill_type"""
        for annotation in eval("annotation_set." + feature):
            annotation_obj = eval("annotation_set." + feature)[annotation]
            seqid = annotation_obj.seqid
            coords = annotation_obj.get_coords()
            try:
                parent = annotation_obj.parent
            except:
                pass
            ID = annotation_obj.ID
            if fill_type == "coords":
                for position in range(coords[0] - 1,coords[1]):
                    self[seqid][position] = eval(fill_with)
            elif fill_type == "start":
                self[seqid][coords[0] - 1] = eval(fill_with)
    
    def count_from_annotations(self, annotation_set, feature):
        count_list = []
        for annotation in eval("annotation_set." + feature):
            annotation_obj = eval("annotation_set." + feature)[annotation]
            seqid = annotation_obj.seqid
            coords = annotation_obj.get_coords()
            featureID = annotation_obj.ID
            try:
                parent = annotation_obj.parent
            except:
                pass
            ID = annotation_obj.ID
            count0 = 0
            count1 = 0
            for position in range(coords[0] - 1,coords[1]):
                if self[seqid][position] == 0:
                    count0 = count0 + 1
                elif self[seqid][position] == 1:
                    count1 = count1 + 1
            count_list.append([featureID,count0,count1])
        return count_list
    

    def at_content(self, genome_sequence):
        for seqid in self:
            for position in range(len(self[seqid])):
                if genome_sequence[seqid][position] in "ATat":
                    self[seqid][position] = 1

    def sliding_window_calculate(self, window_size, window_jump = 1, operation = "sum", output = "dict",
                                 threshold = 1, seqs_to_exclude = []):
        """operation may be "sum", "average", or "set". Output may be "dict","coords", or "annotation_set"."""
        if output == "dict":
            new_dic = {}
        elif output == "annotation_set":
            annotation_set = AnnotationSet()
            annotation_set.region = {}
        elif output == "coords":
            coords_list = []
        for seqid in self:
            if len(self[seqid]) > window_size and not seqid in seqs_to_exclude:
                if output == "dict":
                    new_dic[seqid] = []
                window_in = False
                coords_list = []
                for position in range(len(self[seqid][:-1 * window_size]) / window_jump):
                    window_start = position * window_jump
                    if operation == "set":
                        new_dic[seqid].append(set(self[seqid][window_start:window_start + window_size]))
                        continue
                    window_sum = numpy.sum(self[seqid][window_start:window_start + window_size])
                    if operation == "sum":
                        value = window_sum
                    elif operation == "average":
                        value = window_sum * 1.0 / window_size
                    if output == "dict":
                        new_dic[seqid].append(value)
                    elif output == 'annotation_set':
                        if value >= threshold:
                            if window_in == False:
                                window_in = True
                                if len(coords_list) == 0:
                                    coords_list.append([1 + window_start])
                                elif position > coords_list[-1][1]:
                                    coords_list.append([1 + window_start])
                        else:
                            if window_in == True:
                                window_in = False
                                if len(coords_list[-1]) == 1:
                                    coords_list[-1].append(window_start + window_size)
                                else:
                                    coords_list[-1][1] = window_start + window_size
                    elif output == "coords":
                        if threshold[0] <= value <= thredshold[1]:
                            coords_list.append([seqid, 1 + window_start, window_start + window_size])
                    elif type(output) == file:
                        output.write(seqid + '\t' + str(value) + '\n')
                    if window_start % 10000000 == 0 and verbose:
                        print "processed " + seqid + " to position " + str(position)
                if output == "annotation_set":
                    if len(coords_list) > 0:
                        if len(coords_list[-1]) == 1:
                            coords_list[-1].append(len(self[seqid]))
                        for coords in coords_list:
                            ID = seqid + "-window" + str(coords[0])
                            annotation_set.region[ID] = BaseAnnotation(ID, seqid, tuple(coords), "region", annotation_set = annotation_set)
                if verbose:
                    print "processed " + seqid
                    if output == 'annotation_set':
                        for region in annotation_set.region:
                            if annotation_set.region[region].seqid == seqid:
                                print str(annotation_set.region[region].coords)                    
        if output == "dict":
            return new_dic
        elif output == "annotation_set":
            return copy.deepcopy(annotation_set)
        elif output == "coords":
            return coords_list





