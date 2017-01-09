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


def read_gff3(gff3,annotation_set_to_modify = None,gene_hierarchy = ['gene','mRNA',['CDS','five_prime_UTR','three_prime_UTR']],
              other_hierarchies = [['match','match_part']],features_to_ignore = ['exon']):
    #reads gff3 to string if supplied as a file or file location
    gff_list = read_to_string(gff3).split('\n')
    #checks if annotation_set is given and creates annotation_set if not
    if annotation_set_to_modify == None:
        annotation_set = AnnotationSet()
    else:
        annotation_set = annotation_set_to_modify
    #this dictionary helps generate names if features are passed with parents but not IDs
    generate_ID_from_parent_dict = {}
    #this dictionary helps generate names if features passed with identical ID fields
    generate_new_ID_dict = {}
    #Fills annotiation_set
    for gff_line in gff_list:
        if len(gff_line) > 1:
            if gff_line[0] != '#' and gff_line.count('\t') == 8:
                try:
                    del parent
                except:
                    pass
                try:
                    del ID
                except:
                    pass
                gff_fields = gff_line.split('\t')
                other_attributes = {}
                seqid = gff_fields[0]
                other_attributes['source'] = gff_fields[1]
                feature_type = gff_fields[2]
                other_attributes['score'] = gff_fields[5]
                other_attributes['strand'] = gff_fields[6]
                for additional_attribute in gff_fields[8].split(';'):
                    if '=' in additional_attribute:
                        attr_split = additional_attribute.split('=')
                        if attr_split[0] == 'ID':
                            ID = attr_split[1]
                        elif attr_split[0] == 'Parent':
                            parent = attr_split[1]
                        else:
                            other_attributes[attr_split[0]] = attr_split[1]
                #checks for ID and/or parent in attribute field
                try:
                    ID
                except NameError:
                    try:
                        ID_base = parent + '-' + feature_type
                        if not ID_base in generate_ID_from_parent_dict:
                            generate_ID_from_parent_dict[ID_base] = 0
                        ID = parent + '-' + feature_type + str(generate_ID_from_parent_dict[ID_base])
                        generate_ID_from_parent_dict[ID_base] = generate_ID_from_parent_dict[ID_base] + 1
                    except NameError:
                        print "This gff (or at least one feature therein) seems to have an attributes field without\
                        'ID' or 'Parent' varaibles. This is not yet supported."
                        return
                #checks if feature with same ID already exists, generates a new ID if so. This should only happen with
                #base level features in dumb gff3s, otherwise something is wrong.
                try:
                    annotation_set[ID]
                    if ID in generate_new_ID_dict:
                        generate_new_ID_dict[ID] = generate_new_ID_dict[ID] + 1
                    else:
                        generate_new_ID_dict[ID] = 1
                    ID = ID + '_' + str(generate_new_ID_dict[ID])
                except:
                    pass
                #checks if feature_type in annotation_set, adds if not unless in features_to_ignore
                if not feature_type in annotation_set.__dict__ and not feature_type in features_to_ignore:
                    setattr(annotation_set, feature_type, {})
                #sets parent to None if does not exist
                try:
                    parent
                except NameError:
                    parent = None
                #creates annotations
                if feature_type in gene_hierarchy[-1] :
                    coords_list = [int(gff_fields[3]),int(gff_fields[4])]
                    coords_list.sort()
                    coords=tuple(coords_list)
                    eval('annotation_set.' + feature_type)[ID] = BaseAnnotation(ID,seqid,coords,feature_type,parent,
                                                                              other_attributes = copy.copy(other_attributes),
                                                                              annotation_set = annotation_set)
                elif feature_type in gene_hierarchy:
                    make_feature = True
                    renamed_feature_type = feature_type
                    if len(gene_hierarchy) > 2:
                        if feature_type == gene_hierarchy[1]:
                            renamed_feature_type = 'transcript'
                    if make_feature:
                        eval('annotation_set.' + renamed_feature_type)[ID] = ParentAnnotation(ID, seqid, renamed_feature_type,parent = parent,
                                                                                  annotation_set = annotation_set)
                elif feature_type in features_to_ignore:
                    pass
                else:
                    in_hierarchy = False
                    for other_hierarchy in other_hierarchies:
                        if feature_type in other_hierarchy[-1]:
                            in_hierarchy = True
                            coords_list = [int(gff_fields[3]),int(gff_fields[4])]
                            coords_list.sort()
                            coords=tuple(coords_list)
                            create_parent_chain = other_hierarchy[:-1]
                            create_parent_chain.reverse()
                            eval('annotation_set.' + feature_type)[ID] = BaseAnnotation(ID, seqid, coords, feature_type, parent,
                                                                              other_attributes = copy.copy(other_attributes),
                                                                              annotation_set = annotation_set,
                                                                              create_parent_chain = create_parent_chain)
                        elif feature_type in other_hierarchy:
                            in_hierarchy = True
                            eval('annotation_set.' + feature_type)[ID] = ParentAnnotation(ID, seqid, feature_type, parent = parent,
                                                                                          annotation_set = annotation_set)
                    if not in_hierarchy:
                        coords_list = [int(gff_fields[3]),int(gff_fields[4])]
                        coords_list.sort()
                        coords=tuple(coords_list)
                        eval('annotation_set.' + feature_type)[ID] = BaseAnnotation(ID, seqid, coords, feature_type, parent,
                                                                              other_attributes = copy.copy(other_attributes),
                                                                              annotation_set = annotation_set,
                                                                              create_parents_chain = None)
    if annotation_set_to_modify == None:
        return(annotation_set)


def read_cegma_gff(cegma_gff,annotation_set_to_modify = None):
    """reads gff produced by CEGMA and returns AnnotationSet populated by CEGMA predictions"""
    modified_gff = read_to_string(cegma_gff).replace('\tFirst\t','\tCDS\t').replace('\tInternal\t','\tCDS\t').replace('\tTerminal\t','\tCDS\t').replace('KOG','Parent=KOG')
    annotation_set = read_gff3(modified_gff,annotation_set_to_modify = annotation_set_to_modify, gene_hierarchy=['CDS'])
    if annotation_set_to_modify == None:
        return annotation_set


def write2apollo(annotation_set):
    for gene in annotation_set[gene]:
        pass
    



class AnnotationSet():
    """A set of annotations of a single genome. Each feature type (e.g. gene, transcript, exon, etc.)
    is stored in it's own dictionary as Annotations with their ID as their key (see "Annotation" class).
    The AnnotationSet itself also functions losely as a dictionary, in that any feature can be returned
    by indexing the AnnotationSet with the ID as a key (e.g. my_annotation_set["my_feature_ID"])"""
    def __init__(self):
        self.gene = {}
        self.transcript = {}
        self.CDS = {}
        self.UTR = {}
    def __getitem__(self,item):
        all_dicts = {}
        for attribute in dir(self):
            if type(eval("self."+attribute)) == dict:
                try:
                    all_dicts[item] = eval("self."+attribute)[item]
                except:
                    pass
        return all_dicts[item]
    def read_gff3(self, gff3):
        read_gff3(gff3, annotation_set_to_modify = self)
    

class BaseAnnotation():
    """Bottom-most level annotation on a genome, for example CDS, UTR, Match, etc. Anything that should have no children"""
    def __init__(self, ID, seqid, coords, feature_type, parent = None, other_attributes = {}, annotation_set = None,
                 create_parents_chain=['transcript','gene']):
        #Sets up most attributes
        self.ID = ID
        self.seqid = seqid
        self.coords = coords
        self.feature_type = feature_type
        self.annotation_set = annotation_set
        for attribute in other_attributes:
            setattr(self, attribute, other_attributes[attribute])
        #checks if feature type needs to be added to annotation_set
        if annotation_set != None:
            if not feature_type in annotation_set.__dict__:
                setattr(annotation_set,feature_type, {})
            #checks if parent feature types need to be added to annotation_set
            if create_parents_chain != None:
                for parent_feature_type in create_parents_chain:
                    if not parent_feature_type in annotation_set.__dict__:
                        setattr(annotation_set,parent_feature_type,{})
        #checks if parent features present:
            if parent == None:
                pass
            #checks if parent needs to be created
            elif parent in annotation_set.__dict__[create_parents_chain[0]]:
                self.parent = parent
                annotation_set.__dict__[create_parents_chain[0]][parent].child_list.append(ID)
            #Executes parent creation process if parent needs to be created
            else:
                #sets up hierarchy for parent creation
                hierarchy = {feature_type: create_parents_chain[0]}
                for feature_index in range(len(create_parents_chain))[:-1]:
                    hierarchy[create_parents_chain[feature_index]] = create_parents_chain[feature_index + 1]
                active_feature_type = feature_type
                active_feature_ID = ID
                self.parent = parent + '-' + create_parents_chain[0]
                parent_to_create = parent + '-' + hierarchy[active_feature_type]
                #checks if parent needs to be created and creates parent
                if len(create_parents_chain) > 1:
                    while not parent in annotation_set.__dict__[hierarchy[active_feature_type]] and not parent_to_create in annotation_set.__dict__[hierarchy[active_feature_type]]:
                        parent_to_create = parent + '-' + hierarchy[active_feature_type]
                        #checks whether this is second-to-last round of parent creation
                        if active_feature_type == create_parents_chain[-2]:
                            next_level_parent = parent
                            break
                        elif parent in annotation_set.__dict__[hierarchy[hierarchy[active_feature_type]]]:
                            next_level_parent = parent
                        else:
                            next_level_parent = parent + '-' + hierarchy[hierarchy[active_feature_type]]
                        #creates parent
                        annotation_set.__dict__[hierarchy[active_feature_type]][parent_to_create] = ParentAnnotation(ID = parent_to_create,
                                                                                                                     seqid = seqid,
                                                                                                                     feature_type = hierarchy[active_feature_type],
                                                                                                                     child_list = [active_feature_ID],
                                                                                                                     parent = next_level_parent,
                                                                                                                     annotation_set = annotation_set)
                        #resets active feature type and ID
                        active_feature_type = hierarchy[active_feature_type]
                        active_feature_ID = parent_to_create
                #checks if parent was found or last-level parent needs to be created
                try:
                    annotation_set.__dict__[hierarchy[active_feature_type]][parent]
                    #debug
                    print ID
                    #
                except KeyError:
                    try:
                        annotation_set.__dict__[hierarchy[active_feature_type]][parent_to_create].child_list.append(ID)
                    except KeyError:
                        annotation_set.__dict__[hierarchy[active_feature_type]][parent] = ParentAnnotation(ID = parent, seqid = seqid,
                                                                                                       feature_type = hierarchy[active_feature_type],
                                                                                                       child_list = [active_feature_ID],parent = None,
                                                                                                       annotation_set = annotation_set)


class ParentAnnotation():
    """Parent of any BaseAnnotation. Examples include genes, transcripts, and exons. Suggested hierarchy for genes is
    CDS (as BaseAnnotation) -> exon -> transcript -> gene."""
    def __init__(self, ID, seqid, feature_type, child_list = [], parent = None, annotation_set = None, other_attributes = {}):
        self.ID = ID
        self.seqid = seqid
        self.feature_type = feature_type
        self.child_list = copy.copy(child_list)
        self.parent = parent
        self.annotation_set = annotation_set
        for attribute in other_attributes:
            setattr(self, attribute, other_attributes[attribute])
        if annotation_set != None:
            try:
                annotation_set[parent].child_list.append(ID)
            except:
                pass
    
    def get_coords(self):
        if len(child_list) > 0 and annotation_set != None:
            coords_list = []
            for child in self.child_list:
                child_object = self.annotation_set[child]
                if child_object.__class__ == ParentAnnotation:
                    coords_list = coords_list + list(child_object.get_coords())
                elif child_object.__class__ == BaseAnnotation:
                    coords_list = coords_list + list(child_object.coords)
                else:
                    print "for some reason you have children in ParentAnnotation " + self.ID + " which are neither \
                    ParentAnnotation objects nor BaseAnnotation object. Get your act together"
            return (min(coords_list),max(coords_list))


class Sequence(str):
    """DNA sequence. Will eventually have methods allowing reverse complimenting,
        translating, etc."""
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


class Genome():
    """genome class, which contains sequence and annotations"""
    def __init__(self,genome_sequence = None, annotations = None):
        self.genome_sequence = genome_sequence
        self.annotations = annotations
        
