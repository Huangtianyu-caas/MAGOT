#!/usr/bin/python
#GenomePy: a functional and simple library for genomic analysis
import copy


def vulgar2gff(vulgarstring, feature_types=['match','match_part'],source='exonerate'):
    """takes vulgar alignment string and outputs gff lines. For eventual use with a read_exonerate function"""
    #sets variables to be added to gff lines
    vulgarlist = vulgarstring.split()
    qname = vulgarlist[0]
    qstart = vulgarlist[1]
    qend = vulgarlist[2]
    qstrand = vulgarlist[3]
    tname = vulgarlist[4]
    tstart = vulgarlist[5]
    tend = vulgarlist[6]
    tstrand = vulgarlist[7]
    score = vulgarlist[8]
    vulgartrips = vulgarlist[9:]
    addfeat = False
    IDnum = 1
    if tstrand == "+":
        tposition = int(tstart) + 1
    else:
        tposition = int(tstart)
        tend = str(int(tend)+1)
    #makes top level gff line for hit
    gfflines = ["\t".join([tname,source,feature_types[0],str(tposition),tend,score,tstrand,'.','ID='+qname])]
    #iterates over matches within hit to make bottom level feature gff lines
    for i in range(len(vulgartrips)):
        field = vulgartrips[i]
        if i % 3 == 0:
            if field == 'M' or field == 'S' or field == 'G' or field == 'F':
                if not addfeat:
                    addfeat = True
                    line_to_add = [tname,source,feature_types[1]]
                    coords = [str(tposition)]
                else:
                    pass
            elif addfeat:
                gfflines.append('\t'.join(line_to_add+[str(min(coords)),str(max(coords)),'.',tstrand,'.',
                                                       'ID='+qname+'_'+feature_types[1]+str(IDnum)+';Parent='+qname]))
                IDnum = IDnum + 1
                addfeat = False
            else:
                pass
        if i % 3 == 2:
            if tstrand == "+":
                tposition = tposition + int(field)
            elif tstrand == "-":
                tposition = tposition - int(field)
            if addfeat == True:
                if tstrand == '+':
                    coords.append(str(tposition-1))
                elif tstrand == '-':
                    coords.append(str(tposition+1))
    if addfeat:
        gfflines.append('\t'.join(line_to_add+[str(min(coords)),str(max(coords)),'.',tstrand,'.',
                                                       'ID='+qname+'_'+feature_types[1]+str(IDnum)+';Parent='+qname]))
    return '\n'.join(gfflines)

    
def read_exonerate(exonerate_output,annotation_set_to_modify = None):
    if annotation_set_to_modify == None:
        annotation_set = AnnotationSet()
    else:
        annotation_set = annotation_set_to_modify
    exonerate_lines = read_to_string(exonerate_output).split('\n')
    gfflines = []
    for line in exonerate_lines:
        if line[:8] == "vulgar: ":
            gfflines.append(vulgar2gff(line[8:]))
    read_gff3("\n".join(gfflines), annotation_set_to_modify = annotation_set)
    if annotation_set_to_modify == None:
        return annotation_set
    
    
def write_longform_gff(annotation_set,keep_UTR_features = False):
    """returns gff string formated for compatability with Apollo genome annotation """
    gfflines = []
    fields = ['seqid','source','feature_type','get_coords()[0]','get_coords()[1]','score','strand','phase']
    #adds matches to gff
    if 'match' in annotation_set.__dict__:
        for match in annotation_set.match:
            match_obj = annotation_set.match[match]
            newline_list = []
            for field in fields:
                try:
                    newline_list.append(str(eval('match_obj.' + field)))
                except:
                    newline_list.append('.')
            attribute_list = ['ID=' + match_obj.ID]
            for attribute in match_obj.__dict__:
                if not attribute in fields+['annotation_set','parent','child_list','ID']:
                    attribute_list.append(attribute + '=' + eval('match_obj.' + attribute))
            newline_list.append(';'.join(attribute_list))
            gfflines.append('\t'.join(newline_list))
            for match_part in match_obj.child_list:
                match_part_obj = annotation_set[match_part]
                newline_list = []
                for field in fields:
                    try:
                        newline_list.append(str(eval('match_part_obj.' + field)))
                    except:
                        newline_list.append('.')
                attribute_list = ['ID=' + match_part_obj.ID,'Parent=' + match_part_obj.parent]
                for attribute in match_part_obj.__dict__:
                    if not attribute in fields+['annotation_set','parent','child_list','ID','coords']:
                        attribute_list.append(attribute + '=' + eval('match_part_obj.' + attribute))
                newline_list.append(';'.join(attribute_list))
                gfflines.append('\t'.join(newline_list))
        #adds genes to gff
    if 'gene' in annotation_set.__dict__:
        for gene in annotation_set.gene:
            gene_obj = annotation_set.gene[gene]
            newline_list = []
            for field in fields:
                try:
                    newline_list.append(str(eval('gene_obj.' + field)))
                except:
                    newline_list.append('.')
            attribute_list = ['ID='+gene_obj.ID]
            for attribute in gene_obj.__dict__:
                if not attribute in fields+['annotation_set','parent','child_list','ID']:
                    attribute_list.append(attribute + '=' + eval('gene_obj.' + attribute))
            newline_list.append(';'.join(attribute_list))
            gfflines.append('\t'.join(newline_list))
            for gene_child in gene_obj.child_list:
                gene_child_obj = annotation_set[gene_child]
                if gene_child_obj.feature_type == 'transcript':
                    transcript_obj = gene_child_obj
                    newline_list = []
                    for field in fields:
                        try:
                            newline_list.append(str(eval('transcript_obj.'+field)).replace('transcript','mRNA'))
                        except:
                            newline_list.append('.')
                    attribute_list = ['ID='+transcript_obj.ID,'Parent='+transcript_obj.parent]
                    for attribute in gene_obj.__dict__:
                        if not attribute in fields+['annotation_set','parent','child_list','ID']:
                            attribute_list.append(attribute + '=' + eval('transcript_obj.' + attribute))
                    newline_list.append(';'.join(attribute_list))
                    gfflines.append('\t'.join(newline_list))
                    exondict = {}
                    CDS_UTR_dict = {}
                    for transcript_child in transcript_obj.child_list:
                        transcript_child_obj = annotation_set[transcript_child]
                        line_base_list = []
                        for field in fields:
                            try:
                                line_base_list.append(str(eval('transcript_child_obj.'+ field)))
                            except:
                                line_base_list.append('.')
                        exon_attributes = 'ID=' + transcript_child_obj.ID + '-exon;Parent=' + transcript_child_obj.parent
                        transcript_child_attribute_list = ['ID=' + transcript_child_obj.ID, 'Parent=' + transcript_child_obj.parent]
                        for attribute in transcript_child_obj.__dict__:
                            if not attribute in fields+['annotation_set','parent','child_list','ID','coords']:
                                transcript_child_attribute_list.append(attribute + '=' + eval('transcript_child_obj.' + attribute))
                        transcript_child_attributes = ';'.join(transcript_child_attribute_list)
                        exondict[transcript_child_obj.coords] = '\t'.join(line_base_list).replace('CDS','exon').replace('UTR','exon') + '\t' + exon_attributes
                        CDS_UTR_dict[transcript_child_obj.coords] = '\t'.join(line_base_list) +'\t' + transcript_child_attributes
                    exondict_list = list(exondict)
                    exondict_list.sort()
                    CDS_UTR_dict_list = list(CDS_UTR_dict)
                    CDS_UTR_dict_list.sort()
                    #merges adjacent exons from abbutting UTRs and CDSs and writes exon gff lines
                    for i in range(len(exondict_list) - 1):
                        if exondict_list[i][1] + 1 == exondict_list[i+1][0]:
                            del exondict[exondict_list[i]]
                            new_exon_list = exondict[exondict_list[i+1]].split('\t')
                            exondict[exondict_list[i+1]] = '\t'.join(new_exon_list[:3]+[str(exondict_list[i][0]),str(exondict_list[i+1][1])]+new_exon_list[5:])
                        else:
                            gfflines.append(exondict[exondict_list[i]])
                    gfflines.append(exondict[exondict_list[-1]])
                    for i in CDS_UTR_dict_list:
                        if CDS_UTR_dict[i].split('\t')[2] == 'CDS' or keep_UTR_features:
                            gfflines.append(CDS_UTR_dict[i])
                else:
                    print 'ERROR: currently only accepts AnnotationSets with gene format CDS/UTR -> transcript -> gene'
                    break
    return '\n'.join(gfflines)
    

def read_to_string(potential_file):
    """Tries to read "potential_file" into string- first as file location,
    then file, then string."""
    try:
        output_string = open(potential_file).read().replace('\r','')
    except IOError:
        try:
            output_string = potential_file.read().replace('\r','')
        except AttributeError:
            output_string = potential_file.replace('\r','')
    return output_string


def read_gff3(gff3,annotation_set_to_modify = None,gene_hierarchy = ['gene','mRNA',['CDS','five_prime_UTR','three_prime_UTR']],
              other_hierarchies = [['match','match_part']],features_to_ignore = ['exon'],
              features_to_replace = [('protein_match','match'),('expressed_sequence_match','match')]):
    #reads gff3 to string if supplied as a file or file location
    gff_list = read_to_string(gff3).split('\n')
    #reformates features_to_replace to be used in later command
    for feature in features_to_replace[:]:
        features_to_replace.append(str(feature))
        features_to_replace.remove(feature)
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
                feature_type = eval('gff_fields[2].replace' + '.replace'.join(features_to_replace))
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
                    if feature_type == 'CDS':
                        renamed_feature_type = 'CDS'
                    else:
                        renamed_feature_type = 'UTR'
                    eval('annotation_set.' + renamed_feature_type)[ID] = BaseAnnotation(ID,seqid,coords,renamed_feature_type,parent,
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
                                                                                  annotation_set = annotation_set, other_attributes = copy.copy(other_attributes))
                elif feature_type in features_to_ignore:
                    pass
                else:
                    in_hierarchy = False
                    for other_hierarchy in other_hierarchies:
                        if feature_type == other_hierarchy[-1]:
                            in_hierarchy = True
                            coords_list = [int(gff_fields[3]),int(gff_fields[4])]
                            coords_list.sort()
                            coords=tuple(coords_list)
                            create_parents_chain = other_hierarchy[:-1]
                            create_parents_chain.reverse()
                            eval('annotation_set.' + feature_type)[ID] = BaseAnnotation(ID, seqid, coords, feature_type, parent,
                                                                              other_attributes = copy.copy(other_attributes),
                                                                              annotation_set = annotation_set,
                                                                              create_parents_chain = create_parents_chain)
                        elif feature_type in other_hierarchy:
                            in_hierarchy = True
                            eval('annotation_set.' + feature_type)[ID] = ParentAnnotation(ID, seqid, feature_type, parent = parent,
                                                                                          annotation_set = annotation_set, other_attributes = copy.copy(other_attributes))
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


def read_blast_csv(blast_csv,annotation_set_to_modify = None,hierarchy = ['match','match_part'], source = 'blast'):
    """Reads csv output from blast (-outfmt 10) into an AnnotationSet object. Currently does not string hits together because I'm
    biased towards working on genes in tandem arrays where stringing hits together is annoying. May add option in future."""
    #reads blast_csv from file location, file, or string
    blast_lines = read_to_string(blast_csv).split('\n')
    #checks if annotation_set is given and creates annotation_set if not
    if annotation_set_to_modify == None:
        annotation_set = AnnotationSet()
    else:
        annotation_set = annotation_set_to_modify
    id_generator_dict = {}
    feature_type = hierarchy[-1]
    create_parents_chain = hierarchy[:-1]
    create_parents_chain.reverse()
    if not feature_type in annotation_set.__dict__:
        setattr(annotation_set, feature_type, {})
    for line in blast_lines:
        fields = line.split(',')
        seqid = fields[1]
        tstart = int(fields[8])
        tend = int(fields[9])
        if tstart < tend:
            coords = (tstart,tend)
            strand = '+'
        else:
            coords = (tend,tstart)
            strand = '-'
        score = fields[11]
        IDbase = fields[0]
        if IDbase in eval('annotation_set.' + feature_type):
            ID = IDbase + '-' + str(id_generator_dict[IDbase])
            id_generator_dict[IDbase] = id_generator_dict[IDbase] + 1
            while ID in eval('annotation_set.' + feature_type):
                ID = IDbase + '-' + str(id_generator_dict[IDbase])
                id_generator_dict[IDbase] = id_generator_dict[IDbase] + 1
        else:
            ID = IDbase
            id_generator_dict[IDbase] = 1
        other_attributes = {}
        other_attributes['evalue'] = fields[10]
        parent = ID + '-match'
        eval('annotation_set.' + feature_type)[ID] = BaseAnnotation(ID, seqid, coords, feature_type, parent, other_attributes,
                                                                    annotation_set, create_parents_chain)
    if annotation_set_to_modify == None:
        return annotation_set


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
            if type(eval("self." + attribute)) == dict:
                try:
                    all_dicts[item] = eval("self." + attribute)[item]
                except:
                    pass
        return all_dicts[item]
    
    def read_gff3(self, gff3):
        read_gff3(gff3, annotation_set_to_modify = self)
    
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
                if len(create_parents_chain) > 1:
                    self.parent = parent + '-' + create_parents_chain[0]
                else:
                    self.parent = parent
                parent_to_create = parent + '-' + hierarchy[active_feature_type]
                #checks if parent needs to be created and creates parent
                if len(create_parents_chain) > 1:
                    while not parent in annotation_set.__dict__[hierarchy[active_feature_type]] and not parent_to_create in annotation_set.__dict__[hierarchy[active_feature_type]]:
                        parent_to_create = parent + '-' + hierarchy[active_feature_type]
                        #checks whether this is second-to-last round of parent creation
                        if active_feature_type == create_parents_chain[-2]:
                            parent_to_create = parent
                            break
                        elif parent in annotation_set.__dict__[hierarchy[hierarchy[active_feature_type]]] or hierarchy[active_feature_type] == create_parents_chain[-2]:
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
                except KeyError:
                    try:
                        annotation_set.__dict__[hierarchy[active_feature_type]][parent_to_create].child_list.append(ID)
                    except KeyError:
                        annotation_set.__dict__[hierarchy[active_feature_type]][parent] = ParentAnnotation(ID = parent, seqid = seqid,
                                                                                                       feature_type = hierarchy[active_feature_type],
                                                                                                       child_list = [active_feature_ID],parent = None,
                                                                                                       annotation_set = annotation_set)
    def get_coords(self):
        return self.coords
    


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
                block = locus.split('\n')
                seqid = block[0]
                seqstring = Sequence("".join(block[1:]))
                self[seqid] = seqstring


class Genome():
    """genome class, which contains sequence and annotations"""
    def __init__(self,genome_sequence = None, annotations = None):
        if genome_sequence.__class__.__name__ == 'GenomeSequence' or genome_sequence == None:
            self.genome_sequence = genome_sequence
        else:
            self.genome_sequence = GenomeSequence(genome_sequence)
        self.annotations = annotations
    
    def get_scaffold_fasta(self, seqid):
        return '>' + seqid + '\n' + self.genome_sequence[seqid]
    
    def write_apollo_gff(self, seqid):
        if self.genome_sequence != None and self.annotations != None:
            return write_to_longform_gff(self.annotations.get_seqid(seqid)) + '\n' + self.get_scaffold_fasta(seqid)
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
