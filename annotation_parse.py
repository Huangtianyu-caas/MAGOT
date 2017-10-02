#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas


import copy
import StringIO
from genome import *

def vulgar2gff(vulgarlist, feature_types=['match','match_part'],source='exonerate'):
    """takes vulgar alignment list (e.g. vulgarstring.split() ) and outputs gff lines. For eventual use with a read_exonerate function"""
    #sets variables to be added to gff lines
    qname = vulgarlist[0] + '-against-' + vulgarlist[4]
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
    exonerate_lines = ensure_file(exonerate_output)
    gfflines = []
    IDdic = {}
    qname = ""
    tname = ""
    for original_line in exonerate_lines:
        line = original_line.replace('\r','').replace('\n','')
        if line[:16] == "         Query: ":
            qname = line[16:]
        elif line[:16] == "        Target: ":
            tname = line[16:].replace(':[revcomp]','').replace('[revcomp]','')
            if tname[-1] == " ":
                tname = tname[:-1]
        elif line[:8] == "vulgar: ":
            vulgar_line_list = line[8:].split()
            #trying to makesure IDs are unique
            vulgar_line_list[0] = qname
            vulgar_line_list[4] = tname
            ID = vulgar_line_list[0] + '-against-' + vulgar_line_list[4]
            if ID in IDdic:
                vulgar_line_list[0] = vulgar_line_list[0] + str(IDdic[ID])
                IDdic[ID] = IDdic[ID] + 1
            else:
                IDdic[ID] = 1
            gfflines.append(vulgar2gff(vulgar_line_list))
    read_gff("\n".join(gfflines), annotation_set_to_modify = annotation_set)
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

def write_gff(annotation_set, gff_format = "simple gff3"):
    gff_lines = []
    for attribute in annotation_set.__dict__:
        if type(annotation_set.__dict__[attribute]) == dict:
            if len(annotation_set.__dict__[attribute]) > 0:
                if annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].__class__.__name__ == "ParentAnnotation" or annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].__class__.__name__ == "BaseAnnotation":
                    if annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].parent == None:
                        for annotationID in annotation_set.__dict__[attribute]:
                            annotation = annotation_set.__dict__[attribute][annotationID]
                            gff_lines.append(annotation.get_gff(gff_format))
    return '\n'.join(gff_lines)
                        
    

def read_gff(gff,annotation_set_to_modify = None, base_features = ['CDS','match_part','similarity','region'],features_to_ignore = ['exon'],
    gff_version = "auto", parents_hierarchy = [], features_to_replace = [], IDfield = "ID", parent_field = "Parent", presets = None):
    """Assumptions of file:    
    Relatively consistent, nothing crazy like mixed gff3 and gff2
    
    If a parent's children are base features, they should NOT overlap. E.g. if your file has exons and CDS, only
    give CDS as base_features and give exons as features_to_ignore.
    
    Although parent feature types may be used multiply in hierarchies, any given base feature
    will only be present in one hierarchy
    
    If parent features are not represented by lines (e.g. minimal gff2 format or cegma), they
    can be inferred from defline attributes given in parents_hierarchy. It is expected that the order they
    are given is consistent with the hierarchy (e.g. ['mRNA','gene']). If an underscore is used in
    this defline, the word before the underscore is used as the parent feature type. The "Parent"
    defline attribute is reserved for the imediate parent of the feature- if not given, this will be assumed to
    be the lowest-level parent specified in defline.
    """
    replace_dict = {"\n":"","\r":""}
    presets_dict = {}
    presets_dict["augustus"] = "features_to_ignore = ['gene','transcript','stop_codon','terminal','internal','initial','intron',\
        'start_codon','single']\nparent_field = None\nparents_hierarchy = ['transcript_id','gene_id']\nIDfield = None"
    presets_dict["RepeatMasker"] = "parent_field = None\nIDfield = 'Target'"
    presets_dict["CEGMA"] = "parent_field = ''\nIDfield = None\nparents_hierarchy = ['gene_id\nfeatures_to_replace = \
        [['First','CDS']['Internal','CDS']['Terminal','CDS']['Single','CDS']]\ngff_version = 3"
    if presets in presets_dict:
        exec(presets_dict[presets])
    version = gff_version
    gff_file = ensure_file(gff)
    for feature in features_to_replace:
        replace_dict["\t" + feature[0] + "\t"] = "\t" + feature[1] + "\t"
    if annotation_set_to_modify == None:
        annotation_set = AnnotationSet()
    else:
        annotation_set = annotation_set_to_modify
    #this dictionary helps generate names if features passed with identical ID fields
    generate_new_ID_dict = {}
    #figures out which feature types are base annotations
    base_dict = {}
    #here we go folks
    for original_line in gff_file:
        if original_line[0] != "#" and original_line.count('\t') == 8:
            line = original_line
            for string in replace_dict:
                line = line.replace(string,replace_dict[string])
            fields = line.split('\t')
            if version == "auto":
                if "=" in fields[8]:
                    version = 3
                else:
                    version = 2
            ID = None
            parent = None
            other_attributes = {}
            seqid = fields[0]
            other_attributes['source'] = fields[1]
            feature_type = fields[2]
            if feature_type in features_to_ignore:
                continue
            coords = [int(fields[3]),int(fields[4])]
            coords.sort()
            coords = tuple(coords)
            try:
                other_attributes['score'] = float(fields[5])
            except ValueError:
                pass
            strand = fields[6]
            if fields[7] in ['0','1','2']:
                other_attributes['phase'] = int(fields[7])
            defline_dict = {}
            for defline_field in fields[8].split(';'):
                if defline_field != "":
                    if parent_field == "":
                        defline_dict[""] = defline_field
                    elif version == 2:
                        if '"' in defline_field:
                            defline_dict[defline_field.split()[0]] = defline_field.split('"')[1]
                        else:
                            try: defline_dict[defline_field.split()[0]] = defline_field.split()[1]
                            except:
                                print defline_field
                                return None
                    elif version == 3:
                        defline_dict[defline_field.split('=')[0]] = defline_field.split('=')[1]
            #Tries to figure out parent
            if parent_field != None:
                if parent_field in defline_dict:
                    parent = defline_dict[parent_field]
            elif parents_hierarchy != []:
                for parent_type in parents_hierarchy:
                    if parent_type in defline_dict:
                        parent = defline_dict[parent_type]
                        break
            #Assigns ID
            if IDfield != None:
                try:
                    ID = defline_dict[IDfield]
                except KeyError:
                    if parent != None:
                        ID = parent + '-' + feature_type
            elif parent != None:
                ID = parent + '-' + feature_type
            else:
                ID = seqid + '-' + feature_type + fields[3]
            #checks if ID already in annotation_set and adjusts ID if necessary
            try:
                annotation_set[ID]
                if ID in generate_new_ID_dict:
                    generate_new_ID_dict[ID] = generate_new_ID_dict[ID] + 1
                    ID = ID + "-" + str(generate_new_ID_dict[ID])
                else:
                    generate_new_ID_dict[ID] = 2
                    ID = ID + '2'
            except KeyError:
                pass
            #Creates parents if necessary, adds to parent if exists
            if parent != None:
                child_to_assign = ID
                for parent_feature_index in range(len(parents_hierarchy)):
                    parent_feature = parents_hierarchy[parent_feature_index]
                    if parent_feature in defline_dict:
                        parent_feature_ID = defline_dict[parent_feature]
                        parent_feature_type = parent_feature.split('_')[0]
                        parents_parent = None
                        if parent_feature_index != len(parents_hierarchy) - 1:
                            for parents_parent_feature in parents_hierarchy[parent_feature_index + 1:]:
                                if parents_parent_feature in defline_dict:
                                    parents_parent = defline_dict[parents_parent_feature]                        
                        if not parent_feature_type in annotation_set.__dict__:
                            annotation_set.__dict__[parent_feature_type] = {}
                        if parent_feature_ID in annotation_set.__dict__[parent_feature_type]:
                            if not child_to_assign in annotation_set.__dict__[parent_feature_type][parent_feature_ID].child_list:
                                annotation_set.__dict__[parent_feature_type][parent_feature_ID].child_list.append(child_to_assign)
                        else:
                            annotation_set.__dict__[parent_feature_type][parent_feature_ID] = ParentAnnotation(parent_feature_ID, seqid,
                                                                                            parent_feature_type, child_list = [child_to_assign],
                                                                                            parent = parents_parent, strand = strand,
                                                                                            annotation_set = annotation_set)
                        child_to_assign = parent_feature_ID
                #In case of no parent from parent hierarchy in defline_dict
                try:
                    if not ID in annotation_set[parent].child_list:
                        annotation_set[parent].child_list.append(ID)
                except KeyError:
                    print """It seems that this line has a parent attribute but that that parent doesn't have a line itself nor
                    does this line have a defline attribute that specifies a parent type. I'm afraid this function can't currently
                    deal with that."""
                    print ID
                    print parent
                    return None
            #fills other_attributes from defline
            for defline_attribute in defline_dict:
                if not defline_attribute in [IDfield,parent_field]:
                    other_attributes[defline_attribute] = defline_dict[defline_attribute]
            #And now to create the feature!
            if not feature_type in annotation_set.__dict__:
                annotation_set.__dict__[feature_type] = {}
            if feature_type in base_features:
                annotation_set.__dict__[feature_type][ID] = BaseAnnotation( ID, seqid, coords, feature_type, parent,
                                                                           strand , other_attributes, annotation_set)
            else:
                child_list = []
                annotation_set.__dict__[feature_type][ID] = ParentAnnotation(ID, seqid, feature_type, child_list,
                                                                             parent, strand, annotation_set, other_attributes)
    if annotation_set_to_modify == None:
        return copy.deepcopy(annotation_set)


def read_cegma_gff(cegma_gff,annotation_set_to_modify = None):
    """reads gff produced by CEGMA and returns AnnotationSet populated by CEGMA predictions"""
    annotation_set = read_gff(cegma_gff,annotation_set_to_modify = annotation_set_to_modify, presets = "CEGMA")
    if annotation_set_to_modify == None:
        return annotation_set


def read_blast_csv(blast_csv,annotation_set_to_modify = None,hierarchy = ['match','match_part'], source = 'blast', find_truncated_locname = False):
    """Reads csv output from blast (-outfmt 10) into an AnnotationSet object. Currently does not string hits together because I'm
    biased towards working on genes in tandem arrays where stringing hits together is annoying. May add option in future."""
    #reads blast_csv from file location, file, or string
    blast_file = ensure_file(blast_csv)
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
    if find_truncated_locname:
        if annotation_set.genome == None:
            print '"warning: find_truncated_locname" was set to true, but annotation set has no associated genome object so this cannot be done'
            find_truncated_locname = False
        else:
            genome_seqids = annotation_set.genome.get_seqids()
    for whole_line in blast_file:
        line = whole_line.replace('\r','').replace('\n','')
        fields = line.split(',')
        if len(fields) > 8:
            seqid = fields[1]
            if find_truncated_locname:
                if not seqid in genome_seqids:
                    for genome_seqid in genome_seqids:
                        if seqid == genome_seqid.split()[0]:
                            seqid = genome_seqid
                            break
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
            other_attributes['score'] = score
            parent = ID + '-' + create_parents_chain[0]
            child_to_set = ID
            for parent_index in range(len(create_parents_chain)):
                parent_feature = create_parents_chain[parent_index]
                if not parent_feature in annotation_set.__dict__:
                    setattr(annotation_set, parent_feature, {})
                if parent_index != len(create_parents_chain) - 1:
                    parent_to_set = ID + '-' + create_parents_chain[parent_index + 1]
                else:
                    parent_to_set = None
                annotation_set.__dict__[parent_feature][ID + '-' + parent_feature] = ParentAnnotation(ID + '-' + parent_feature, seqid, parent_feature,
                                                                                                    [child_to_set], parent_to_set, strand,
                                                                                                    annotation_set, other_attributes = {})
                child_to_set = ID + '-' + parent_feature
            eval('annotation_set.' + feature_type)[ID] = BaseAnnotation(ID, seqid, coords, feature_type, parent, strand,
                                                                        other_attributes, annotation_set)

            
    if annotation_set_to_modify == None:
        return annotation_set