#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas

#This library will eventually house interesting functions for annotations- e.g. set opperations, well, not sure what else


import copy
import StringIO
from genome import *

#Very experimental, eventually should be developed into a full set operation function (e.g. unions, differences, etc.)
def annotation_overlap(annotation_set_features1, annotation_set_features2):
    annotation_coords_dict = {}
    contained_list = []
    for annotation in annotation_set_features1:
        annotation_obj = annotation_set_features1[annotation]
        seqid = annotation_obj.seqid
        if seqid in annotation_coords_dict:
            annotation_coords_dict[seqid].append(annotation_obj.get_coords())
        else:
            annotation_coords_dict[seqid] = [annotation_obj.get_coords()]
    for annotation in annotation_set_features2:
        annotation_obj = annotation_set_features2[annotation]
        seqid = annotation_obj.seqid
        coords2 = annotation_obj.get_coords()
        if seqid in annotation_coords_dict:
            for coords1 in annotation_coords_dict[seqid]:
                if coords1[0] <= coords2[0] <= coords1[1] or coords1[0] <= coords2[1] <= coords1[1]:
                    contained_list.append(annotation_obj.ID)
    return contained_list