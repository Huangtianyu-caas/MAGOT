#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas

#This whole library is neglected and may not work with well with the current MAGOT version.

#Theoretically, this library will deal with genomic variants (SNPs, indels, maybe eventually structural variants),
#performing tasks like heterozygosity calculations or genetic distance calculations.


import copy
import StringIO
from genome import *


def calculate_heterozygosity(variant_set,locus_list = None, window = 10000):
    """calculates heterozygosity at each locus specified in the vcf_header associated with the variant set,
    or by a locus_list. Loci are specified in the locus list as either "seqid,seq_length" for entire
    contigs/scaffold, or for regions within contigs/scaffold, "seqid,start_position,stop_position" """
    outlist = []
    if locus_list == None:
        #generates locus list from vcf header if no locus list passed
        if variant_set.vcf_headers != None:
            locus_list = []
            for vcf_header_line in variant_set.vcf_headers[0].split('\n'):
                if "##contig=<ID=" in vcf_header_line:
                    seqid = vcf_header_line.split(',')[0][13:]
                    seqlen = int(vcf_header_line.split(',')[1][7:-1])
                    seqchunk_number = seqlen / window
                    for seqchunk in range(seqchunk_number):
                        locus = seqid + ',' + str(1 + window * seqchunk) + ',' + str(window * (seqchunk + 1))
                        locus_list.append(locus)
                    locus = seqid + ',' + str(1 + window * seqchunk_number) + ',' + str(seqlen)
                    locus_list.append(locus)
        else:
            print "function run with out locus_list, but variant_set has no vcf_header. One day I'll add a \
            function to get loci from a GenomeSequence; if you really want this, pester me about it on github 'Issues'"
            return None
    locus_list = list(set(locus_list))
    print "built locus list"    
    print len(locus_list)
    locus_het_count_dic = {}
    loci_dic = {}
    for locus in locus_list:
            locus_het_count_dic[locus] = {}
            for genotype in variant_set[list(variant_set)[0]].genotypes:
                locus_het_count_dic[locus][genotype] = 0
            seqid = locus.split(',')[0]
            if seqid in loci_dic:
                loci_dic[seqid].append(locus)
            else:
                loci_dic[seqid] = [locus]
    for variant_code in variant_set:
        variant_code_split = variant_code.split('<:>')
        seqid = variant_code_split[0]
        position = int(variant_code_split[1])
        for locus in loci_dic[seqid]:
            locus_split = locus.split(',')
            if len(locus.split(',')) > 2:
                if int(locus_split[1]) <= position <= int(locus_split[2]):
                    for genotype in variant_set[variant_code].genotypes:
                        genotypeGT = variant_set[variant_code].genotypes[genotype].GT
                        if genotypeGT[0] != genotypeGT[-1]:
                            locus_het_count_dic[locus][genotype] = locus_het_count_dic[locus][genotype] + 1
            else:
                if seqid == locus_split[0]:
                    for genotype in variant_set[variant_code].genotypes:
                        genotypeGT = variant_set[variant_code].genotypes[genotype].GT
                        if genotypeGT[0] != genotypeGT[-1]:
                            locus_het_count_dic[locus][genotype] = locus_het_count_dic[locus][genotype] + 1
    standard_genotype_list = None
    for locus in locus_list:
        if standard_genotype_list == None:
            standard_genotype_list = list(locus_het_count_dic[locus])
            print standard_genotype_list
        locus_split = locus.split(',')
        if len(locus.split(',')) > 2:
            seqlen = int(locus_split[2]) - int(locus_split[1]) + 1
        else:
            seqlen = int(locus_split[1])
        if seqlen < 1:
            pass
        else:
            hetrate_list = []
            for genotype in standard_genotype_list:
                hetrate_list.append(str( 100.0 * locus_het_count_dic[locus][genotype] / seqlen))
            outlist.append(locus + "," + ",".join(hetrate_list))
    return '\n'.join(outlist)


def read_vcf(vcf, variant_set_to_modify = None):
    """reads vcf file into variant_set"""
    vcf_list = read_to_string(vcf).split('\n')
    vcf_headerlines = []
    if variant_set_to_modify == None:
        variant_set = VariantSet()
    else:
        variant_set = variant_set_to_modify
    if variant_set.vcf_headers == None:
        variant_set.vcf_headers = []
    vcf_header_index = len(variant_set.vcf_headers)
    for line in vcf_list:
        if line == "":
            pass
        elif line[:2] == "##":
            vcf_headerlines.append(line)
        elif line[0] == "#":
            field_names = line[1:].split('\t')
            if len(field_names) > 8:
                sample_names = field_names[9:]
        else:
            fields = line.split('\t')
            seqid = fields[0]
            position = fields[1]
            ref_allele = fields[3]
            alt_alleles = fields[4].split(',')
            variant_ID = seqid + "<:>" + position
            variant_set[variant_ID] = Variant(ref_allele, alt_alleles, variant_set = variant_set)
            variant_set[variant_ID].ID = fields[2]
            variant_set[variant_ID].quality = fields[5]
            variant_set[variant_ID].filter_status = fields[6]
            variant_set[variant_ID].info = fields[7]
            variant_set[variant_ID].vcf_header_index = vcf_header_index
            if len(fields) > 8:
                genotypes = {}
                sample_features = fields[8].split(':')
                for sample in fields[9:]:
                    sample_name = sample_names[fields.index(sample) - 9]
                    genotypes[sample_name] = Genotype()
                    this_sample_features = sample.split(':')
                    for feature in sample_features:
                        try:
                            setattr(genotypes[sample_name], feature, this_sample_features[sample_features.index(feature)])
                        except IndexError:
                            pass
                variant_set[variant_ID].genotypes = copy.deepcopy(genotypes)
    variant_set.vcf_headers.append('\n'.join(vcf_headerlines))
    if variant_set_to_modify == None:
        return variant_set
