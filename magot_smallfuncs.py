#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas


import copy
import StringIO
from genome import *

def tab2fasta(tab_file):
    tabs = ensure_file(tab_file)
    newlines = []
    for line in tabs:
        fields = line.replace('\n','').replace('\r','').split('\t')
        newlines.append('>'+fields[0]+'\n'+fields[1])
    return "\n".join(newlines)


def fasta2tab(fasta_file):
    fasta = ensure_file(fasta_file)
    newlines = []
    for line in fasta:
        if line[0] == '>':
            newlines.append(line[1:].replace('\n','').replace('\r','') + '\t')
        else:
            newlines[-1] = newlines[-1] + line.replace('\n','').replace('\r','')
    return '\n'.join(newlines)


def ensure_file(potential_file):
    """takes an input that can be either a file location, opened file, or string, and returns an opened file or file like object.
    Used to make downstream applications robust to different input types"""
    if potential_file == None:
        return None
    elif type(potential_file).__name__ == "file":
        return potential_file
    else:
        try:
            return open(potential_file)
        except IOError:
            return StringIO.StringIO(potential_file)


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

def apollo2genome(apollo_gff):
    apollo_list = read_to_string(apollo_gff).split('>')
    gff3 = apollo_list[0]
    fasta = '>' + apollo_list[1]
    return Genome(fasta,gff3,annotation_format='gff3')


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
