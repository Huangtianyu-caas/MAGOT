#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas
#See license on github (https://github.com/biorover/MAGOT/blob/master/LICENSE)
#Comments, feature requests, and friendly encouragement can be offered via github (https://github.com/biorover/MAGOT/issues)

#Known issues:
#Currently UTRs are not handled well. Most ways to import gff3 with UTRs encoded by exon feature regions which don't
#overlap with CDS feature regions will result in the UTRs being lost (because exon features are usually summarily ignored).
#Additionally, even if you manage to import UTRs, most output will likely ignore them. Bug us about it on a github ticket
#and we may well fix it.
#
#Importing GFFs still quite slow, likely having to do with ID assignment. Will try to fix soon.

try:
    import genome
    import genome_tools_config as config
except:
    print "it appears that genome_tools.py is not in the same directory as genome.py and genome_tools_config.py"
import sys
import subprocess
import numpy


def main():
    program = sys.argv[1]
    arguments = sys.argv[2:]
    command = program + "("
    if program == '-h' or program == '-help' or program == '--help' or program == 'help':
        help_func()
        return None
    elif sys.argv[2] in ['-h','--help','-help','help','--h']:
        program_help_func(program)
        return None
    for argument in arguments:
        if '=' in argument:
            argsplit = argument.split('=')
            command = command + argsplit[0] + '="'+ argsplit[1] +'",'
        else:
            command = command + '"' + argument + '",'
    if command[-1] == ',':
        command = command[:-1] + ')'
    else:
        command = command + ')'    
    eval(command)

def program_help_func(program):
    p_func = eval(program)
    p_args = p_func.func_code.co_varnames[:p_func.func_code.co_argcount]
    p_defs = p_func.func_defaults
    if p_func.func_doc:
        p_doc = p_func.func_doc
    else:
        p_doc = ""
    if p_defs:
        def_offset = len(p_defs) - len(p_args)
    else:
        def_offset = -1 * len(p_args)
    arg_string = ""
    for arg_num in range(len(p_args)):
        if p_defs:
            def_index = arg_num + len(p_defs) - len(p_args)
            if def_index >= 0:
                arg_string = arg_string + " " + p_args[arg_num] + "=" + str(p_defs[def_index])
            else:
                arg_string = arg_string + " " + p_args[arg_num]
        else:
            arg_string = arg_string + " " + p_args[arg_num]
    print program + " from the MAGOT genomics toolkit\n\n" + p_doc + "\n\nUsage: " + sys.argv[0] + " " + program + arg_string


def sanitize_pathname(pathname):
    return pathname.replace('|','').replace('<','').replace('>','').replace(':','').replace(';','')


def help_func():
    func_list = []
    for attribute in globals():
        if type(globals()[attribute]).__name__ == "function":
            func_list.append(attribute)
    func_list.remove('main')
    good_boys = ['exclude_from_fasta','coords2fasta','tab2fasta','fasta2tab','gff2fasta','cds2pep','convert_gff','exonerate2gff']
    good_boys.sort()
    experimental = list(set(func_list) - set(good_boys))
    experimental.sort()
    print "\ngenome_tools script from MAGOT.\n\nUsage: " + sys.argv[0] + \
    " function [option1=<option1 choice> ...]\n\nRelatively well-behaved functions:\n    " + '\n    '.join(good_boys) + \
    '\n\nExperimental functions\n    ' +'\n    '.join(experimental)


#debug
def print_input(*arg):
    """This is for debuging and testing. It just prints all arguments"""
    subprocess.call(
    """
    echo """ + '"' + " ".join(arg) + '"',
    shell = True
    )
    #/debug


def nucmer_plot(qgenome_file_loc,tgenome_file_loc):
    subprocess.call('\n'.join([
        config.nucmer +" -l 100 -c 1000 " + tgenome_file_loc + " " + qgenome_file_loc,
        config.dnadiff + " -d out.delta",
        config.mummerplot + " --small --fat --postscript out.1delta",
        config.ps2pdf + " out.ps out.pdf"
        ]),shell = True)


def fqstats(fastq_location):
    """computes basic summary stats for fastq file"""
    lenlist = []
    fastq = open(fastq_location)
    fastq_line = fastq.readline()
    counter = 0
    while fastq_line != "":
        if counter % 4 == 1:
            lenlist.append(len(fastq_line) - 1)
        counter = counter + 1
        fastq_line = fastq.readline()
    read_number = len(lenlist)
    basepairs = sum(lenlist)
    mean_length = basepairs / read_number
    lenlist.sort()
    lenlist.reverse()
    median_length = lenlist[read_number / 2]
    percentile25 = lenlist[read_number / 4]
    percentile75 = lenlist[read_number * 3 / 4]
    nsum = 0
    n25 = False
    n50 = False
    n75 = False
    for length in lenlist:
        nsum = nsum + length
        if nsum > basepairs / 4 and not n25:
            n25 = length
        if nsum > basepairs / 2 and not n50:
            n50 = length
        if nsum > basepairs * 3 / 4 and not n75:
            n75 = length
    print str(read_number) + " reads"
    print "total basepairs=" + str(basepairs)
    print "median length=" + str(median_length)
    print "mean length=" + str(mean_length)
    print "25% of reads >" + str(percentile25)
    print "75% of reads >" + str(percentile75)
    print "n25=" + str(n25)
    print "n50=" + str(n50)
    print "n75=" + str(n75)
    
    
def genewise_wrapper(query_file,genome_file,hmm = False):
    genome = open(genome_file).read().split('>')[1:]
    if hmm:
        hmmopt = '-hmmer '
    else:
        hmmopt = ''
    subprocess.call('mkdir temp', shell = True)
    for sequence in genome:
        seqid = sequence.split('\n')[0].split()[0]
        out = open('temp/' + seqid + '.fa','w')
        out.write('>'+sequence)
        out.close()
        command = config.genewise + hmmopt + query_file + ' temp/' + seqid + '.fa > ' + seqid + '.out'
        subprocess.call(command, shell = True)
    subprocess.call('cat temp/*.out > genwise.out', shell = True)
    subprocess.call('rm -r temp', shell = True)


def dna2orfs(fasta_location,output_file,from_atg = False,longest = False):
    """takes a dna sequence in fasta format and returns ORFs found therein"""
    dna = genome.Genome(fasta_location)
    out = open(output_file, 'w')
    if longest:
        orf_list = []
    for seq in dna.genome_sequence:
        if longest:
            candidate_list = []
            longest_orf_len = 0
        for frame in [0,1,2]:
            for strand in ['-','+']:
                translated_seq_list = dna.genome_sequence[seq].translate(frame=frame,strand=strand).split('*')
                if strand == '+':
                    length_list = [frame] #populates with lengths of orfs so that I can trace back orf positions
                else:
                    length_list = [len(dna.genome_sequence[seq]) - frame]
                for orf in translated_seq_list:
                    orf_start = sum(length_list)
                    length_list.append((1 + len(orf)) * 3)
                    if from_atg:
                        try:
                            output_orf = 'M' + ''.join(orf.split('M')[1:])
                        except IndexError:
                            continue
                    else:
                        output_orf = orf
                    if longest:
                        if len(output_orf) > longest_orf_len:
                            candidate_list.append('>'+seq+'_longestORF\n'+output_orf+'\n')
                            longest_orf_len = len(output_orf)
                    else:
                        out.write('>'+seq+'-pos:'+str(orf_start)+'\n'+output_orf+'\n')
        if longest:
            out.write(candidate_list[-1])
    out.close()


def prep4apollo(genome_sequence, suppress_fasta = "False", output_directory = 'apollo_gffs', exon_fasta = None, full_length_seqs = None,
                               exon_blast_csv = None, exonerate_output = None, starjuncs = None, other_gff = None, other_gff_format = 'gff3',
                               blast_evalue = '0.01', exonerate_percent = '50',output_empty_scaffolds = "False",
                               exonerate_intron_steps = "2000,5000,200000", mapping_threads = "1"):
    """takes evidence inputs and returns gff files to open in apollo"""
    subprocess.call("mkdir -p " + output_directory, shell = True)
    subprocess.call("mkdir -p " + output_directory + "/temp", shell = True)
    mapping_cmds = []
    blast_run = False
    exonerate_run = False
    suppress_fasta = eval(suppress_fasta)
    output_empty_scaffolds = eval(output_empty_scaffolds)
    if exon_fasta != None:
        subprocess.call(config.makeblastdb + ' -in ' + genome_sequence + ' -out ' + output_directory
                        + '/temp/tempdb -dbtype nucl', shell = True)
        mapping_cmds.append(config.tblastn + ' -query ' + exon_fasta + ' -db ' + output_directory + '/temp/tempdb -evalue '
                        + blast_evalue + " -out " + output_directory + "/exon_tblastn.csv -outfmt 10")
        blast_run = True
    if full_length_seqs != None:
        exonerate_intron_lengths = exonerate_intron_steps.split(',')
        for intron_length in exonerate_intron_lengths:        
            mapping_cmds.append(config.exonerate + ' --model protein2genome --percent ' + exonerate_percent + ' --maxintron '
                            + intron_length + ' ' + full_length_seqs + ' ' + genome_sequence + ' > ' + output_directory
                            + '/exonerate_output_' + intron_length + 'bp_introns.txt')
        exonerate_run = True
    running_cmds = []
    if mapping_cmds != []:
        if blast_run:
            print "mapping exons with tblastn"
        if exonerate_run and blast_run:
            print "       and"
        if exonerate_run:
            print "mapping full length sequences with exonerate"
    for cmd_index in range(len(mapping_cmds)):
        running_cmds.append(subprocess.Popen(mapping_cmds[cmd_index],shell = True))
        if (cmd_index + 1) % int(mapping_threads) == 0 or cmd_index == (len(mapping_cmds) - 1):
            for cmd in running_cmds:
                cmd.wait()
            running_cmds = []
    if blast_run:
        if exon_blast_csv != None:
            subprocess.call('cat ' + exon_blast_csv + ' ' + output_directory + '/exon_tblastn.csv > ' + output_directory
                            + '/cat_exon_tblastn.csv', shell = True)
            exon_blast_csv = output_directory + '/cat_exon_tblastn.csv'
        else:
           exon_blast_csv = output_directory + '/exon_tblastn.csv'
    if exonerate_run:
        if exonerate_output != None:
            subprocess.call('cat ' + exonerate_output + ' ' + output_directory + '/exonerate_output* > ' + output_directory
                            + '/cat_exonerate_output.txt', shell = True)
        else:
            subprocess.call('cat ' + output_directory + '/exonerate_output* > ' + output_directory + '/cat_exonerate_output.txt', shell = True)
        exonerate_output = output_directory + '/cat_exonerate_output.txt'  
    print "building apollo gffs"
    my_genome = genome.Genome(genome_sequence,other_gff,annotation_format = other_gff_format)
    if exon_blast_csv != None:
        my_genome.read_blast_csv(exon_blast_csv, find_truncated_locname = True)
    if exonerate_output != None:
        my_genome.read_exonerate(exonerate_output)
    if output_empty_scaffolds:
        seqids = my_genome.get_seqids()
    else:
        seqids = my_genome.annotations.get_all_seqids()
    if starjuncs != None:
        starjunc_dic = {}
        starjunc_list = genome.starjunc2gff(starjuncs,output = "list")
        for junc in starjunc_list:
            seqid = junc.split('\t')[0]
            if seqid in starjunc_dic:
                starjunc_dic[seqid].append(junc)
            else:
                starjunc_dic[seqid] = [junc]
    for seqid in seqids:
        out = open(output_directory + '/' + sanitize_pathname(seqid) + '.gff','w')
        if starjuncs != None:
            if seqid in starjunc_dic:
                out.write('\n'.join(starjunc_dic[seqid]) + '\n')
        out.write(my_genome.write_apollo_gff(seqid, suppress_fasta = suppress_fasta).replace('\ttranscript\t','\tmRNA\t'))
        out.close()
    subprocess.call('rm -rf ' + output_directory + '/temp', shell = True)


def blast_csv2fasta(genome_sequence,blast_csv):
    my_genome = genome.Genome(genome_sequence)
    my_genome.read_blast_csv(blast_csv)
    outfasta = []
    for match in my_genome.annotations.match:
        outfasta.append(my_genome.annotations.match[match].get_fasta())
    print '\n'.join(outfasta)

  
def exonerate2fasta(genome_sequence,exonerate_file):
    my_genome = genome.Genome(genome_sequence)
    my_genome.read_exonerate(exonerate_file)
    outfasta = []
    for match in my_genome.annotations.match:
        outfasta.append(my_genome.annotations.match[match].get_fasta())
    print '\n'.join(outfasta)


def get_CDS_peptides(genome_sequence,gff,output_location,gene_name_filters = [], gene_length_filter = None, names_from = "CDS"):
    my_genome = genome.Genome(genome_sequence)
    my_genome.read_gff3(gff)
    out = open(output_location,'w')
    for gene in my_genome.annotations.gene:
        gene_obj = my_genome.annotations.gene[gene]
        keepgene = True
        for name_filter in gene_name_filters:
            if name_filter in gene_obj.ID:
                keepgene = False
        if gene_length_filter != None:
            seqlen = len(gene_obj.get_fasta().split('\n')[1])
            if seqlen < int(gene_length_filter):
                keepgene = False
        if keepgene:
            for transcript in gene_obj.child_list:
                CDSdict = {}
                transcript_obj = my_genome.annotations.transcript[transcript]
                for CDS in transcript_obj.child_list:
                    CDS_obj = my_genome.annotations.CDS[CDS]
                    CDSdict[CDS_obj.coords] = (CDS_obj.ID,CDS_obj.get_seq().get_orfs(longest = True))
                CDSlist = list(CDSdict)
                CDSlist.sort()
                if transcript_obj.strand == "-":
                    CDSlist.reverse()
                counter = 1
                for CDS in CDSlist:
                    if names_from == 'CDS':
                        pep_name = CDSdict[CDS][0]
                    elif names_from == 'transcript':
                        pep_name = transcript_obj.ID + '-CDS' + str(counter)
                        counter = counter + 1
                    elif names_from == 'gene':
                        pep_name = gene_obj.ID + '-CDS' + str(counter)
                        counter = counter + 1
                    else:
                        print "invalid option for 'names_from' argument"
                        break
                    out.write('>' + pep_name + '\n' + CDSdict[CDS][1] + '\n')


def gff2fasta(genome_sequence,gff,from_exons = "False",seq_type = "nucleotide", longest = "False", genomic = "False",
              truncate_names = "False", name_from = "transcript"):
    my_genome = genome.Genome(genome_sequence, truncate_names = eval(truncate_names))
    if from_exons == "True":
        my_genome.read_gff(gff, features_to_ignore = "CDS", features_to_replace = [('exon','CDS')])
    else:
        my_genome.read_gff(gff)
    if name_from == "transcript":
        print my_genome.annotations.get_fasta('gene',seq_type = seq_type, longest=eval(longest), genomic = eval(genomic))
    elif name_from == "gene":
        print my_genome.annotations.get_fasta('gene',seq_type = seq_type, longest=eval(longest), genomic = eval(genomic), name_from = "Parent")


def starjunc2gff(starjunc_file, output = 'stdout'):
    if output == 'stdout':
        outopt = 'print'
    else:
        out = open(output,'w')
        outopt = 'string'
    stargff = genome.starjunc2gff(starjunc_file,output = outopt)
    if outopt == 'string':
        out.write(stargff)
        out.close()


def tab2fasta(tab_file):
    print genome.tab2fasta(tab_file)
    

def fasta2tab(fasta_file):
    print genome.fasta2tab(fasta_file)


def multithread_exonerate(query_fasta, database_fasta, threads, exonerate_options="--model protein2genome", tempdir = "temp_multithread_exonerate"):
    """splits query_fasta into smaller files (equivalent to # threads) and runs them seperately, then combines the output"""
    query_list = open(query_fasta).read().split('>')[1:]
    n_query_seqs = len(query_list)
    n_per_file = n_query_seqs / int(threads) + 1
    subprocess.call('mkdir '+tempdir, shell = True)
    running_cmds = []
    for i in range(int(threads)):
        chunkstart = i * n_per_file
        chunkstop = (i+1) * n_per_file
        outname = tempdir + '/chunk' + str(i) + '.fasta'
        out=open(outname,'w')
        out.write(">" + ">".join(query_list[chunkstart:chunkstop]))
        out.close()
        running_cmds.append(subprocess.Popen(config.exonerate + " " + exonerate_options + " " + outname + " " + database_fasta + " > "
                                             + outname + ".exonerate",shell = True))
    for cmd in running_cmds:
                cmd.wait()
    subprocess.call("cat " + tempdir + "/*.exonerate",shell = True)
    subprocess.call('rm -rf ' + tempdir, shell = True)

    


def exclude_from_fasta(fasta, exclude_list, just_firstword = "False"):
    """excludes specific fasta entries from fasta file. "exclude_list" can be either comma
    seperated names or name of file with names on each line"""
    my_fasta = genome.Genome(fasta)
    try:
        exlist = open(exclude_list).read().replace('\r','').split('\n')
    except:
        exlist = exclude_list.split(',')
    for seqid in my_fasta.genome_sequence:
        if just_firstword == "True":
            seqid_fixed = seqid.split()[0]
        else:
            seqid_fixed = seqid
        if not seqid_fixed in exlist:
            print '>' + seqid + '\n' + my_fasta.genome_sequence[seqid]
    

def mask_from_gff(genome_sequence,gff,mask_type="soft", overwrite_softmask="True", feature_type="CDS"):
    genome_dict = {}
    genome_file = open(genome_sequence)
    working_seqid = ""
    for line in genome_file:
        if line[0] == ">":
            working_seqid = line.split()[0][1:].replace('\r','').replace('\n','')
            genome_dict[working_seqid] = []
        else:
            if overwrite_softmask in ["True","T","true","t","TRUE"]: 
                genome_dict[working_seqid].extend(list(line.replace('\r','').replace('\n','').upper()))
            elif overwrite_softmask in ["False","F","false","f","FALSE"]:
                genome_dict[working_seqid].extend(list(line.replace('\r','').replace('\n','')))
            else:
                print "Invalid option for 'overwrite_softmask', argument accepts 'True' or 'False'"
                return None
    gff_file = open(gff)
    for line in gff_file:
        if line.count('\t') > 5:
            fields = line.split('\t')
            if fields[2] == feature_type:
                start = int(fields[3])
                stop = int(fields[4])
                if mask_type == "soft":
                    genome_dict[fields[0]][start - 1:stop] = list("".join(genome_dict[fields[0]][start - 1:stop]).lower())
                elif mask_type == "hard":
                    nlist = []
                    for i in range(1 + stop - start):
                        nlist.append('N')
                    genome_dict[fields[0]][start - 1:stop] = nlist
                else:
                    print "Invalid option for mask_type, argument accepts 'soft' and 'hard'"
                    return None
    for seqid in genome_dict:
        print ">" + seqid + '\n' + "".join(genome_dict[seqid])
    


def repeatmasker2augustushints(repeatmasker_gff):
    gff = open(repeatmasker_gff)
    for line in gff:
        if line.count('\t') > 5:
            fields = line.split('\t')
            print "\t".join([fields[0],fields[1],"nonexonpart",fields[3],fields[4],fields[5],".",".","pri=2;src=RM"])


def extract_upstream_downstream(genome_sequence,gff,sequence_length,stream,feature_type = "gene",namefrom = "ID", truncate_names = "True"):
    sequence_dict = genome.GenomeSequence(genome_sequence, truncate_names = eval(truncate_names))
    output_seqs = []
    for line in open(gff):
        if line.count('\t') > 5 and line[0] != "#":
            fields = line.split('\t')
            if fields[2] == feature_type:
                name = None
                coords = [int(fields[3]), int(fields[4])]
                coords.sort()
                for attribute in fields[-1].split(';'):
                    if namefrom == attribute.split('=')[0]:
                        name = attribute.split('=')[1].replace('\r','').replace('\n','')
                if name == None:
                    name = 'seq' + str(len(output_seqs))
                if stream == "up" and fields[6] == "+" or stream == "down" and fields[6] == "-":
                    stop = coords[0] - 1
                    sequence = sequence_dict[fields[0]][stop - int(sequence_length):stop]
                elif stream == "down" and fields[6] == "+" or stream == "up" and fields[6] == "-":
                    start = coords[1]
                    sequence = genome.Sequence(sequence_dict[fields[0]][start:start + int(sequence_length)]).reverse_compliment()
                if len(sequence) == int(sequence_length):
                    output_seqs.append('>' + name + '\n' + sequence)
    print "\n".join(output_seqs)
    

def get_seq_from_fasta(genome_sequence, seq_name, truncate_names = "False"):
    my_genome = genome.Genome(genome_sequence, truncate_names = eval(truncate_names))
    print my_genome.get_scaffold_fasta(seq_name)


def composition_by_site(fasta_alignment):
    seq_dict = genome.GenomeSequence(fasta_alignment)
    position_list = [["a","t","c","g"]]
    for site in range(len(seq_dict[list(seq_dict)[0]])):
        count_dict = {'a':0,'t':0,'c':0,'g':0}
        allcounts = 0
        list_entry = []
        for seq in seq_dict:
            if seq_dict[seq][site].lower() in count_dict:
                count_dict[seq_dict[seq][site].lower()] = count_dict[seq_dict[seq][site].lower()] + 1
                allcounts = allcounts + 1
        for nucleotide in position_list[0]:
            list_entry.append(str(count_dict[nucleotide] * 1.0 / allcounts))
        position_list.append(list_entry)
    for position in position_list:
        print "\t".join(position)


def protMSA2codonMSA(protein_msa_fasta, unaligned_cdna_fasta):
    """takes a protein multiple sequence alignment and the corresponding cDNA (both in fasta format) and returns the cDNA
    aligned accorinding to the protein alignment. Names must be the same in both files. Note: does not perform sanity check.
    If cDNA is more than 3x longer than protein seq, it will be truncated, while if it is less extra gaps will be output at the end."""
    msa_dict = {}
    pep_order =  []
    for seq in open(protein_msa_fasta).read().replace('\r','').split('>')[1:]:
        seqsplit = seq.split('\n')
        msa_dict[seqsplit[0]] = "".join(seqsplit[1:])
        pep_order.append(seqsplit[0])
    cdna_dict = {}
    for seq in open(unaligned_cdna_fasta).read().replace('\r','').split('>')[1:]:
        seqsplit = seq.split('\n')
        cdna_dict[seqsplit[0]] = "".join(seqsplit[1:])
    for seqheader in pep_order:
        print ">" + seqheader
        sequence = ""
        char_num = 0
        for char in msa_dict[seqheader]:
            if char == "-" or char_num >= len(cdna_dict[seqheader]) / 3:
                sequence = sequence + "---"
            else:
                sequence = sequence + cdna_dict[seqheader][char_num * 3: char_num * 3 + 3]
                char_num += 1
        print sequence

   
def at_content_from_fasta(fasta):
    fasta_file = open(fasta)
    firstline = True
    for line in fasta_file:
        if line == "\n":
            continue
        elif '>' in line:
            if firstline:
                name = line[1:].replace('\n','').replace('\r','') + '\t'
                firstline = False
            else:
                print name + str(ats) + '\t' + str(gcs)
                name = line[1:].replace('\n','').replace('\r','') + '\t'
            ats = 0
            gcs = 0
        else:
            ats = ats + line.upper().count('A') + line.upper().count('T')
            gcs = gcs + line.upper().count('G') + line.upper().count('C')
    print name + str(ats) + '\t' + str(gcs)
    

def convert_gff(gff, input_format, output_format):
    """Converts gff from any of the many formats handled by this program to any format this program can output to.
    Currently accepts as input: gff3 (with parent and ID attributes), gtf, augustus, RepeatMasker, CEGMA
    Currently accepts as output: gff3, gtf, apollo_gff3"""
    if input_format == 'gff3':
        presets = None
    else:
        presets = input_format
    if output_format == 'gff3':
        gff_format = "simple gff3"
    elif output_format == 'gtf':
        gff_format = 'gtf'
    elif output_format == "apollo_gff3":
        gff_format = "apollo gff3"
    else:
        print "currently only writes 'gff3' and 'gtf' format"
        return None
    annotations = genome.read_gff(gff, presets = presets)
    print genome.write_gff(annotations, gff_format)


def exonerate2gff(exonerate_output, gff_format='gff3'):
    """Takes the output from an exonerate run (default output format) and returns a gff file formatted as specified"""
    if gff_format == 'gff3':
        write_format = "simple gff3"
    elif gff_format == 'gtf':
        write_format = 'gtf'
    elif gff_format == "apollo_gff3":
        write_format = "apollo gff3"
    else:
        print "currently only writes 'gff3' and 'gtf' format"
        return None
    annotations = genome.read_exonerate(exonerate_output)
    print genome.write_gff(annotations, write_format)

def gff2bed(gff, input_format, columns=12):
    """Takes a gff file and returns a bed file. Number of bed columns can be specified with "columns". \
    Unless your gff defline has "thickStart", "thickEnd", and/or "itemRGB" in it, columns 7-9 will be returned at "." """
    if input_format == 'gff3':
        presets = None
    else:
        presets = input_format
    annotations = genome.read_gff(gff, presets = presets)
    print genome.write_bed(annotations,columns=eval(columns))



def purge_overlaps(gff1, gff_to_purge):
    purge_dict = {}
    for line in open(gff1):
        if line.count('\t') > 6:
            x=line.split('\t')
            if x[0] in purge_dict:
                purge_dict[x[0]].append([int(x[3]),int(x[4])])
            else:
                purge_dict[x[0]] = [[int(x[3]),int(x[4])]]
    for line in open(gff_to_purge):
        printline = True
        if line.count('\t') > 6:
            x=line.split('\t')
            if x[0] in purge_dict:
                coords = [int(x[3]),int(x[4])]
                for purge_coords in purge_dict[x[0]]:
                    if purge_coords[0] < coords[0] < purge_coords[1] or purge_coords[0] < coords[1] < purge_coords[1] or coords[0] < purge_coords[0] < coords[1]:
                        printline = False
                        continue
            if printline:
                print line[:-1]



def besthits_from_psl(psl):
    score_dict = {}
    for line in open(psl):
        if line != "":
            if line[0] not in ['p','\n','m','-','\t',' ']:
                fields = line.split('\t')
                try:
                    qID = fields[9]
                    score = int(fields[0]) - int(fields[1])
                    if qID in score_dict:
                        if score > score_dict[qID][0]:
                            score_dict[qID] = [score,line[:-1]]
                    else:
                        score_dict[qID] = [score,line[:-1]]
                except:
                    print line
                    return None
            else:
                print line[:-1]
    for qID in score_dict:
        print score_dict[qID][1]
    

def depth_from_gff(depth_file, gff, features = "['CDS','intron','upstream','downstream']", stream_length = '1000'):
    my_annotations = genome.read_gff(gff)
    missing_list = []
    features_list = eval(features)
    depth_array_lens = {}
    for line in open(depth_file):
        seqid = line.split()[0]
        if seqid in depth_array_lens:
            depth_array_lens[seqid] = depth_array_lens[seqid] + 1
        else:
            depth_array_lens[seqid] = 1
    depth_arrays = {}
    for seqid in depth_array_lens:
        depth_arrays[seqid] = numpy.zeros(depth_array_lens[seqid],dtype = int)
    for line in open(depth_file):
        fields = line.split()
        if len(fields) > 1:
            depth_arrays[fields[0]][int(fields[1]) - 1] = int(fields[2])
    if "CDS" in features_list:
        print "#CDS counts"
        transcript_CDS_dict = {}
        for CDSID in my_annotations.CDS:
            CDS = my_annotations.CDS[CDSID]
            if CDS.seqid in depth_arrays:
                covsum = sum(depth_arrays[CDS.seqid][CDS.coords[0]-1:CDS.coords[1]])
            else:
                missing_list.append(CDS.parent + ' on ' + CDS.seqid)
            if CDS.parent in transcript_CDS_dict:
                transcript_CDS_dict[CDS.parent] = transcript_CDS_dict[CDS.parent] + covsum
            else:
                transcript_CDS_dict[CDS.parent] = covsum
        for transcript in transcript_CDS_dict:
            print transcript + '\t' + str(transcript_CDS_dict[transcript])
    if 'upstream' in features_list:
        print "#upstream " + stream_length + "bp counts"
        transcript_dict = my_annotations.__dict__[my_annotations[my_annotations.CDS[list(my_annotations.CDS)[0]].parent].feature_type]
        for transcriptID in transcript_dict:
            transcript = transcript_dict[transcriptID]
            if transcript.strand == "+":
                stop = transcript.get_coords()[0] - 1
                start = stop - int(stream_length)
                if start < 0:
                    start = 0
            elif transcript.strand == "-":
                start = transcript.get_coords()[1]
                stop = start + int(stream_length)
            if transcript.seqid in depth_arrays:
                covsum = sum(depth_arrays[transcript.seqid][start:stop])
            else:
                missing_list.append(transcriptID + ' on ' + transcript.seqid)
            print transcriptID + '\t' + str(covsum)
    if len(missing_list) > 0:
        print "#missing from depth file"
        for transcript_info in list(set(missing_list)):
            print transcript_info
        
    
    
def coords2fasta(fasta_file,seqid,start,stop,truncate_names = "False"):
    """prints fasta-format sequence between coordinates (1-based, as in gff-format) within
    a specific entry in a fasta file. truncate_names="True" can be used if you only want to provide
    the first word after the ">" as the seqid (assuming it's unique of course)"""
    print ">" + seqid + ":" + start + "-" + stop
    print genome.Genome(fasta_file, truncate_names=eval(truncate_names)).genome_sequence[seqid][int(start) - 1:int(stop)]

    
def cds2pep(fasta_file, longest = False):
    working_string = ""
    for original_line in open(fasta_file):
        line = original_line.replace('\n','').replace('\r','')
        if line[0] == '>':
            if working_string != "":
                if longest:
                    print genome.Sequence(working_string).get_orfs(longest = True)
                else:
                    print genome.Sequence(working_string).translate()
                working_string = ""
            print line
        else:
            working_string = working_string + line
    print genome.Sequence(working_string).translate()

    
    
    

if __name__ == "__main__":
    main()
