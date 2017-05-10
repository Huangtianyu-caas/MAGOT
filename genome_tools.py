#!/usr/bin/python
try:
    import genome
    import genome_tools_config as config
except:
    print "it appears that genome_tools.py is not in the same directory as genome.py and genome_tools_config.py"
import sys
import subprocess


def main():
    program = sys.argv[1]
    arguments = sys.argv[2:]
    command = program + "("
    if program == '-h' or program == '-help' or program == '--help' or program == 'help':
        help_func()
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



def sanitize_pathname(pathname):
    return pathname.replace('|','').replace('<','').replace('>','').replace(':','').replace(';','')


def help_func():
    func_list = []
    for attribute in globals():
        if type(globals()[attribute]).__name__ == "function":
            func_list.append(attribute)
    func_list.remove('main')
    print "\ngenome_tools script from GenomePy.\n\nUsage: python genome_tools.py function [option1=<option1 choice> ...] \n\nFunctions:\n\
    " + '\n    '.join(func_list)


#debug
def print_input(*arg):
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
                out.write('\n'.join(starjunc_dic[seqid]))
        out.write(my_genome.write_apollo_gff(seqid, suppress_fasta = suppress_fasta))
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


def gff2fasta(genome_sequence,gff,from_exons = "False",seq_type = "nucleotide", longest = "False"):
    my_genome = genome.Genome(genome_sequence)
    if from_exons == "True":
        my_genome.read_gff3(gff, features_to_ignore = "CDS", features_to_replace = [('exon','CDS')])
    else:
        my_genome.read_gff3(gff)
    print my_genome.annotations.get_fasta('gene',seq_type = seq_type, longest=eval(longest))


def starjunc2gff(starjunc_file, output = 'stdout'):
    if output == 'stdout':
        outopt = 'print'
    else:
        out = open(output,'w')
        outopt = 'string'
    stargff = genome.starjunc2gff(starjunc_file,output = outopt)
    if outopt == 'string':
        out.write(starff)
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
    


if __name__ == "__main__":
    main()
