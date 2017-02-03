#!/usr/bin/python
try:
    import genome
    import genome_tools_config
except:
    print "it appears that genome_tools.py is not in the same directory as genome.py and genome_tools_config.py"
import sys
import subprocess

config = genome_tools_config


def sanitize_pathname(pathname):
    return pathname.replace('|','').replace('<','').replace('>','').replace(':','').replace(';','')


#debug
def print_input(*arg):
    subprocess.call(
    """
    echo """ + '"' + " ".join(arg) + '"',
    shell = True
    )
    #/debug


def main():
    program = sys.argv[1]
    arguments = sys.argv[2:]
    command = program + "("
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


def nucmer_plot(qgenome_file_loc,tgenome_file_loc):
    subprocess.call('\n'.join([
        config.nucmer +" -l 100 -c 1000 " + tgenome_file_loc + " " + qgenome_file_loc,
        config.dnadif + " -d out.delta",
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


def  dna2orfs(fasta_location,output_file,from_atg = False,longest = False):
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
                        out.write('>'+seq+'-pos:'+str(orfstart)+'\n'+output_orf+'\n')
        if longest:
            out.write(candidate_list[-1])
    out.close()


def prep4apollo(genome_sequence, output_directory = 'apollo_gffs', exon_fasta = None, full_length_seqs = None,
                               exon_blast_csv = None, exonerate_output = None, other_gff = None, other_gff_format = 'gff3',
                               blast_evalue = '0.01', exonerate_percent = '50',output_empty_scaffolds = False, exonerate_intron_steps="2000,5000,200000"):
    """takes evidence inputs and returns gff files to open in apollo"""
    subprocess.call("mkdir -p " + output_directory, shell = True)
    subprocess.call("mkdir -p " + output_directory + "/temp", shell = True)
    if exon_fasta != None:
        print "Running tblastn to map exons from exon_fasta to genome"
        subprocess.call(config.makeblastdb + ' -in ' + genome_sequence + ' -out ' + output_directory
                        + '/temp/tempdb -dbtype nucl', shell = True)
        subprocess.call(config.tblastn + ' -query ' + exon_fasta + ' -db ' + output_directory + '/temp/tempdb -evalue '
                        + blast_evalue + " -out " + output_directory + "/exon_tblastn.csv -outfmt 10", shell = True)
        if exon_blast_csv != None:
            subprocess.call('cat ' + exon_blast_csv + ' ' + output_directory + '/exon_tblastn.csv > ' + output_directory
                            + '/cat_exon_tblastn.csv', shell = True)
            exon_blast_csv = output_directory + '/cat_exon_tblastn.csv'
        else:
            exon_blast_csv = output_directory + '/exon_tblastn.csv'
    if full_length_seqs != None:
        print "Running exonerate to map full sequences to genome"
        exonerate_intron_lengths = exonerate_intron_steps.split(',')
        for intron_length in exonerate_intron_lengths:        
            subprocess.call(config.exonerate + ' --model protein2genome --percent ' + exonerate_percent + ' --maxintron '
                            + intron_length + ' ' + full_length_seqs + ' ' + genome_sequence + ' > ' + output_directory
                            + '/exonerate_output_' + intron_length + 'bp_introns.txt', shell = True)
        if exonerate_output != None:
            subprocess.call('cat ' + exonerate_output + ' ' + output_directory + '/exonerate_output* > ' + output_directory
                            + '/cat_exonerate_output.txt', shell = True)
        else:
            subprocess.call('cat ' + output_directory + '/exonerate_output* > ' + output_directory + '/cat_exonerate_output.txt', shell = True)
        exonerate_output = output_directory + '/cat_exonerate_output.txt'
    print "building apollo gffs"
    my_genome = genome.Genome(genome_sequence,other_gff,annotation_format = other_gff_format)
    if exon_blast_csv != None:
        my_genome.read_blast_csv(exon_blast_csv)
    if exonerate_output != None:
        my_genome.read_exonerate(exonerate_output)
    if output_empty_scaffolds:
        seqids = my_genome.get_seqids()
    else:
        seqids = my_genome.annotations.get_all_seqids()
    for seqid in seqids:
        out = open(output_directory + '/' + sanitize_pathname(seqid) + '.gff','w')
        out.write(my_genome.write_apollo_gff(seqid))
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
    my_genome = genome.Genome(genome_sequence,gff,annotation_format = 'gff3')
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

    







if __name__ == "__main__":
    main()
