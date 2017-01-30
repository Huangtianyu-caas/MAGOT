#!/usr/bin/python
try:
    import genome
    import genome_tools_config
except:
    print "it appears that genome_tools.py is not in the same directory as genome.py and genome_tools_config.py"
import sys
import subprocess

config = genome_tools_config

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
    command = program + "('" + "','".join(arguments) + "')"
    eval(command)


def nucmer_plot(qgenome_file_loc,tgenome_file_loc):
    subprocess.call(
"""
nucmer -l 100 -c 1000 """ + tgenome_file_loc + " " + qgenome_file_loc + """
dnadif -d out.delta
mummerplot --small --fat --postscript out.1delta
ps2pdf out.ps out.pdf""", shell = True
    )


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
        command = 'genewise ' + hmmopt + query_file + ' temp/' + seqid + '.fa > ' + seqid + '.out'
        subprocess.call(command, shell = True)
    subprocess.call('cat temp/*.out > genwise.out', shell = True)
    subprocess.call('rm -r temp', shell = True)

def  dna2orfs(fasta_location,output_file,from_atg = False,longest = False):
    """takes a dna sequence in fasta format and returns ORFs found therein"""
    dna = genome.Genome(fasta_location)
    out = open(output_file, 'w')
    counter = 0
    if longest:
        orf_list = []
    for seq in dna.genome_sequence:
        if longest:
            candidate_list = []
            longest_orf_len = 0
        for frame in [0,1,2]:
            for strand in ['-','+']:
                translated_seq_list = dna.genome_sequence[seq].translate(frame=frame,strand=strand).split('*')
                for orf in translated_seq_list:
                    if from_atg:
                        try:
                            output_orf = 'M' + ''.join(orf.split('M')[1:])
                        except IndexError:
                            continue
                    else:
                        output_orf = orf
                    if longest:
                        if len(ouput_orf) > longest_orf_len:
                            candidate_list.append('>'+seq+'_longestORF\n'+output_orf+'\n')
                            longest_orf_len = len(output_orf)
                    else:
                        out.write('>'+seq+'-seq'+str(counter)+'\n'+output_orf+'\n')
                    counter=counter+1
        if longest:
            out.write(candidate_list[-1])
    out.close()













if __name__ == "__main__":
    main()
