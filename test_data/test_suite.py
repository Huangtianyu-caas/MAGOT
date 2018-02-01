#!/usr/bin/env python
#Runs a suite of tests to make sure MAGOT is functioning properly

import sys
import subprocess

command_list = [
    ('exclude_from_fasta O.biroi_refseqGenomeSubset.fasta NW_011924881.1 > temp.test', "1797510917 256187 temp.test")
    ,('coords2fasta C14.fasta Chromosome14 10000 12780 > temp.test','1726369361 2808 temp.test')
    ,('tab2fasta Align.pep.tab > temp.test','925279487 14901 temp.test')
    ,('fasta2tab Align.pep.fa > temp.test','1028007007 14869 temp.test')
#    ,('gff2fasta C14.fasta StandardGTF.gtf > temp.test','2836090577 690750 temp.test')
#    ,('gff2fasta C14.fasta StandardGTF.gtf seq_type=protein > temp.test','111942461 233762 temp.test')
    ,('cds2pep CDSannotations.cds > temp.test','111942461 233762 temp.test')
#    ,('convert_gff minimalGFF3.gff gff3 gtf > temp.test','1904390924 226225 temp.test')
#    ,('convert_gff StandardGTF.gtf gtf gff3 > temp.test','2459709915 283380 temp.test')
#    ,('convert_gff StandardGTF.gtf gtf apollo_gff3 > temp.test','2937090515 531005 temp.test')
#    ,('convert_gff transcriptlessGTF.gtf gtf apollo_gff3 > temp.test','3402448568 564629 temp.test')
#    ,('convert_gff cegma.gff CEGMA apollo_gff3 > temp.test','1350100166 21428 temp.test')
#    ,('convert_gff AugustusGTF.gtf augustus gff3 > temp.test','3088821175 6039 temp.test')
    ,('exonerate2gff Exonerate.txt apollo_gff3 > temp.test','268733671 136536 temp.test')
    ,('gff2bed minimalGFF3.gff gff3 12 > temp.test','')
    ]

for command in command_list:
    subprocess.call('python ../genome_tools.py ' + command[0], shell = True)
    subprocess.call('cksum temp.test > temp.cksum',shell = True)
    cksum = open('temp.cksum').read()[:-1]
    if cksum == command[1]:
        print 'Command: "' + command[0] + '" successful! :)'
    else:
        print 'Command: "' + command[0] + '" failed :(\n\tcksum was:\t\t' + cksum + "\n\tcksum should have been:\t" + command[1]