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


















if __name__ == "__main__":
    main()
