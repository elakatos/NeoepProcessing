
import subprocess

def CheckPeptideNovelty(line):
    peptide = line.split('\t')[-13]
    with open('peptidematch.tmp.log', 'w') as logFile:
        cmd = ['java', '-jar', '/data/home/hfx365/Software/PeptideMatchCMD_1.0.jar', '-a', 'query', '-i', '/data/home/hfx365/Reference/Ensembl/index/','-q', peptide, '-o', 'tmp_peptidematch.out']
        runcmd = subprocess.Popen(cmd, stdout=logFile)
        runcmd.wait()

    with open('tmp_peptidematch.out', 'r') as pmFile:
        lines = pmFile.readlines()
    match = lines[2].strip('\n').split('\t')[1]
    novel = int(match =='No match')
    return(novel)

def ProcessPepmatch(file):
    with open(file, 'r') as pmFile:
        pmFile.readline()
        pmFile.readline() #read first two header lines
        pmDict = {line.split('\t')[0] : line.split('\t')[1].rstrip('\n') for line in pmFile.readlines() }

    return(pmDict)
