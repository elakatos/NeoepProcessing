
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

def RunPepmatch(lines, pepmatchJar, refIndex, pmfileName):
    with open('peptidematch.tmp.input', 'w') as pmInput:
        for line in lines:
            linespl = line.split('\t')
            pmInput.write('>'+linespl[10]+';'+linespl[2]+'\n'+linespl[2]+'\n')

    with open('peptidematch.tmp.log', 'w') as logFile:
        cmd = ['java', '-jar', pepmatchJar, '-a', 'query', '-i', refIndex,'-Q', 'peptidematch.tmp.input', '-o', pmfileName]
        runcmd = subprocess.Popen(cmd, stdout=logFile)
        runcmd.wait()



def ProcessPepmatch(pmfileName, epLines):
    with open(pmfileName, 'r') as pmFile:
        pmFile.readline()
        pmFile.readline() #read first two header lines
        pmDict = {line.split('\t')[0] : line.split('\t')[1].rstrip('\n') for line in pmFile.readlines() }
    appendedLines = []
    for line in epLines:
        epkey = line.split('\t')[10]+';'+line.split('\t')[2]
        novel = int(pmDict[epkey]=='No match')
        appendedLines.append(line+'\t'+str(novel))

    return(appendedLines)
