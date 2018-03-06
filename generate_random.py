from numpy import random


def generateAAsequence(outFileName):
    aaOneLetter='GPAVLIMCFYWHKRQNEDST'
    for fileInd in range(20):
        mutNum = random.randint(100, 800)
        #outFileName = '/data/BCI-EvoCa2/eszter/Neoepitopes/Random/random_epitopes_' + str(fileInd) + '.fasta'
        with open(outFileName, 'w') as outFile:
            for i in range(mutNum):
                fastaLine = ''.join(random.choice(list(aaOneLetter), 19, replace=True))
                outFile.write('>line' + str(i) + '\n' + fastaLine + '\n')

def processProteomeFasta(fastaName, outName):
    outFile = open(outName, 'w')
    with open(fastaName, 'r') as fastaFile:
        line = fastaFile.readline()
        line = fastaFile.readline()
        while line:
            buffer = ''
            while line[0]!='>':
                buffer = buffer + line.rstrip('\n')
                line = fastaFile.readline()
                if not line:
                    break
            outFile.write(buffer+'\n')
            line=fastaFile.readline()

def codonTable():
    codonDict = {'F': 'SYCLIV', 'L':'FIMVSPHQRW', 'I': 'FLMVTNSKR', 'M':'LIVTKR', 'V':'FLIMADEG',
    'S':'FLYCWPTARGNI', 'P':'STAHQRL', 'T':'SPANKRIM', 'A':'VDEGSPT', 'Y':'FSCHND', 'H':'LPRYNDQ',
    'Q':'LPRKEH', 'N':'YHDITSK', 'K':'QEIMTRN', 'D':'YHNVAGE', 'E':'VAGQKD', 'C':'FSYWRSG',
    'W':'CLSRG', 'R':'CWSGLPHQIMTK', 'G':'CWRSVADE', 'U':'SLCRG', 'O':'LSWQKEY'}
    return(codonDict)

def sampleProteome(protFastaName, outName):
    aaOneLetter='GPAVLIMCFYWHKRQNEDST'
    codonDict=codonTable()
    with open(protFastaName, 'r') as pfasta:
        proteome = pfasta.readlines()
    N = len(proteome)
    with open(outName, 'w') as outFile:
        mutNum = random.randint(100, 800)
        #mutNum = 10
        for i in range(mutNum):
            pID = random.randint(0, N)
            prot = proteome[pID].rstrip('\n')
            if len(prot) < 19:
                pass
            else:
                aaID = random.randint(9, len(prot)-9)
                peptide = prot[(aaID-9):(aaID+10)]
                if '*' in peptide:
                    pass
                else:
                    possibleSubs=codonDict[peptide[9]]
                    new_peptide = peptide[0:9]+str(random.choice(list(possibleSubs)))+peptide[10:]
                    outFile.write('>line'+str(i)+';protein:'+str(pID)+';aa:'+str(aaID)+';mut:'+peptide[9]+'->'+new_peptide[9]+'\n' + new_peptide + '\n')
