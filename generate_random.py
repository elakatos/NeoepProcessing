from numpy import random


def generateAAsequence(outFileName):
    for fileInd in range(20):
        mutNum = random.randint(100, 800)
        fileName = '/data/BCI-EvoCa2/eszter/Neoepitopes/Random/random_epitopes_' + str(fileInd) + '.fasta'
        with open(fileName, 'w') as outFile:
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