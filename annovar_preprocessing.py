
def fillInfo(line):
    infoarray = line.split('\t')[:5]
    {infoarray.append('.') for i in range(11)}
    infoarray.append('NR:NV') # the genotype info pairs are read number and variant read number
    return(infoarray)


def getPairedNumbers(line):
    readsarray = line.split('\t')[5:]
    nsample = int(len(readsarray)/2)
    pairedNumbers = []
    for i in range(nsample):
        pairedNumbers.append(readsarray[i]+':'+readsarray[i+nsample])
    return(pairedNumbers)

def normalToZero(numberarray, normalindex):
    numberarray.insert(0, numberarray.pop(normalindex))
    return(numberarray)

def processTableFile(tableFile, normalindex):
    outputLines = []
    with open(tableFile, 'r') as tableInfile:
        for line in tableInfile:
            infoArray = fillInfo(line.strip('\n'))
            tableArray = getPairedNumbers(line.strip('\n'))
            rearArray = normalToZero(tableArray, normalindex)
            outputLines.append(('\t').join(infoArray)+'\t'+('\t').join(rearArray))
    return(outputLines)

def processAllFiles(sampleListFile):
    with open(sampleListFile, 'r') as slFile:
        for line in slFile:
            tableFile = line.strip('\n').split('\t')[0]
            normIndex = int(line.strip('\n').split('\t')[1])
            print("Processing file %s with normal sample at %s" %(tableFile, normIndex))
            processedLines = processTableFile(tableFile, normIndex)
            outfileName = tableFile.split('_')[0]+'_processed.avinput'
            outfile = open(outfileName, 'w')
            outfile.write(('\n').join(processedLines))
            outfile.close()
            print("Processed file saved at %s" %(outfileName))

