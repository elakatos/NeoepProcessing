
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
            tableFile = line.rstrip('\n').split('\t')[0]
            normIndex = int(line.rstrip('\n').split('\t')[1])
            print("Processing file %s with normal sample at %s" %(tableFile, normIndex))
            processedLines = processTableFile(tableFile, normIndex)
            outfileName = tableFile.replace('.annoVarInput.txt','.avinput')
            with open(outfileName, 'w') as outfile:
                outfile.write(('\n').join(processedLines))
            print("Processed file saved at %s" %(outfileName))

            #Create decoy VCF file
            vcffileName = outfileName.replace('.avinput', '.vcf')
            with open(vcffileName, 'w') as decoyvcf:
                decoyvcf.write("#VCF - none\n Annovar ready sample: "+outfileName)

