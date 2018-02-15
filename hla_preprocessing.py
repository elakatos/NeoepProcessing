def readInHLAwinners(hladir):
    hlaFileName = hladir.rstrip('/')+'/winners.hla.txt'
    with open(hlaFileName, 'r') as hlafile:
        hlaList = []
        for line in hlafile.readlines():
            lineArray = line.rstrip('\n').split('\t')
            hlaList.append(lineArray[1])
            if lineArray[2]!=lineArray[1]:
                hlaList.append(lineArray[2])
            else:
                hlaList.append('NA')
    return(hlaList)

def composeHLAFile(sampleListFile):
    with open('hlatypes.txt', 'w') as outFile:
        outFile.write('Patient\tHLA-A_1\tHLA-A_2\tHLA-B_1\tHLA-B_2\tHLA-C_1\tHLA-C_2\n')
        with open(sampleListFile, 'r') as slFile:
            for sample in slFile.readlines():
                hlaDir = sample.rstrip('\n').split('\t')[0]
                sampleID = sample.rstrip('\n').split('\t')[1]
                hlaList = readInHLAwinners(hlaDir)
                outFile.write(sampleID+'\t'+ ('\t').join(hlaList)+'\n')

