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