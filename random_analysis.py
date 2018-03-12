#!/usr/bin/env python

import sys
import os
import subprocess

def DigestAllSamples(listFile, outFile):
    with open(listFile, 'r') as predlist:
        preds = predlist.readlines()
        pm = {'peptidematch_jar':'/data/home/hfx365/Software/PeptideMatchCMD_1.0.jar', 'reference_index': '/data2/home/hfx365/Reference/Ensembl/index'}
        allEps = DigestSample(preds, False, pm)
    with open(outFile, 'w') as outf:
        outf.write('peptide_pos\thla\tpeptide\tpep_core\tOf\tGp\tGl\tIp\tIl\tIcore\tIdentity\tScore\tRank\tCandidate\tBindLevel\tPatIndex\tNovelty\n')
        outf.write(('\n').join(allEps))


def DigestSample(toDigest, checkPeptides, pepmatchPaths):
    '''
    Filters the resulting file and strips all information within it down to individual calls.
    :param toDigest: A list of files to be digested for an individual patient.
    :param patName: Patient/sample identifier
    :return: All Neoantigen Prediction lines free of other information in prediction files.
    '''
    # temp_files = None
    # output_file = "%s%s.digested.txt" % (FilePath, toDigest[0].split('/')[len(toDigest[0].split('/')) - 1].split('.epitopes.')[0])

    lines = []
    pmInputFile = 'tmp/normal.peptidematch.input'
    pmInput = open(pmInputFile,'w')
    for epFile in toDigest:
        print("INFO: Digesting neoantigens for %s" % (patName))
        with open(epFile.rstrip('\n'), 'r') as digest_in:
            patName = epFile.rstrip('\n').replace('.netMHCout', '')
            for line in digest_in:
                line = line.rstrip('\n')
                try:
                    if line.strip()[0].isdigit():
                        linespl = line.split()
                        if '<=' not in linespl:
                            linespl.append('<=\tN')
                        lines.append('\t'.join(linespl)+'\t'+patName)
                        if checkPeptides:
                            pmInput.write('>' + linespl[10] + ';' + linespl[2] + '\n' + linespl[2] + '\n')
                except IndexError as e:
                    pass
    pmInput.close()
    if checkPeptides:
        pmOutFile = 'tmp/normal.peptidematch.out'
        print('Checking peptides against proteome')
        RunPepmatch(pmInputFile, pepmatchPaths['peptidematch_jar'], pepmatchPaths['reference_index'], pmOutFile)
        lines = ProcessPepmatch(pmOutFile, lines)
    print("INFO: Object size of neoantigens: %s Kb"%(sys.getsizeof(lines)))
    return(lines)


def RunPepmatch(pmInput, pepmatchJar, refIndex, pmfileName):
    with open('logForPeptideMatch.tmp', 'a') as logFile:
        cmd = ['java', '-jar', pepmatchJar, '-a', 'query', '-i', refIndex,'-Q', pmInput, '-o', pmfileName]
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

DigestAllSamples('wt_file_list.txt', 'WT.neoantigens.txt')
