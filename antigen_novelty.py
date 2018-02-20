
import subprocess

oneline = 'Set.09.Proximal.snv\t1\t1\t1\t1\t1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\tline15\t1\t1403821\tG\tA\tATAD3C:NM_001039211\t4\tHLA-A*01:01\tMMDACMQDF\tMMDACMQDF\t0\t0\t0\t0\t0\tMMDACMQDF\tline15_NM_00103\t0.1579180\t1.5976\t<=\tWB'

def CheckPeptideNovelty(line):
    peptide = line.split('\t')[-13]
    cmd = ['java', '-jar', '~/Software/PeptideMatchCMD_1.0.jar', '-a', 'query', '-i', '~/Reference/Ensembl/index/',
           '-q', peptide, '-o', 'tmp_peptidematch.out']
    runcmd = subprocess.Popen(cmd, stdout='peptidematch.tmp.log')
    runcmd.wait()

    with open('tmp_peptidematch.out', 'r') as pmFile:
        lines = pmFile.readlines()
    match = lines[2].strip('\n').split('\t')[1]
    novel = int(match =='No match')
    print(match)
    return(novel)