import unittest, os

from annovar_preprocessing import getPairedNumbers, normalToZero, fillInfo, processTableFile, processAllFiles
from hla_preprocessing import readInHLAwinners, composeHLAFile
from antigen_novelty import ProcessPepmatch, RunPepmatch, RetrieveWT
from generate_random import processProteomeFasta


class TestProcessing(unittest.TestCase):

    def test_read_in_numbers(self):
        oneline = "1\t2\t3\tC\tT\t20\t30\t40\t50\t60\t11\t9\t6\t5\t0"
        self.assertEqual(['20:11', '30:9', '40:6', '50:5', '60:0'], getPairedNumbers(oneline))

    def test_normal_rearrange(self):
        pairs = ['20:11', '30:9', '40:6', '50:5', '60:0']
        self.assertEqual(['60:0', '20:11', '30:9', '40:6', '50:5'], normalToZero(pairs, 4))

    def test_normal_rearrange(self):
        pairs = ['20:11', '30:9', '40:6', '50:5']
        self.assertEqual(pairs, normalToZero(pairs, 4))

    def test_fill_in_info(self):
        oneline = "1\t2\t3\tC\tT\t20\t30\t40\t50\t60\t11\t9\t6\t5\t0"
        self.assertEqual('NR:NV', fillInfo(oneline)[16])

    def test_one_file(self):
        onefile = "test/example1.annoVarInput.txt"
        outputlines = processTableFile(onefile, 4)
        self.assertEqual("54:0", outputlines[3].split('\t')[17])
        correctLine = "1	6659494	6659494	C	T	.	.	.	.	.	.	.	.	.	.	.	NR:NV	56:0	79:16	80:21	42:0	68:8"
        self.assertEqual(correctLine, outputlines[0])

    def test_multiple_sample(self):
        os.system('rm test/example*.avinput')
        os.system('rm test/example*.vcf')
        samplelist = "test/sample_list.tsv"
        processAllFiles(samplelist)
        with open("test/example2.avinput", 'r') as testof:
            lines = testof.readlines()

        self.assertEqual('NR:NV' , lines[0].split('\t')[16])
        self.assertEqual('30:0', lines[1].split('\t')[17])
        self.assertEqual(True, os.path.isfile("test/example1.vcf"))

    def test_read_in_hla(self):
        onefolder = "test/hla_1"
        correctAlleleList = ['hla_a_01_01_01_01','hla_a_29_01_01_01','hla_b_38_01_01','hla_b_14_02_01', 'hla_c_12_03_01_01', 'hla_c_08_02_01']
        self.assertEqual(correctAlleleList, readInHLAwinners(onefolder))

    def test_read_in_hla_NA(self):
        onefolder = "test/hla_2"
        correctAlleleList = ['hla_a_01_01_01_01','NA','hla_b_38_01_01','hla_b_14_02_01', 'hla_c_12_03_01_01', 'NA']
        self.assertEqual(correctAlleleList, readInHLAwinners(onefolder))

    def test_make_hlafile(self):
        os.system("rm hlatypes.txt")
        samplelist = "test/hla_sample_list.tsv"
        composeHLAFile(samplelist, "./")
        with open("hlatypes.txt", 'r') as testof:
            lines = testof.readlines()

        correctLine = "Test2\thla_a_01_01_01_01\tNA\thla_b_38_01_01\thla_b_14_02_01\thla_c_12_03_01_01\tNA"
        self.assertEqual(correctLine, lines[2].rstrip('\n') )

    def test_run_pepmatch(self):
        with open('test/example.eplines', 'r') as epfile:
            eplines = [ line.rstrip('\n') for line in epfile.readlines()]
        RunPepmatch(eplines, '/data/home/hfx365/Software/PeptideMatchCMD_1.0.jar', '/data/home/hfx365/Reference/Ensembl/index/', 'tmp_pepmatch.out')
        with open('peptidematch.tmp.input', 'r') as pminput:
            pmlines = pminput.readlines()

        appendedlines = ['6\tHLA-C*07:02\tTLASKITGM\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t1',
                         '6\tHLA-C*07:02\tASKITGMLL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB\t0',
                         '6\tHLA-C*07:02\tSKITGMLLE\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB\t0',
                         '6\tHLA-C*07:02\tRLFPLIQAL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline196_NM_0025\t0.1744960\t1.6035\t<=\tWB\t1']

        self.assertEqual( ('>line196_NM_0025;RLFPLIQAL','RLFPLIQAL'), (pmlines[6].rstrip('\n'), pmlines[7].rstrip('\n')) )
        self.assertEqual( appendedlines, ProcessPepmatch('tmp_pepmatch.out', eplines) )


    def test_read_pepmatch_file(self):
        pmfileName = 'test/example_pepmatch.out'
        with open('test/example.eplines', 'r') as epfile:
            eplines = [ line.rstrip('\n') for line in epfile.readlines()]
        appendedlines = ['6\tHLA-C*07:02\tTLASKITGM\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t1',
                         '6\tHLA-C*07:02\tASKITGMLL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB\t0',
                         '6\tHLA-C*07:02\tSKITGMLLE\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline195_NM_0025\t0.1744960\t1.6035\t<=\tWB\t0',
                         '6\tHLA-C*07:02\tRLFPLIQAL\tTLASKITGM\t0\t0\t0\t0\t0\tTLASKITGM\tline196_NM_0025\t0.1744960\t1.6035\t<=\tWB\t1']
        self.assertEqual(appendedlines, ProcessPepmatch(pmfileName, eplines))

    def test_proteome_process(self):
        infile = 'test/test.hs.fasta'
        resultfile = 'test/test.hs.processed.fasta'
        processProteomeFasta(infile, resultfile)
        with open(resultfile) as rf:
            lines = rf.readlines()
        correctLines = ['MQRISSLIHLSLFWAGVMSAIELVPEHQTVPVSIGVPATLRCSMKGEAIGNYYINWYRKTQGNTMTFIYREKDIYGPGFKDNFQGDIDIAKNLAVLKILAPSERDEGSYYCACDT\n',
                     'MLSLLHTSTLAVLGALCVYGAGHLEQPQISSTKTLSKTARLECVVSGITISATSVYWYRERPGEVIQFLVSISYDGTVRKESGIPSGKFEVDRIPETSTSTLTIHNVEKQDIATYYCALWEV\n',
                     'MRWALLVLLAFLSPASQKSSNLEGGTKSVTRPTRSSAEITCDLTVINAFYIHWYLHQEGK\n']
        self.assertEqual(correctLines, lines)

    def test_wt_from_tumor(self):
        infile = 'test/example_tumorpep.fasta'
        outfile = 'test/example_normalpep.fasta'
        RetrieveWT(infile, outfile)
        with open(outfile) as rf:
            lines = rf.readlines()
        wtPeps = ['PQSSALTEGDYVPDSPALS\n', 'DHLDAASLQRFLQVEQKMA\n', 'FGRKMDRISSSSGLGCKVL\n']
        self.assertEqual(wtPeps, lines[1::2])



if __name__ == '__main__':
    unittest.main()
