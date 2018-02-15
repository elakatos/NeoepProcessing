import unittest, os

from annovar_preprocessing import getPairedNumbers, normalToZero, fillInfo, processTableFile, processAllFiles
from hla_preprocessing import readInHLAwinners, composeHLAFile


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
        composeHLAFile(samplelist)
        with open("hlatypes.txt", 'r') as testof:
            lines = testof.readlines()

        correctLine = "Test2\thla_a_01_01_01_01\tNA\thla_b_38_01_01\thla_b_14_02_01\thla_c_12_03_01_01\tNA"
        self.assertEqual(correctLine, lines[2].rstrip('\n') )





if __name__ == '__main__':
    unittest.main()
