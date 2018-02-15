import unittest, os

from annovar_preprocessing import getPairedNumbers, normalToZero, fillInfo, processTableFile, processAllFiles


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
        self.assertEqual("1	6659494	6659494	C	T	.	.	.	.	.	.	.	.	.	.	.	NR:NV	56:0	79:16	80:21	42:0	68:8", outputlines[0])

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





if __name__ == '__main__':
    unittest.main()
