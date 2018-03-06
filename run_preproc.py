import annovar_preprocessing as prep
import hla_preprocessing as hla
import generate_random as gr

#prep.processAllFiles("/data2/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp/sample_list.tsv")
#hla.composeHLAFile("/data2/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp/hla_sample_list.tsv", "/data2/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp/")

#prep.processAllFilesSomatic('Samples/sample_list.tsv')

#gr.processProteomeFasta('/data/home/hfx365/Reference/Ensembl/Homo_sapiens.GRCh38.pep.all.fa', '/data/home/hfx365/Reference/Ensembl/Extracted_proteome.fasta')

for i in range(25):
    gr.sampleProteome('/data/home/hfx365/Reference/Ensembl/Extracted_proteome.fasta','/data2/BCI-EvoCa2/eszter/Neoepitopes/Random/random_proteome_'+str(i+1)+'.fasta')
