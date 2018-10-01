import annovar_preprocessing as prep
import hla_preprocessing as hla

#prep.processAllFiles("/data2/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp/sample_list.tsv")
hla.composeHLAFile("~/hla_sample_list.tsv", "~/")

#prep.processAllFilesSomatic('Samples/sample_list.tsv')

#gr.processProteomeFasta('/data/home/hfx365/Reference/Ensembl/Homo_sapiens.GRCh38.pep.all.fa', '/data/home/hfx365/Reference/Ensembl/Extracted_proteome.fasta')
