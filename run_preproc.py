import annovar_preprocessing as prep
import hla_preprocessing as hla

#prep.processAllFiles("/data2/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp/sample_list.tsv")
#hla.composeHLAFile("/data2/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp/hla_sample_list.tsv", "/data2/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp/")

prep.processAllFilesSomatic('Samples/sample_list.tsv')