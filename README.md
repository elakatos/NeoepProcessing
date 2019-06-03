# Summary

Helper scripts to do pre- and postprocessing in neoantigen discovery and evaluation. Scripts include preparation of files to be processed by the antigen prediction pipeline of [NeoPredPipe](https://github.com/MathOnco/NeoPredPipe), processing of antigen output tables from NeoPredPipe, processing of HLA modication predictions, and further downstream analysis of immunoediting and immune escape.

(These files should not be run as they are, but rather in parts after appropriate adjustment for file paths and versions.)

## Pre-processing scripts

- **annovar_preprocessing** : Generate files to be directly inputted into the annovar annotation step of [NeoPredPipe](https://github.com/MathOnco/NeoPredPipe), from specific pre-processed mutation tables, generated from Platypus with NR (number of reads at locus) and NV (number of variant reads at locus) information. Outputs a file with '.avinput' extension that can be used in NeoPredPipe.

- **generate_random** : Generate random mutated peptide sequences for neoantigen analysis in fasta format, by introducing mutations to the normal peptidome and extracting the surrounding sequences. Corresponding HLA types are also randomly sampled from the provided file listing HLA types. Outputs a set of fasta files and an hlatypes table-file.

- **random_analysis** : Process neoantigens called using [netMHCpan](http://www.cbs.dtu.dk/services/NetMHCpan/) from randomly introduced mutations of the peptidome. Outputs a neoantigen table with all binding neoantigens.

- **antigen_novelty** : Check whether mutated peptides (as outputted by netMHCpan and processed by NeoPredPipe) can be found in the healthy peptidome using [PeptideMatch](https://research.bioinformatics.udel.edu/peptidematch/index.jsp). NOTE: updated versions of these functions are now incorporated in NeoPredPipe.

- **vcf_filtering** : Functions to process a list of vcf files, filter out non-passing variants and those suspected to be not somatic. Outputs a filtered vcf file for each input file.

## Post-processing and analysis scripts

- **r_functions_postprocessing** : Collection of functions used in processing neoantigen calls and immune escape predictions. 

- **hla_escape_processing** : Process HLA calls to generate appropriate input for the standalone version of netMHCpan (using netMHCpan's nearest neighbour association), subset samples according to HLA type or generate randomly shuffled HLA types. Process alterations of the HLA locus: loss of heterozygosity (as obtained from [lohhla](https://bitbucket.org/mcgranahanlab/lohhla/src/master/)) and somatic mutation predictions (obtained from [polysolver](https://software.broadinstitute.org/cancer/cga/polysolver)). Analyse immune escape status of patients and compile summaries/plots of a patient set.

- **tcga_crc_analysis** : Analysis of immuno-genotypes, immuno-editing and immune evasion of TCGA colorectal cancer samples, as reported in [Lakatos et al., biorXiv, 2019.](https://www.biorxiv.org/content/10.1101/536433v1) Use calls from NeoPredPipe (including NeoRecoPo), clinical information and immune escape calls to explore the interactions between neoantigen burden, hyper-mutator phenotype, immune escape, immune infiltration level and the clonal structure of antigenic mutations. Generate plots and statistics to compare with simulation results obtained with [CloneGrowthSimulation](https://github.com/elakatos/CloneGrowthSimulation).


