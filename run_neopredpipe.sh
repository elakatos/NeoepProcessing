#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=3:20:0
#$ -l h_vmem=4G

module load python/2.7.14
export PYTHONPATH=/data/home/hfx365/lib/python2.7/site-packages
neopredDir=/data/home/hfx365/Software/NeoPredPipe
dataDir=/data/BCI-EvoCa2/eszter/Neoepitopes/CRCmseq_Polyp

cd $dataDir
#mkdir avready
#mkdir VCF
mv ./raw_annovar/*.avinput ./avready/
mv ./raw_annovar/*.vcf ./VCF

python $neopredDir/main_netMHCpan_pipe.py -I VCF/ -H hlatypes.txt -o Neopred_results -n CRCmseq -c 1 2 3 4 5 6 -E 9 10 -d
