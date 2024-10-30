#!/bin/sh
#$ -S /bin/sh
#$ -pe smp 1
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N nf_rna

#### Error Outputfile
#$ -e /home/boeck/myScratch/Pseudomonas/error_outputfile
#$ -o /home/boeck/myScratch/Pseudomonas/error_outputfile

#### Resubmit
#$ -r y


sh /home/boeck/myScratch/Pseudomonas/work/run-nf-core-rnaseq.sh