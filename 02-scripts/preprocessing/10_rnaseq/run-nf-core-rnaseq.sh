#!/bin/bash

 nextflow run nf-core/rnaseq -r 3.10.1 \
   -profile icbi,singularity \
   -w $(readlink -f /home/boeck/myScratch/Pseudomonas/work) \
   -params-file params.yml -resume