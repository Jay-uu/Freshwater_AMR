#!/bin/bash

#environments necessary for the AMR-finding nextflow script

mamba activate nextflow
mamba activate --stack checkm_env
mamba activate --stack rgi_env
mamba activate --stack python_venn
mamba activate --stack gtdbtk-2.1.1
mamba activate --stack abricate_env

