#!/bin/bash

snakemake \
    -j 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -c {cluster.cpus} --mem={cluster.mem} -t {cluster.time} -J {cluster.name} --constraint={cluster.constraint}" \
    --latency-wait 30 \
    --use-envmodules \
    --use-conda \
    --conda-prefix /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_RBD_MAP/env/ \
    -R make_summary
