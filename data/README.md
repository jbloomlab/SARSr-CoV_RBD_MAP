# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - [PacBio_amplicons.gb](PacBio_amplicons.gb): the amplicons we previously sequenced by PacBio.
     Note that there are a variety of possible amplicons based on the different unmutated RBD sequences.
     This file is used for loading in the barcode-variant association table for the `count_variants.ipynb` notebook using `dms_variants`

   - [feature_parse_specs.yaml](feature_parse_specs.yaml): how to parse the amplicon when handling the PacBio amplicons

   - [barcode_variant_table.csv](barcode_variant_table.csv): the barcode:homolog look up table, determined in the `SARSr-CoV_homolog_survey` repo. The file was edited to only include unmutated homolog barcode variants, which is the sub-pool that we used for these experiments.
   
   - [barcode_expression.csv](barcode_expression.csv): previously measured per-barcode expression from experiments analyzed in the `SARSr-CoV_homolog_survey` repo. We use these per-barcode expression measurements to filter in our escape score calculations.

   - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples.

   - [RBD_aa_aligned.fasta](RBD_aa.aligned.fasta_): alignment of all of the homologs, as determined in the `SARSr-CoV_homolog_survey` repo.