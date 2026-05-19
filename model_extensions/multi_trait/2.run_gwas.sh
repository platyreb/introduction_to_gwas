#!/bin/bash

## prepare data for GWAS

genotypef=../../6.steps/rice_imputed.raw
snpf=../../6.steps/rice_imputed.map
phenotypef=rice_phenotypes_multi.txt
reff=../../cross_reference/rice_group.reference
traits='SL,SW'
covariates='population'

## stand-alone script
# Rscript --vanilla gwas_sommer_multitrait_single_snp.R genotype_file=$genotypef snp_map=$snpf phenotype_file=$phenotypef traits=$traits covariates=$covariates thin=200
Rscript --vanilla gwas_sommer_multitrait_gblup.R genotype_file=$genotypef snp_map=$snpf phenotype_file=$phenotypef traits=$traits covariates=$covariates

