# Landrace control GWAS on Hazel

Working directory: `/share/maize/ntanduk/GrainN/GWAS_env`

## What to copy from the laptop first

    scp final/scripts/server/*            ntanduk@login.hpc.ncsu.edu:/share/maize/ntanduk/GrainN/GWAS_env/
    scp ~/Documents/Research/data/Phenotypes/control_phenotypes_maize_landrace.csv \
        ntanduk@login.hpc.ncsu.edu:/share/maize/ntanduk/GrainN/GWAS_env/
    scp ~/Documents/Github/Grain_Nitrogen/results/PCA/PCA_maize_romerro.csv \
        ntanduk@login.hpc.ncsu.edu:/share/maize/ntanduk/GrainN/GWAS_env/

The VCF is already on the server. Only `all_maize2.vcf.gz` is used: 3,539
accessions, 495,755 markers, the Romero-Navarro landrace panel.

## Steps

    # 0. check the environment has what GAPIT needs
    Rscript install_requirements.R check

    # 1. VCF -> numeric matrix (one-off, ~2 min)
    #    edit the VCF path at the top of 03_prepare.bsub first
    bsub < 03_prepare.bsub

    # 2. soil N first, as the pipeline check
    bsub -J "gwas[1]" < 02_run_gwas.bsub

    #    expect, at Bonferroni 1.009e-07:  MLM 1,  BLINK 38
    grep Bonferroni out/N_0_5/GWAS_models_N_0_5.log

    # 3. the seven control traits
    bsub -J "gwas[2-8]" < 02_run_gwas.bsub

    # 4. consensus, per trait
    GN_DIR_GAPIT=out/latitude \
    GN_GAPIT_TMPL='raw_GWAS_%s_3PC_latitude.txt' \
    GN_FARMCPU_TXT=out/latitude/raw_GWAS_FarmCPU_3PC_latitude.txt \
    GN_OUTDIR=out/latitude/consensus \
    Rscript soilN_crossmodel_consensus.R

## Notes

`01_GWAS_models.R` differs from the laptop version only in how genotypes are
read: HapMap files do not exist on the server, so it slices the numeric matrix
written by `00_vcf_to_matrix.py`. MLM and BLINK run per chromosome; FarmCPU runs
once genome-wide via a filebacked big.matrix, which is the only way it can draw
cofactors from every chromosome.

All traits are restricted by default to the 3,107 accessions that have soil N
(`GN_MATCH_TAXA=N_0_5`), so every scan uses the same panel, kinship and PCs and
the only thing differing between them is the phenotype. Set `GN_MATCH_TAXA=`
to disable.

To redo only the collection step after a scan has finished (seconds, no re-run):

    Rscript 01_GWAS_models.R latitude collect
