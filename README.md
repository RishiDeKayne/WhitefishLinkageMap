# WhitefishLinkageMap
# This repository will contain all custom scripts used in De-Kayne & Feulner 2018

#All scripts are numbered with output/intermediate files included

#Large files required: i.e. VCF file, lepmap .linkage input file will be uploaded to a data repository
#as part of the paper submission process

#the files are split into two sets. The first, up to the production of the VCF file containings:
# *01_LinkageMapScript_commands.txt*
#this includes all commands run for the processing of raw RAD data, production of the denovo parental reference and genotype calling
            
#and the second for everything following VCF file production i.e. lepmap input preparation, lep map commands and analysis of output
# *02_VCF_to_LepMap_PythonINPUT.R*
#this includes scripts to filter the genotypes in the VCF file by missing data and segregation types in R
# *03_PrepareLepmapInput_fromSexinfo_and_genotypeinfo.ipynb*
#this is an ipython notebook/jupyter notebook file to join the sex-info file made by hand and the genotype file and to encode genotypes as required by lepmap
# *04_LepMapCommands.txt*
#this contains all bash/command line scripts used to carry out linkage mapping with lepmap3
# *05_lepmap_map_organizer.R*
#this is an R script to analyse sex averaged and sex specific maps from lepmap3 mapping
#it also includes all scripts to carry out synteny analysis and plot circlize figure 2 and boxplot figure 3
