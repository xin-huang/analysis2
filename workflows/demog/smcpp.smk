"""
Snakefile for running smc++ on stdpopsim.
"""

import sys
import os
import numpy as np
import stdpopsim


# ###############################################################################
# KNOBS -
# ###############################################################################

configfile: "workflows/config/snakemake/config.yaml"
output_dir = os.path.abspath(config["output_dir"])

# A seed to replicate results
# TODO mutation rates
np.random.seed(config["seed"])

# The number of replicates of each analysis you would like to run
# For now leaving it a 1 just to get results quickly
replicates = config["replicates"]

# The analysis species
species = stdpopsim.get_species(config["species"])

# This is the number of samples to simulate for within each population
# for each replicate
# TODO double check this is up to date with stdpopsim backend
population_id = config["population_id"]
num_sampled_genomes_per_replicate = config["num_sampled_genomes_per_replicate"]

# The DFE id used for selection analyses
dfe_id = config["dfe_list"][0] # need to generalize to more than one...

# The names of all chromosomes to simulate, separated by commas
# Use "all" to simulate all chromosomes for the genome
chrm_list = [chrom.id for chrom in species.genome.chromosomes]
if "chrY" in chrm_list:
    chrm_list.remove("chrY")
if(config["chrm_list"] != "all"):
    chrm_list = [chr for chr in config["chrm_list"]]

# The genetic maps you would like to use
# if value None is given default_recombination_rates are
# used with a flat map
genetic_map_id = config["genetic_map"]
genmap_list = config['genetic_map']

# This grabs the default mr from the first chromosome,
# Ultimitely This needs to be replaced with the weighted average
# of all chromosomes: This should be done in stdpopsim.
mutation_rate = species.genome.mean_mutation_rate

# The specific demographic model you would like to run
demo_model_array =  config["demo_models"]
demo_model_ids = [x["id"] for x in demo_model_array]
demo_sample_size_dict = {}
for x in demo_model_array:
    demo_sample_size_dict[x["id"]] = x["num_samples_per_population"]

# Select DFE model from catalog
dfe_list = config["dfe_list"]
annotation_list = config["annotation_list"]

# ###############################################################################
# GENERAL RULES & GLOBALS
# ###############################################################################

seed_array = np.random.random_integers(1,2**31,replicates)

rule all:
    input:
        expand(output_dir + "/inference/demog/smcpp/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.smc.Ne.png",
            chrms=chrm_list,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            seeds=seed_array,
            ids=[0,1,2],
        )


rule clone_smcpp:
    output:
        "ext/smcpp/pyproject.toml"
    message: "Cloning SMC++"
    shell:
        """
        cd ext/
        git clone https://github.com/popgenmethods/smcpp.git
        cat smc_setup_stdpopsim_patch > smcpp/setup.py
        cat smc_pyproject_stdpopsim_patch > smcpp/pyproject.toml
        cd smcpp/
        pip install .
        cd ..
        """


rule vcf2smc:
    input:
        vcf = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.vcf.gz",
        ind = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.ind.list",
        smc = rules.clone_smcpp.output, 
    output:
        smc = output_dir + "/inference/demog/smcpp/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.smc.gz",
    shell:
        """
        smc++ vcf2smc {input.vcf} {output.smc} {wildcards.chrms} pop{wildcards.ids}:`awk 'BEGIN{{ORS=","}}{{print $0}}' {input.ind} | sed 's/,$//'`
        """


rule run_smcpp:
    input:
        smc = rules.vcf2smc.output.smc,
    output:
        res = output_dir + "/inference/demog/smcpp/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.smc.final.json",
    params:
        path = output_dir + "/inference/demog/smcpp/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/",
        base = "sim_{chrms}.pop{ids}.smc", 
        nonseg_cutoff = 50000,  # See https://github.com/popgenmethods/smcpp/issues/165
        timepoints = "1e2 1e5", # See https://github.com/popgenmethods/smcpp/issues/165
        mut_rate = lambda wildcards: species.genome.mean_mutation_rate if wildcards.demog == 'Constant' else species.get_demographic_model(wildcards.demog).mutation_rate,
    resources:
        cpus = 8,
    shell:
        """
        cd {params.path}
        smc++ estimate --timepoints {params.timepoints} --base {params.base} --nonseg-cutoff {params.nonseg_cutoff} --cores {resources.cpus} {params.mut_rate} {input.smc}
        """


rule plot_smcpp:
    input:
        res = rules.run_smcpp.output.res, 
    output:
        png = output_dir + "/inference/demog/smcpp/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.smc.Ne.png",
    params:
        generation_time = lambda wildcards: species.generation_time if wildcards.demog == 'Constant' else species.get_demographic_model(wildcards.demog).generation_time,
    shell:
        """
        smc++ plot {output.png} {input.res} -g {params.generation_time} -c
        """