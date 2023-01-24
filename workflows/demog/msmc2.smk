"""
Snakefile for running msmc2 on stdpopsim.
"""

import pathlib
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
genetic_map_downloaded_flag= ".genetic_map_downloaded"

rule all:
    input:
        expand(output_dir + "/inference/demog/msmc2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.msmc2.final.demes.png",
            chrms=chrm_list,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            seeds=seed_array,
            ids=[0,1,2],
        )

rule clone_msmc2:
    output: 
        directory("ext/msmc2"),
    message: 
        "Cloning msmc2",
    shell:
        """
        cd ext
	git clone https://github.com/stschiff/msmc2.git
	cat msmc2_makefile_stdpopsim_patch > msmc2/Makefile
	cd msmc2
	make
        cd ../../
        """

rule clone_msmc_tools:
    output:
        directory("ext/msmc-tools"),
    message:
        "Cloning msmc-tools",
    shell:
        """
        cd ext
        git clone https://github.com/stschiff/msmc-tools.git
        cd ../../
        """

rule vcf2multihetsep:
    input:
        vcf = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.vcf.gz",
        msmc2 = rules.clone_msmc2.output,
        msmc_tools = rules.clone_msmc_tools.output,
    output:
        multihetsep = output_dir + "/inference/demog/msmc2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.msmc2.multihetsep.txt",
    shell:
        """
        python ext/msmc-tools/generate_multihetsep.py --chr {wildcards.chrms} {input.vcf} > {output.multihetsep}
        """

rule run_msmc2:
    input:
        multihetsep = rules.vcf2multihetsep.output.multihetsep,
    output:
        final = output_dir + "/inference/demog/msmc2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.msmc2.final.txt",
    params:
        output_prefix = output_dir + "/inference/demog/msmc2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.msmc2",
        time_seg_pattern = "1*2+1*5+1*2+1*3", # See https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md#estimating-the-effective-population-size and https://github.com/stschiff/msmc2/issues/21
        msmc2_exec = os.path.abspath(config["msmc2_exec"]),
    resources:
        cpus = 8,
    shell:
        """
        {params.msmc2_exec} -t {resources.cpus} -p {params.time_seg_pattern} -o {params.output_prefix} {input.multihetsep}
        """

rule plot_msmc2:
    input:
        res = rules.run_msmc2.output.final,
    output:
        yaml = output_dir + "/inference/demog/msmc2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.msmc2.final.demes.yaml", 
        png = output_dir + "/inference/demog/msmc2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.msmc2.final.demes.png",
    params:
        mut_rate = lambda wildcards: species.genome.mean_mutation_rate if wildcards.demog == 'Constant' else species.get_demographic_model(wildcards.demog).mutation_rate,
    shell:
        """
        python ext/msmc-tools/convert_msmc_to_demes.py {input.res} {params.mut_rate} > {output.yaml}
        demesdraw size_history --log-time {output.yaml} {output.png}
        """