"""
Snakefile for running dadi for demographic inference on stdpopsim.
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
        expand(output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.{types}.final.demes.png",
            chrms=chrm_list,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            seeds=seed_array,
            ids=[0,1,2],
            types=["folded", "unfolded"],
        )


rule run_dadi:
    input:
        fs = output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.{types}.fs",
    output:
        bestfits = output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.{types}.two_epoch.InferDM.bestfits",
    resources:
        cpus=8,
    params:
        demog = 'two_epoch',
        demog_p0 = '5 5',
        demog_ubounds = '10 1',
        demog_lbounds = '10e-3 10e-3',
        grid_size = '300 400 500',
        prefix = output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.{types}.two_epoch",
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {params.demog} --p0 {params.demog_p0} --ubounds {params.demog_ubounds} --lbounds {params.demog_lbounds} --output-prefix {params.prefix} --grids {params.grid_size} --cpus {resources.cpus} --force-convergence --nomisid
        """


rule plot_dadi:
    input:
        bestfits = rules.run_dadi.output.bestfits,
    output:
        yaml = output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.{types}.final.demes.yaml",
        png = output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.{types}.final.demes.png",
    params:
        mut_rate = lambda wildcards: species.genome.mean_mutation_rate if wildcards.demog == 'Constant' else species.get_demographic_model(wildcards.demog).mutation_rate,
        seq_len = lambda wildcards: species.get_contig(wildcards.chrms, genetic_map=wildcards.genmap).recombination_map.sequence_length,
    shell:
        """
        grep Log {input.bestfits} -A 1 | head -n 2 | tail -n 1 | \
        awk -v mut={params.mut_rate} -v len={params.seq_len} \
            'BEGIN{{print "description:"; \
                print "  Piecewise-constant population size model as inferred by dadi"; \
                print "time_units: generations"; \
                print "demes:"; \
                print "  - name: generic"; \
                print "    epochs:" \
            }}{{Na=$4/(4*mut*len); \
                print "    - {{end_time: "$3*2*Na", start_size: "Na"}}"; \
                print "    - {{end_time: 0, start_size: "$2*Na"}}"; \
            }}' > {output.yaml}
        demesdraw size_history --log-time {output.yaml} {output.png}
        """