"""
Snakefile for running GONE on stdpopsim.
"""

import pathlib
import sys
import os
import numpy as np
import stdpopsim

configfile: "workflows/config/snakemake/config.yaml"

np.random.seed(config["seed"])


# ###############################################################################
# KNOBS -
# ###############################################################################


# The number of replicates of each analysis you would like to run
# For now leaving it a 1 just to get results quickly
replicates = config["replicates"]
seed_array = np.random.random_integers(1,2**31,replicates)


# Where you would like all output files from analysis to live
output_dir = os.path.abspath(config["output_dir"])


# The analysis species
species = stdpopsim.get_species(config["species"])


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


# The specific demographic model you would like to run
demo_model_array =  config["demo_models"]
demo_model_ids = [x["id"] for x in demo_model_array]
demo_sample_size_dict = {}
for x in demo_model_array:
    demo_sample_size_dict[x["id"]] = x["num_samples_per_population"]


# Select DFE model from catalog
dfe_list = config["dfe_list"]


# Select annotation file from stdpopsim
annotation_list = config["annotation_list"]


# ###############################################################################
# GENERAL RULES & GLOBALS
# ###############################################################################


rule all:
    input:
        expand(output_dir + "/inference/demog/gone/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.gone.final.demes.png",
            chrms=chrm_list,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            seeds=seed_array,
            ids=[0,1,2],
        )


rule clone_gone:
    output:
        directory("ext/GONE"),
    message:
        "Cloning GONE repo",
    shell:
        """
        cd ext/
        git clone https://github.com/esrud/GONE.git
        chmod u+x GONE/Linux/PROGRAMMES/*
        cd ..
        """


rule run_gone:
    input:
        ped = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.ped",
        map = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.map",
        gone = rules.clone_gone.output,
    output:
        gone = output_dir + "/inference/demog/gone/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/Output_Ne_sim_{chrms}.pop{ids}.gone",
    params:
        dir = output_dir + "/inference/demog/gone/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/",
        gone_dir = os.path.abspath("ext/GONE/Linux"),
        input_prefix = "sim_{chrms}.pop{ids}.gone",
        input_parameters_file = os.path.abspath("workflows/config/GONE/") + "/INPUT_PARAMETERS_FILE",
    resources:
        cpus = 8,
    shell:
        """
        cd {params.dir}
        ln -sf {params.input_parameters_file} INPUT_PARAMETERS_FILE
        ln -sf {params.gone_dir}/script_GONE.sh script_GONE.sh
        ln -sf {params.gone_dir}/PROGRAMMES PROGRAMMES
        bash script_GONE.sh {params.input_prefix}
        """


rule plot_gone:
    input:
        gone = rules.run_gone.output.gone,
    output:
        yaml = output_dir + "/inference/demog/gone/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.gone.final.demes.yaml",
        png = output_dir + "/inference/demog/gone/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.gone.final.demes.png",
    shell:
        """
        sed '1,2d' {input.gone} | sort -rnk 1,1 | \
        awk 'BEGIN{{print "description:"; \
                    print "  Piecewise-constant population size model as inferred by GONE"; \
                    print "time_units: generations"; \
                    print "demes:"; \
                    print "  - name: generic"; \
                    print "    epochs:" \
            }}{{print "    - {{end_time: "$1", start_size: "$2"}}"}}' > {output.yaml}
        demesdraw size_history --log-time {output.yaml} {output.png}
        """
