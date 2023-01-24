"""
Snakefile for running stairwayplot2 on stdpopsim.
"""

import pathlib
import sys
import os
import numpy as np
import stdpopsim
import tskit
import ts2fs

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
        expand(output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.final.demes.png",
            chrms=chrm_list,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            seeds=seed_array,
            ids=[0,1,2],
        )

rule clone_stairwayplot2:
    output:
        directory("ext/stairway-plot-v2"),
    message:
        "Cloning stairwayplot2",
    shell:
        """
        cd ext
        git clone https://github.com/xiaoming-liu/stairway-plot-v2.git
        cd stairway-plot-v2
        unzip stairway_plot_v2.1.1.zip
        cd ../../
        """

rule ts2blueprint:
    input:
        ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.trees",
        stairwayplot2 = rules.clone_stairwayplot2.output,
    output:
        blueprint = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.blueprint",
    params:
        popid = "pop{ids}", # id of the population
        nseq = 100, # number of sequences
        L = lambda wildcards: species.get_contig(wildcards.chrms, genetic_map=wildcards.genmap).length, # total number of observed nucleic sites, including polymorphic and monomorphic
        whether_folded = "false", # whethr the SFS is folded (true or false)
        pct_training = 0.67, # percentage of sites for training
        nrand = "7 15 22 28", # number of random break points for each try (separated by white space)
        project_dir = "demography", # project directory
        stairway_plot_dir = "stairway_plot_es", # directory to the stairway plot files
        ninput = 200, # number of input files to be created for each estimation
        mu = lambda wildcards: 1.29e-8 if wildcards.demog == 'Constant' else species.get_demographic_model(wildcards.demog).mutation_rate, # assumed mutation rate per site per generation
        year_per_generation = lambda wildcards: species.generation_time if wildcards.demog == 'Constant' else species.get_demographic_model(wildcards.demog).generation_time, # assumed generation time (in years)
        plot_title = "sim_{chrms}.pop{ids}", # title of the plot
        xrange = "0.1,10000", # Time (1k year) range; format: xmin,xmax; "0,0" for default
        yrange = "0,0", # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
        xspacing = 2, # X axis spacing
        yspacing = 2, # Y axis spacing
        fontsize = 12, # Font size
    run:
        ts = tskit.load(input.ts)
        index = int(wildcards.ids)
        samps = ts.samples(population=index)
        if len(samps) == 0: samps = ts.samples(population=0)
        fs = ts2fs.generate_fs(ts, samps, output=None, format='stairwayplot2')
        fs = " ".join([str(x) for x in fs[1:-1]])

        with open(output.blueprint, 'w') as o:
            o.write(f"popid: {params.popid}\n")
            o.write(f"nseq: {params.nseq}\n")
            o.write(f"L: {params.L}\n")
            o.write(f"whether_folded: {params.whether_folded}\n")
            o.write(f"SFS: {fs}\n")
            o.write(f"pct_training: {params.pct_training}\n")
            o.write(f"nrand: {params.nrand}\n")
            o.write(f"project_dir: {params.project_dir}\n")
            o.write(f"stairway_plot_dir: {params.stairway_plot_dir}\n")
            o.write(f"ninput: {params.ninput}\n")
            o.write(f"mu: {params.mu}\n")
            o.write(f"year_per_generation: {params.year_per_generation}\n")
            o.write(f"plot_title: {params.plot_title}\n")
            o.write(f"xrange: {params.xrange}\n")
            o.write(f"yrange: {params.yrange}\n")
            o.write(f"xspacing: {params.xspacing}\n")
            o.write(f"yspacing: {params.yspacing}\n")
            o.write(f"fontsize: {params.fontsize}\n")      
        
rule run_stairwayplot2:
    input:
        blueprint = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.blueprint",
    output:
        summary = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/demography/sim_{chrms}.pop{ids}.final.summary",
    params:
        dir = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/",
        stairway_plot_es_dir = os.path.abspath("ext/stairway-plot-v2/stairway_plot_v2.1.1/stairway_plot_es"),
    shell:
        """
        cd {params.dir}
        ln -sf {params.stairway_plot_es_dir} stairway_plot_es
        java -cp stairway_plot_es Stairbuilder {input.blueprint}
        bash {input.blueprint}.sh
        """

rule plot_stairwayplot2:
    input:
        summary = rules.run_stairwayplot2.output.summary,
    output:
        yaml = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.final.demes.yaml",
        png = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.final.demes.png",
    params:
        generation_time = lambda wildcards: int(species.generation_time) if wildcards.demog == 'Constant' else int(species.get_demographic_model(wildcards.demog).generation_time),
    shell:
        """
        sed '1d' {input.summary} | sort -rnk 6,6 | \
        awk -v gen={params.generation_time} '{{print int($6/gen)"\\t"$7}}' | sort -urnk 1,1 | \
        awk 'BEGIN{{print "description:"; \
                print "  Piecewise-constant population size model as inferred by stairway-plot-v2"; \
                print "time_units: generations"; \
                print "demes:"; \
                print "  - name: generic"; \
                print "    epochs:" \
        }}{{print "    - {{end_time: "$1", start_size: "$2"}}"}}' > {output.yaml}
        demesdraw size_history --log-time {output.yaml} {output.png}
        """
