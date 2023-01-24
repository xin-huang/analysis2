"""
Snakefile for running fastsimcoal2 for demographic inference on stdpopsim.
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
        expand(output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.final.demes.png",
            chrms=chrm_list,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            seeds=seed_array,
            ids=[0,1,2],
        )


rule download_fsc2:
    output:
        directory("ext/fsc27_linux64"),
    message:
        "Downloading fsc2",
    shell:
        """
        cd ext
        wget -c http://cmpg.unibe.ch/software/fastsimcoal27/downloads/fsc27_linux64.zip
        unzip http://cmpg.unibe.ch/software/fastsimcoal27/downloads/fsc27_linux64.zip
        cd fsc27_linux64
        chmod u+x fsc27093
        cd ../..
        """


rule run_fsc2:
    input:
        fs = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2_DAFpop0.obs",
        tpl = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.tpl",
        est = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.est",
        fsc2 = rules.download_fsc2.output,
    output:
        bestfits = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.bestfits",
    resources:
        cpus=8,
    params:
        dir = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}",
        repeats = 100, # runs of replicates using fastsimcoal2
        n = 100000, # runs of simulation need to be done to estimate the expected SFS
        L = 40, # runs of ECM cycles
        fsc2 = os.path.abspath("ext/fsc27_linux64/fsc27093"),
    shell:
        """
        cd {params.dir}
        for i in {{1..{params.repeats}}}
        do
            if [ ! -d run$i ]; then mkdir run$i; fi
            cd run$i
            ln -sf {input.fs} sim_{wildcards.chrms}.pop{wildcards.ids}.fastsimcoal2_DAFpop0.obs
            ln -sf {input.tpl} sim_{wildcards.chrms}.pop{wildcards.ids}.fastsimcoal2.tpl
            ln -sf {input.est} sim_{wildcards.chrms}.pop{wildcards.ids}.fastsimcoal2.est 
            {params.fsc2} -t sim_{wildcards.chrms}.pop{wildcards.ids}.fastsimcoal2.tpl -n {params.n} -d -e sim_{wildcards.chrms}.pop{wildcards.ids}.fastsimcoal2.est -M -L {params.L} -q -c {resources.cpus}
            cd ..
        done
        cat run*/sim_{wildcards.chrms}.pop{wildcards.ids}.fastsimcoal2/sim_{wildcards.chrms}.pop{wildcards.ids}.fastsimcoal2.bestlhoods | 
        grep -v MaxObsLhood | sort -rnk 5,5 | sed '1iNPOP\\tTEXP\\tRESIZE\\tANCSIZE\\tMaxEstLhood\\tMaxObsLhood' > {output.bestfits} 
        """


rule plot_fsc2:
    input:
        bestfits = rules.run_fsc2.output.bestfits,
    output:
        yaml = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.final.demes.yaml",
        png = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.final.demes.png",
    shell:
        """
        head -n 2 {input.bestfits} | tail -n 1 |
        awk 'BEGIN{{print "description:"; \
            print "  Piecewise-constant population size model as inferred by fastsimcoal2"; \
            print "time_units: generations"; \
            print "demes:"; \
            print "  - name: generic"; \
            print "    epochs:" \
        }}{{
            print "    - {{end_time: "$2", start_size: "$4"}}"; \
            print "    - {{end_time: 0, start_size: "$1"}}"; \
        }}' > {output.yaml}
        demesdraw size_history --log-time {output.yaml} {output.png}
        """
