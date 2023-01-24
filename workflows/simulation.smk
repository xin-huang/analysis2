"""
Snakefile for running simulation on stdpopsim.
"""

import os
import numpy as np
import stdpopsim
import tskit

configfile: "workflows/config/snakemake/config.yaml"

np.random.seed(config["seed"])


# ###############################################################################
# GENERAL RULES & GLOBALS
# ###############################################################################


# The number of replicates of each analysis you would like to run
replicates = config["replicates"]
seed_array = np.random.random_integers(1,2**31,replicates)


# Where you would like all output files from analysis to live
output_dir = os.path.abspath(config["output_dir"])


# The analysis species
species = stdpopsim.get_species(config["species"])


# The names of all chromosomes to simulate, separated by commas
# Use "all" to simulate all chromsomes for the genome
chrm_list = [chrom.id for chrom in species.genome.chromosomes]
if "chrY" in chrm_list:
    chrm_list.remove("chrY")
if(config["chrm_list"] != "all"):
    chrm_list = [chr for chr in config["chrm_list"]]


# The genetic maps you would like to use
genmap_list = config['genetic_map']


# The specific demographic model you would like to run
demo_model_array =  config["demo_models"]
demo_sample_size_dict = {}
demo_ind_start_index_dict = {}
for x in demo_model_array:
    demo_sample_size_dict[x["id"]] = x["num_samples_per_population"]
    demo_ind_start_index_dict[x["id"]] = (np.cumsum(demo_sample_size_dict[x["id"]]) - demo_sample_size_dict[x["id"]][0])/2
demo_model_ids = [x["id"] for x in demo_model_array] 


# Select DFE model from catalog
dfe_list = config["dfe_list"]


# Select annotation file for simulation
annotation_list = config["annotation_list"]


# ###############################################################################
# GENERAL RULES & GLOBALS
# ###############################################################################


rule all:
    input:
        expand(output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.vcf.gz", 
            seeds=seed_array, 
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            chrms=chrm_list,
            ids=[0,1,2]
        )


rule simulation:
    input:
    output:
        ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.trees",
    resources: time=3000, mem_mb=10000,
    run:
        if wildcards.demog == 'Constant': 
            model = stdpopsim.PiecewiseConstantSize(species.population_size)
            mutation_rate = species.genome.mean_mutation_rate
            #mutation_rate = 1.29e-08 # where is this from?
            samples = model.get_samples(*demo_sample_size_dict[wildcards.demog])
        else: 
            model = species.get_demographic_model(wildcards.demog)
            mutation_rate = model.mutation_rate
            samples = model.get_samples(*demo_sample_size_dict[wildcards.demog])  # YRI, CEU, CHB

        genetic_map_id = config["genetic_map"]

        contig = species.get_contig(wildcards.chrms, genetic_map=genetic_map_id)

        if wildcards.dfes != "none":
            # Load dfe only if provided
            dfe = species.get_dfe(wildcards.dfes)
        if wildcards.annots == "all_sites":
            # Adding selection to the whole contig
            contig.add_dfe(intervals=np.array([[0, int(contig.length)]]), DFE=dfe)
        elif wildcards.annots == "none":
            contig = species.get_contig(wildcards.chrms, genetic_map=genetic_map_id)
        else:
            ## Adding annotation only seletion on exon region
            annot = species.get_annotations(wildcards.annots)
            annot_intervals = annot.get_chromosome_annotations(wildcards.chrms)
            contig.add_dfe(intervals=annot_intervals, DFE=dfe)
        
        contig.mutation_rate = mutation_rate
        engine = stdpopsim.get_engine("slim")
        ts = engine.simulate(
            model,
            contig,
            samples,
            slim_scaling_factor=config["slim_scaling_factor"],
            slim_burn_in=config["slim_burn_in"],
            seed=wildcards.seeds,
        )
        ts.dump(output.ts)