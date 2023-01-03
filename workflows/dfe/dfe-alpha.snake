"""
DFE-alpha pipeline for running DFE benchmark on stdpopsim.
"""

import os
import numpy as np
import stdpopsim
import tskit
import ts2fs

configfile: "workflows/config/snakemake/config.yaml"

np.random.seed(config["seed"])


# ###############################################################################
# GENERAL RULES & GLOBALS
# ###############################################################################

# The number of replicates of each analysis you would like to run
replicates = config["replicates"]
seed_list = np.random.random_integers(1,2**31,replicates)

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
demo_model_id_list = [x["id"] for x in demo_model_array] 
demo_sample_size_dict = {}
for x in demo_model_array:
    demo_sample_size_dict[x["id"]] = x["num_samples_per_population"]

# Select DFE model from catalog  
dfe_list = config["dfe_list"]   
annotation_list = config["annotation_list"]


# ###############################################################################
# DFE-alpha
# ###############################################################################


rule all:
    input:
        output_dir + "/inference/dfe/bestfits/DFE-alpha.bestfits.txt",


rule generate_dfe_alpha_fs:
    input:
        ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.trees",
    output:
        neu_config = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.DFE-alpha.neu.config",
        nonneu_config = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.DFE-alpha.nonneu.config",
        fs = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.DFE-alpha.fs",
    params:
        output = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}",
        dfe_alpha_data_path_1 = os.path.abspath(config["dfe_alpha_data_path_1"]), # Path for the data files for the one and two epoch models from DFE-alpha 
        dfe_alpha_data_path_2 = os.path.abspath(config["dfe_alpha_data_path_2"]), # Path for the data files for the three epoch model from DFE-alpha
    threads: 1,
    run:
        if wildcards.genmap == 'Uniform_GeneticMap': contig = species.get_contig(wildcards.chrms)
        else: contig = species.get_contig(wildcards.chrms, genetic_map=wildcards.genmap)
        total_len = contig.recombination_map.sequence_length
        neu_prop = 0.3
        nonneu_prop = 0.7

        ts = tskit.load(input.ts)
        index = int(wildcards.ids)
        samps = ts.samples(population=index)
        if len(samps) == 0: samps = ts.samples(population=0)

        if (wildcards.annots != "all_sites") and (wildcards.annots != "none"):
            annot = species.get_annotations(wildcards.annots)
            annot_intervals = annot.get_chromosome_annotations(wildcards.chrms)
            exon_len = np.sum(annot_intervals[:,1]-annot_intervals[:,0])
            ts2fs.generate_fs(ts, samps, output, format='DFE-alpha', intervals=annot_intervals, 
                              is_folded=True, data_path_1=params.dfe_alpha_data_path_1, data_path_2=params.dfe_alpha_data_path_2,
                              sfs_input_file=output.fs, est_dfe_results_dir=params.output, est_dfe_demography_results_file=params.output+"/neu/est_dfe.out")
        else: 
            ts2fs.generate_fs(ts, samps, output, format='DFE-alpha', 
                              is_folded=True, data_path_1=params.dfe_alpha_data_path_1, data_path_2=params.dfe_alpha_data_path_2,
                              sfs_input_file=output.fs, est_dfe_results_dir=params.output, est_dfe_demography_results_file=params.output+"/neu/est_dfe.out")


rule run_dfe_alpha:
    input:
        neu_config = rules.generate_dfe_alpha_fs.output.neu_config,
        nonneu_config = rules.generate_dfe_alpha_fs.output.nonneu_config,
    output:
        res = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/nonneu/est_dfe.out",
    params:
        dfe_alpha_exec = os.path.abspath(config["dfe_alpha_exec"]),
        working_dir = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/", # Avoid input file name too long
        neu_config = "pop{ids}.DFE-alpha.neu.config",
        nonneu_config = "pop{ids}.DFE-alpha.nonneu.config",
    threads: 1,
    shell:
        """
        d=$PWD
        cd {params.working_dir}
        {params.dfe_alpha_exec} -c {params.neu_config} 
        sleep 60
        {params.dfe_alpha_exec} -c {params.nonneu_config}
        cd $d
        """


rule get_dfe_alpha_bestfits:
    input:
        res = rules.run_dfe_alpha.output.res,
    output:
        bestfit = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.DFE-alpha.bestfit",
    threads: 1,
    shell:
        """
        cat {input.res} | awk -v seeds={wildcards.seeds} -v genmap={wildcards.genmap} -v demog={wildcards.demog} -v dfe={wildcards.dfes} -v annot={wildcards.annots} -v chr={wildcards.chrms} -v id={wildcards.ids} 'BEGIN{{OFS="\\t"}}{{print seeds,genmap,demog,dfe,annot,chr,"pop"id,$16,$12,$10}}' > {output.bestfit}
        """


rule dfe_alpha_summary:
    input:
        bestfits = expand(output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.DFE-alpha.bestfit",
                          genmap=genmap_list, demog=demo_model_id_list, dfes=dfe_list, annots=annotation_list, seeds=seed_list, chrms=chrm_list, ids=[0,1,2]),
    output:
        summary = output_dir + "/inference/dfe/bestfits/DFE-alpha.bestfits.txt",
    threads: 1,
    shell:
        """
        cat {input.bestfits} | sed '1iseed\\tgenetic_map\\tdemography\\tDFE\\tannotation\\tchromosome\\tpopulation\\tlikelihood\\tEs\\tshape' > {output.summary} 
        """
