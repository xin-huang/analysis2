"""
grapes pipeline for running DFE benchmark on stdpopsim.
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
# grapes
# ###############################################################################


rule all:
    input:
        output_dir + "/inference/dfe/bestfits/grapes.bestfits.txt",


rule generate_grapes_fs:
    input:
        ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.trees",
    output:
        fs = output_dir + "/inference/dfe/grapes/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.grapes.fs",
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

        header = species.common_name + " " + wildcards.chrms
        data_description = wildcards.annots

        if (wildcards.annots != "all_sites") and (wildcards.annots != "none"):
            annot = species.get_annotations(wildcards.annots)
            annot_intervals = annot.get_chromosome_annotations(wildcards.chrms)
            exon_len = np.sum(annot_intervals[:,1]-annot_intervals[:,0])
            ts2fs.generate_fs(ts, samps, output.fs, format='grapes', intervals=annot_intervals,
                              header=header, data_description=data_description,
                              seq_len=exon_len, neu_prop=neu_prop, nonneu_prop=nonneu_prop, sample_size=len(samps))
        else:
            ts2fs.generate_fs(ts, samps, output.fs, format='grapes',
                              header=header, data_description=data_description,
                              seq_len=total_len, neu_prop=neu_prop, nonneu_prop=nonneu_prop, sample_size=len(samps))


rule run_grapes:
    input:
        fs = rules.generate_grapes_fs.output,
    output:
        ept_fs = output_dir + "/inference/dfe/grapes/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.grapes.ept.fs",
        res = output_dir + "/inference/dfe/grapes/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.grapes.out",
    params:
        grapes_exec = os.path.abspath(config["grapes_exec"]),
    threads: 1,
    shell:
        """
        {params.grapes_exec} -in {input.fs} -out {output.ept_fs} -model GammaZero -no_div_param > {output.res}
        """


rule get_grapes_bestfits:
    input:
        grapes_out = output_dir + "/inference/dfe/grapes/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.grapes.out",
    output:
        bestfit = output_dir + "/inference/dfe/grapes/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.grapes.bestfit",
    threads: 1,
    shell:
        """
        grep Separate -A 4 {input.grapes_out} | sed '1d' | awk -F ": " '{{print $2}}' | awk -v seeds={wildcards.seeds} -v genmap={wildcards.genmap} -v demog={wildcards.demog} -v dfe={wildcards.dfes} -v annot={wildcards.annots} -v chr={wildcards.chrms} -v id={wildcards.ids} 'BEGIN{{RS="\\t";FS="\\n";OFS="\\t"}}{{print seeds,genmap,demog,dfe,annot,chr,"pop"id,$1,$3,$2,$4}}' > {output.bestfit}
        """


rule grapes_summary:
    input:
        bestfits = expand(output_dir + "/inference/dfe/grapes/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.grapes.bestfit",
                          genmap=genmap_list, demog=demo_model_id_list, dfes=dfe_list, annots=annotation_list, seeds=seed_list, chrms=chrm_list, ids=[0,1,2]),
    output:
        summary = output_dir + "/inference/dfe/bestfits/grapes.bestfits.txt",
    threads: 1,
    shell:
        """
        cat {input.bestfits} | sed '1iseed\\tgenetic_map\\tdemography\\tDFE\\tannotation\\tchromosome\\tpopulation\\tlikelihood\\tEs\\tshape\\ttheta' > {output.summary} 
        """
