"""
polyDFE pipeline for running DFE benchmark on stdpopsim.
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
# polyDFE
# ###############################################################################


rule all:
    input:
        output_dir + "/inference/dfe/bestfits/polyDFE.bestfits.txt",


rule generate_polydfe_fs:
    input:
        ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.trees",
    output:
        fs = output_dir + "/inference/dfe/polyDFE/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.polyDFE.fs",
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
        if len(samps) > 20: samps = samps[:20]

        if (wildcards.annots != "all_sites") and (wildcards.annots != "none"):
            annot = species.get_annotations(wildcards.annots)
            annot_intervals = annot.get_chromosome_annotations(wildcards.chrms)
            exon_len = np.sum(annot_intervals[:,1]-annot_intervals[:,0])
            ts2fs.generate_fs(ts, samps, output.fs, format='polyDFE', intervals=annot_intervals, 
                              seq_len=exon_len, neu_prop=neu_prop, nonneu_prop=nonneu_prop, sample_size=len(samps))
        else: 
            ts2fs.generate_fs(ts, samps, output.fs, format='polyDFE', 
                              seq_len=total_len, neu_prop=neu_prop, nonneu_prop=nonneu_prop, sample_size=len(samps))
         

rule run_polydfe:
    input: 
        fs = rules.generate_polydfe_fs.output,
    output: 
        res = output_dir + "/inference/dfe/polyDFE/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.polyDFE.out",
    params:
        poly_dfe_exec = os.path.abspath(config["poly_dfe_exec"]),
    threads: 1,
    shell:
        """
        {params.poly_dfe_exec} -d {input.fs} -m C -i workflows/config/polyDFE/polyDFE_init_models.txt 1 -e > {output.res}
        """


rule get_polydfe_bestfits:
    input:
        polydfe_out = rules.run_polydfe.output,
    output:
        bestfit = output_dir + "/inference/dfe/polyDFE/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.polyDFE.bestfit",
    threads: 1,
    shell:
        """
        paste <(echo pop{wildcards.ids}) \
              <(grep "Best joint likelihood" {input.polydfe_out} | awk '{{print $6}}') \
              <(grep eps_an {input} -A 3 | tail -1 | awk 'BEGIN{{OFS="\\t"}}{{print $2,$3}}') \
              <(grep eps_an {input} -A 1 | tail -1 | awk '{{print $3}}') | awk -v seeds={wildcards.seeds} -v genmap={wildcards.genmap} -v demog={wildcards.demog} -v dfe={wildcards.dfes} -v annot={wildcards.annots} -v chr={wildcards.chrms} -v id={wildcards.ids} 'BEGIN{{OFS="\\t"}}{{print seeds,genmap,demog,dfe,annot,chr,$0}}' > {output.bestfit}
        """


rule polydfe_summary:
    input:
        bestfits = expand(output_dir + "/inference/dfe/polyDFE/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.polyDFE.bestfit",
                          genmap=genmap_list, demog=demo_model_id_list, dfes=dfe_list, annots=annotation_list, seeds=seed_list, chrms=chrm_list, ids=[0,1,2])
    output:
        summary = output_dir + "/inference/dfe/bestfits/polyDFE.bestfits.txt",
    threads: 1
    shell:
        """
        cat {input.bestfits} | sed '1iseed\\tgenetic_map\\tdemography\\tDFE\\tannotation\\tchromosome\\tpopulation\\tlikelihood\\tEs\\tshape\\ttheta_bar' > {output.summary}
        """
