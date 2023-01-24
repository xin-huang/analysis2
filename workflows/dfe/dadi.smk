"""
dadi pipeline for running DFE benchmark on stdpopsim.
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
# dadi
# ###############################################################################


dadi_dfe_list = ['gamma']
dadi_dfe_params_list = {
    'gamma': {
        'p0': '5 500',
        'ubounds': '100 5000',
        'lbounds': '10e-3 10e-3',
    },
    'lognormal': {
        'p0': '5 5',
        'ubounds': '100 100',
        'lbounds': '10e-3 10e-3',
    },
}


rule all:
    input:
        expand(output_dir + "/inference/dfe/bestfits/dadi.{dadi_dfe}.bestfits.txt", dadi_dfe=dadi_dfe_list)


rule generate_dadi_fs:
    input:
        ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.trees"
    output:
        neu_fs = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.dadi.neu.fs",
        nonneu_fs = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.dadi.nonneu.fs",
    threads: 1
    resources: time_min=60, mem_mb=5000, cpus=1
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
            ts2fs.generate_fs(ts, samps, [output.neu_fs, output.nonneu_fs], intervals=annot_intervals, format='dadi')
        else:
            ts2fs.generate_fs(ts, samps, [output.neu_fs, output.nonneu_fs], format='dadi')


rule dadi_infer_dm:
    input: 
        fs = rules.generate_dadi_fs.output.neu_fs,
    output:
        bestfit = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.two_epoch.InferDM.bestfits",
    params:
        demog = 'two_epoch',
        demog_p0 = '5 5',
        demog_ubounds = '10 1',
        demog_lbounds = '10e-3 10e-3',
        grid_size = '300 400 500',
        prefix = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.two_epoch",
    resources: time_min=3000, mem_mb=5000, cpus=8
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {params.demog} --p0 {params.demog_p0} --ubounds {params.demog_ubounds} --lbounds {params.demog_lbounds} --output-prefix {params.prefix} --grids {params.grid_size} --cpus {resources.cpus} --force-convergence --nomisid
        """


rule dadi_generate_cache:
    input:
        dm_bestfit = rules.dadi_infer_dm.output.bestfit,
    output:
        cache = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.two_epoch.spectra.bpkl",
    params:
        demog = 'two_epoch_sel',
        sample_size = 100,
        grid_size = '800 1000 1200',
        gamma_pts = 2000,
    resources: time_min=3000, mem_mb=5000, cpus=8
    shell:
        """
        dadi-cli GenerateCache --model {params.demog} --demo-popt {input.dm_bestfit} --sample-size {params.sample_size} --output {output.cache} --cpus {resources.cpus} --grids {params.grid_size} --gamma-pts {params.gamma_pts}
        """


rule dadi_infer_dfe:
    input:
        fs = rules.generate_dadi_fs.output.nonneu_fs,
        cache = rules.dadi_generate_cache.output.cache,
        dm_bestfit = rules.dadi_infer_dm.output.bestfit,
    output:
        bestfit = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.two_epoch.{dadi_dfe}.InferDFE.bestfits",
    params:
        dfe = lambda wildcards: wildcards.dadi_dfe, 
        dfe_p0 = lambda wildcards: dadi_dfe_params_list[wildcards.dadi_dfe]['p0'],
        dfe_lbounds = lambda wildcards: dadi_dfe_params_list[wildcards.dadi_dfe]['lbounds'],
        dfe_ubounds = lambda wildcards: dadi_dfe_params_list[wildcards.dadi_dfe]['ubounds'],
        ratio = 2.31,
        prefix = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.two_epoch.{dadi_dfe}",
    resources: time_min=3000, mem_mb=5000, cpus=8
    shell:
        """
        dadi-cli InferDFE --fs {input.fs} --cache1d {input.cache} --demo-popt {input.dm_bestfit} --output-prefix {params.prefix} --pdf1d {params.dfe} --p0 {params.dfe_p0} --ubounds {params.dfe_ubounds} --lbounds {params.dfe_lbounds} --ratio {params.ratio} --cpus {resources.cpus} --force-convergence --nomisid
        """


rule get_dadi_dfe_bestfits:
    input:
        bestfit = rules.dadi_infer_dfe.output.bestfit,
    output:
        res = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.dadi.{dadi_dfe}.bestfit",
    threads: 1
    resources: time_min=60, mem_mb=5000, cpus=1
    shell:
        """
        grep 'Converged results' {input.bestfit} -A 2 | tail -1 | awk -v seeds={wildcards.seeds} -v genmap={wildcards.genmap} -v demog={wildcards.demog} -v dfe={wildcards.dfes} -v annot={wildcards.annots} -v chr={wildcards.chrms} -v id={wildcards.ids} 'BEGIN{{OFS="\\t"}}{{print seeds,genmap,demog,dfe,annot,chr,"pop"id,$0}}' > {output}
        """


rule dadi_summary:
    input:
        gamma_bestfits = expand(output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/{chrms}/pop{ids}/pop{ids}.dadi.gamma.bestfit",
                                genmap=genmap_list, demog=demo_model_id_list, dfes=dfe_list, annots=annotation_list, seeds=seed_list, chrms=chrm_list, ids=[0,1,2]),
    output:
        gamma_summary = output_dir + "/inference/dfe/bestfits/dadi.gamma.bestfits.txt",
    threads: 1,
    shell:
        """
        cat {input.gamma_bestfits} | sed '1iseed\\tgenetic_map\\tdemography\\tDFE\\tannotation\\tchromosome\\tpopulation\\tlikelihood\\tshape\\tscale\\ttheta' > {output.gamma_summary} 
        """