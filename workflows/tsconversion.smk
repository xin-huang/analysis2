"""
Snakefile for converting simulated tree-sequences from stdpopsim to different formats.
"""

import os
import numpy as np
import stdpopsim
import tskit
import masks
import ts2fs

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
genmap_list = [map for map in config['genetic_map']]


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


# Select annotation file from stdpopsim
annotation_list = config["annotation_list"]


# ###############################################################################
# GENERAL RULES & GLOBALS
# ###############################################################################


rule all:
    input:
        expand(output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.{extensions}",
            seeds=seed_array,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            chrms=chrm_list,
            ids=[0,1,2],
            extensions=['vcf.gz', 'ped', 'map'],
        ),
        expand(output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.{types}.fs",
            seeds=seed_array,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            chrms=chrm_list,
            ids=[0,1,2],
            types=['folded', 'unfolded'],
        ),
        expand(output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2{extensions}",
            seeds=seed_array,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            chrms=chrm_list,
            ids=[0,1,2],
            extensions=['_DAFpop0.obs', '.tpl', '.est'],
        ),
        expand(output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.{types}.blueprint",
            seeds=seed_array,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            chrms=chrm_list,
            ids=[0,1,2],
            types=['folded', 'unfolded'],
        ),
        expand(output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.{types}.fs",
            seeds=seed_array,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            chrms=chrm_list,
            ids=[0,1,2],
            types=['neu', 'nonneu'],
        ),
        expand(output_dir + "/inference/dfe/{tools}/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.{tools}.fs",
            seeds=seed_array,
            genmap=genmap_list,
            demog=demo_model_ids,
            dfes=dfe_list,
            annots=annotation_list,
            chrms=chrm_list,
            ids=[0,1,2],
            tools=['polyDFE', 'DFE-alpha', 'grapes'],
        ),


rule ts2vcf:
    input:
        ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.trees",
        mask_file = "workflows/masks/HapmapII_GRCh37.mask.bed",
    output:
        masked_ts = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.masked.trees",
        masked_vcf = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/sim_{chrms}.masked.vcf",
    run:
        ts = tskit.load(input.ts)
        if wildcards.annots == 'all_sites':
            mask_intervals = masks.get_combined_masks(species.id, input.mask_file, wildcards.chrms)
        else:
            mask_intervals = masks.get_combined_masks(species.id, input.mask_file, wildcards.chrms, chrom_annotation=wildcards.annots)

        masked_ts = ts.delete_intervals(mask_intervals)
        masked_ts.dump(output.masked_ts)

        with open(output.masked_vcf, 'w') as o:
            masked_ts.write_vcf(o, contig_id=wildcards.chrms)


rule get_ind_list:
    input:
    output:
        ind = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.ind.list",
    run:
        if wildcards.demog == 'Constant':
            samples = [i for i in range(int(demo_sample_size_dict[wildcards.demog][0]/2))]
        else:
            samples = [i+int(demo_ind_start_index_dict[wildcards.demog][int(wildcards.ids)]) for i in range(int(demo_sample_size_dict[wildcards.demog][int(wildcards.ids)]/2))]

        with open(output.ind, 'w') as o:
            for s in samples:
                o.write(f'tsk_{s}\n')


rule extract_biallelic_variants:
    input:
        vcf = rules.ts2vcf.output.masked_vcf,
        ind = rules.get_ind_list.output.ind,
    output:
        vcf = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.vcf.gz",
    shell:
        """
        bgzip -c {input.vcf} | bcftools view -m 2 -M 2 -S {input.ind} | bcftools view -i 'INFO/AC>0&INFO/AC<INFO/AN' | awk -F "\\t" 'BEGIN{{OFS="\\t";ORS=""}}{{if($0~/^#/){{print $0"\\n"}}else{{print $1,$2,$3,"A","T",".","PASS",$8";AA=A\\t";for(i=9;i<NF;i++){{print $i"\\t"}};print $NF"\\n"}}}}' | sed '/##FORMAT/a##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">"' | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule vcf2plink:
    input:
        vcf = rules.extract_biallelic_variants.output.vcf,
    output:
        ped = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.ped",
        map = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked.map",
    params:
        output_prefix = output_dir + "/simulated_data/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.biallelic.masked",
    shell:
        """
        plink --vcf {input.vcf} --const-fid pop{wildcards.ids} --recode tab 12 --geno 0.0 --out {params.output_prefix} --set-missing-var-ids @:#
        """


rule ts2fs_demog: # AFS for demographic inference
    input:
        ts = rules.ts2vcf.output.masked_ts,
    output:
        dadi_folded_fs = output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.folded.fs",
        dadi_unfolded_fs = output_dir + "/inference/demog/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.unfolded.fs",
        fastsimcoal2_dsfs = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2_DAFpop0.obs",
        fastsimcoal2_tpl = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.tpl",
        fastsimcoal2_est = output_dir + "/inference/demog/fastsimcoal2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.fastsimcoal2.est",
        stairwayplot2_folded_fs = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.folded.blueprint",
        stairwayplot2_unfolded_fs = output_dir + "/inference/demog/stairwayplot2/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.stairwayplot2.unfolded.blueprint",
    params:
        popid = "pop{ids}", # id of the population
        nseq = 100, # number of sequences
        L = lambda wildcards: species.get_contig(wildcards.chrms, genetic_map=wildcards.genmap).length, # total number of observed nucleic sites, including polymorphic and monomorphic
        pct_training = 0.67, # percentage of sites for training
        nrand = "7 15 22 28", # number of random break points for each try (separated by white space)
        project_dir = "demography", # project directory
        stairway_plot_dir = "stairway_plot_es", # directory to the stairway plot files
        ninput = 200, # number of input files to be created for each estimation
        mu = lambda wildcards: species.genome.mean_mutation_rate if wildcards.demog == 'Constant' else species.get_demographic_model(wildcards.demog).mutation_rate, # assumed mutation rate per site per generation
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

        # dadi
        ts2fs.generate_fs(ts, samps, [output.dadi_folded_fs], format='dadi', is_folded=True, is_separated=False)
        ts2fs.generate_fs(ts, samps, [output.dadi_unfolded_fs], format='dadi', is_folded=False, is_separated=False)

        # fastsimcoal2
        ts2fs.generate_fs(ts, samps, output.fastsimcoal2_dsfs, format='fastsimcoal2')
        with open(output.fastsimcoal2_tpl, 'w') as o:
            o.write("//Number of population samples (demes)\n")
            o.write("1\n")
            o.write("//Population effective sizes (number of genes)\n")
            o.write("NPOP\n")
            o.write("//Sample sizes\n")
            o.write(f"{params.nseq}\n")
            o.write("//Growth rates : negative growth implies population expansion\n")
            o.write("0\n")
            o.write("//Number of migration matrices : 0 implies no migration between demes\n")
            o.write("0\n")
            o.write("//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix\n")
            o.write("1 historical event\n")
            o.write("TEXP 0 0 0 RESIZE 0 0\n")
            o.write("//Number of independent loci [chromosome]\n")
            o.write("1 0\n")
            o.write("//Per chromosome: Number of linkage blocks\n")
            o.write("1\n")
            o.write("//per Block: data type, num loci, rec. rate and mut rate + optional parameters\n")
            o.write(f"FREQ 1 0 {params.mu}\n")
    
        with open(output.fastsimcoal2_est, 'w') as o:
            o.write("// Priors and rules file\n")
            o.write("// *********************\n")
            o.write("[PARAMETERS]\n")
            o.write("//#isInt? #name #dist. #min #max\n")
            o.write("//all N are in number of haploid individuals\n")
            o.write("1 NPOP logunif 100 100000 output bounded\n")
            o.write("1 TEXP logunif 100 5000 output\n")
            o.write("0 RESIZE logunif 1e-3 1000 output\n")
            o.write("[COMPLEX PARAMETERS]\n")
            o.write("1 ANCSIZE = NPOP*RESIZE output\n")

        # stairwayplot2
        ts2fs.generate_fs(ts, samps, output.stairwayplot2_folded_fs, format='stairwayplot2', is_folded=True,
                          popid=params.popid, nseq=params.nseq, L=params.L, pct_training=params.pct_training,
                          nrand=params.nrand, project_dir=params.project_dir, stairway_plot_dir=params.stairway_plot_dir,
                          ninput=params.ninput, mu=params.mu, year_per_generation=params.year_per_generation,
                          plot_title=params.plot_title, xrange=params.xrange, yrange=params.yrange,
                          xspacing=params.xspacing, yspacing=params.yspacing, fontsize=params.fontsize)
        ts2fs.generate_fs(ts, samps, output.stairwayplot2_unfolded_fs, format='stairwayplot2', is_folded=False,
                          popid=params.popid, nseq=params.nseq, L=params.L, pct_training=params.pct_training,
                          nrand=params.nrand, project_dir=params.project_dir, stairway_plot_dir=params.stairway_plot_dir,
                          ninput=params.ninput, mu=params.mu, year_per_generation=params.year_per_generation,
                          plot_title=params.plot_title, xrange=params.xrange, yrange=params.yrange,
                          xspacing=params.xspacing, yspacing=params.yspacing, fontsize=params.fontsize)


rule ts2fs_dfe: # AFS for DFE inference
    input:
        ts = rules.ts2vcf.output.masked_ts,
    output:
        dadi_neu_fs = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.neu.fs",
        dadi_nonneu_fs = output_dir + "/inference/dfe/dadi/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.dadi.nonneu.fs",
        polydfe_fs = output_dir + "/inference/dfe/polyDFE/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.polyDFE.fs",
        dfe_alpha_fs = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.DFE-alpha.fs",
        dfe_alpha_neu_config = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.DFE-alpha.neu.config",
        dfe_alpha_nonneu_config = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.DFE-alpha.nonneu.config",
        grapes_fs = output_dir + "/inference/dfe/grapes/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}/sim_{chrms}.pop{ids}.grapes.fs",
    params:
        dfe_alpha_output = output_dir + "/inference/dfe/DFE-alpha/{genmap}/{demog}/{dfes}/{annots}/{seeds}/pop{ids}",
        dfe_alpha_data_path_1 = os.path.abspath(config["dfe_alpha_data_path_1"]), # Path for the data files for the one and two epoch models from DFE-alpha
        dfe_alpha_data_path_2 = os.path.abspath(config["dfe_alpha_data_path_2"]), # Path for the data files for the three epoch model from DFE-alpha
    run:
        ts = tskit.load(input.ts)
        index = int(wildcards.ids)
        samps = ts.samples(population=index)
        if len(samps) == 0: samps = ts.samples(population=0)

        if wildcards.genmap == 'Uniform_GeneticMap': contig = species.get_contig(wildcards.chrms)
        else: contig = species.get_contig(wildcards.chrms, genetic_map=wildcards.genmap)
        total_len = contig.recombination_map.sequence_length
        neu_prop = 0.3
        nonneu_prop = 0.7

        # dadi 
        ts2fs.generate_fs(ts, samps, [output.dadi_neu_fs, output.dadi_nonneu_fs], format='dadi')

        # DFE-alpha
        ts2fs.generate_fs(ts, samps, [output.dfe_alpha_neu_config, output.dfe_alpha_nonneu_config, output.dfe_alpha_fs], format='DFE-alpha', 
                          is_folded=True, data_path_1=params.dfe_alpha_data_path_1, data_path_2=params.dfe_alpha_data_path_2,
                          sfs_input_file=output.dfe_alpha_fs, est_dfe_results_dir=params.dfe_alpha_output, est_dfe_demography_results_file=params.dfe_alpha_output+"/neu/est_dfe.out")

        # grapes
        header = species.common_name + " " + wildcards.chrms
        data_description = wildcards.annots
        ts2fs.generate_fs(ts, samps, output.grapes_fs, format='grapes', 
                          header=header, data_description=data_description,
                          seq_len=total_len, neu_prop=neu_prop, nonneu_prop=nonneu_prop, sample_size=len(samps))

        # polyDFE
        if len(samps) > 20: samps = samps[:20]
        ts2fs.generate_fs(ts, samps, output.polydfe_fs, format='polyDFE', 
                          seq_len=total_len, neu_prop=neu_prop, nonneu_prop=nonneu_prop, sample_size=len(samps))