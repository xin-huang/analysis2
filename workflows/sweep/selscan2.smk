chr_name_list = [c+1 for c in range(22)]
population_list = ['A', 'B']
maf_list = [0.01, 0.05]
method_list = ['ihs', 'nsl']
flag_list = ['skip', 'keep']
threshold_list = [2, 3, 4]


rule all:
    input:
        expand("results/{population}/{method}_{maf}/{population}.chr{chr_name}.{flag}.low.freq.{method}.threshold_{threshold}.candidates",
               population=population_list, maf=maf_list, chr_name=chr_name_list, method=method_list, flag=flag_list, threshold=threshold_list)


rule extract_samples:
    input:
        vcf = "results/chr{chr_name}.biallelic.snps.vcf.gz",
        samples = "data/populations/{population}.list",
    output:
        vcf = "results/{population}/{population}.chr{chr_name}.vcf.gz",
        map = "results/{population}/{population}.chr{chr_name}.map",
    shell:
        """
        bcftools view {input.vcf} -S {input.samples} --force-samples | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        bcftools query -f "%CHROM\\t%CHROM:%POS:%REF:%ALT\\t%POS\\t%POS\\n" {output.vcf} > {output.map}
        """


rule estimate_stats:
    input:
        vcf = rules.extract_samples.output.vcf,
        map = rules.extract_samples.output.map,
    output:
        out = "results/{population}/{method}_{maf}/{population}.chr{chr_name}.{flag}.low.freq.{method}.out",
        log = "results/{population}/{method}_{maf}/{population}.chr{chr_name}.{flag}.low.freq.{method}.log",
    resources:
        cpus = 8,
    params:
        output_prefix = lambda wildcards: f"results/{wildcards.population}/{wildcards.method}_{wildcards.maf}/{wildcards.population}.chr{wildcards.chr_name}.{wildcards.flag}.low.freq",
    shell:
        """
        ext/selscan/bin/linux/selscan --vcf {input.vcf} --map {input.map} --{wildcards.method} --out {params.output_prefix} --threads {resources.cpus} --{wildcards.flag}-low-freq --maf {wildcards.maf}
        """


rule normalize_stats:
    input:
        files = expand("results/{population}/{method}_{maf}/{population}.chr{chr_name}.{flag}.low.freq.{method}.out", chr_name=chr_name_list, allow_missing=True),
    output:
        files = expand("results/{population}/{method}_{maf}/{population}.chr{chr_name}.{flag}.low.freq.{method}.out.100bins.norm", chr_name=chr_name_list, allow_missing=True),
        log = "results/{population}/{method}_{maf}/{population}.{flag}.low.freq.normalized.{method}.log",
    shell:
        """
        ext/selscan/bin/linux/norm --files {input.files} --log {output.log} --{wildcards.method}
        """


rule get_candidates:
    input:
        file = "results/{population}/{method}_{maf}/{population}.chr{chr_name}.{flag}.low.freq.{method}.out.100bins.norm",
    output:
        file = "results/{population}/{method}_{maf}/{population}.chr{chr_name}.{flag}.low.freq.{method}.threshold_{threshold}.candidates",
    shell:
        """
        awk '($7>{wildcards.threshold})||($7<-{wildcards.threshold})' {input.file} > {output.file}
        """
