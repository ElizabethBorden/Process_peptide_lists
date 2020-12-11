# Setting up filenames here:
from os.path import join
configfile: "carreno.config.json"

# Files
sample = config["all_samples"]

# Path to files
peptide_path = "/scratch/eknodel/carreno_data/"

rule all:
    input:
        expand(peptide_path + "{sample}/{sample}_9_netctl.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_10_netctl.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_netCTL.out", sample=sample),
        expand(peptide_path + "{sample}/{sample}_9_netmhc.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_10_netmhc.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_netMHC.out", sample=sample),
        expand(peptide_path + "{sample}/{sample}_9_netmhcstab.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_10_netmhcstab.xsl", sample=sample),
        expand(peptide_path + "{sample}/{sample}_netMHCstab.out", sample=sample),
        expand(peptide_path + "{sample}/{sample}.epitopes.annotated.tsv", sample=sample),
        expand(peptide_path + "{sample}/{sample}.WTmatch.in", sample=sample),
        expand(peptide_path + "{sample}/{sample}.WTmatch_9_netMHC.out", sample=sample),
        expand(peptide_path + "{sample}/{sample}.WTmatch_10_netMHC.out", sample=sample),
        expand(peptide_path + "{sample}/{sample}.WTmatch_netMHC.out", sample=sample)

rule netCTL:
    input:
        peptides = os.path.join(peptide_path, "{sample}/{sample}_peptides.in")
    output:
        netCTL_9  = os.path.join(peptide_path, "{sample}/{sample}_9_netctl.xsl"),
        netCTL_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netctl.xsl"),
        netCTL_11 = os.path.join(peptide_path, "{sample}/{sample}_11_netctl.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        netCTLpan -a {params.hla} -f {input.peptides} -s -xls -l 9  -xlsfile {output.netCTL_9};
        netCTLpan -a {params.hla} -f {input.peptides} -s -xls -l 10 -xlsfile {output.netCTL_10};
        netCTLpan -a {params.hla} -f {input.peptides} -s -xls -l 11 -xlsfile {output.netCTL_11};
        """

rule combine_netCTL:
    input:
        netCTL_9 = os.path.join(peptide_path, "{sample}/{sample}_9_netctl.xsl"),
        netCTL_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netctl.xsl"),
        netCTL_11 = os.path.join(peptide_path, "{sample}/{sample}_11_netctl.xsl"),
    output:
        netCTL = os.path.join(peptide_path, "{sample}/{sample}_netCTL.out")
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        cat {input.netCTL_11} | sed '/Allele/d' | awk {{'print $2, "{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netCTL};
        cat {input.netCTL_10} | sed '/Allele/d' | awk {{'print $2, "{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netCTL}
        cat {input.netCTL_9} | sed '/Allele/d'  | awk {{'print $2, "{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netCTL} 
        """

rule netMHC:
    input:
        peptides = os.path.join(peptide_path, "{sample}/{sample}_peptides.in")
    output:
        netMHC_9  = os.path.join(peptide_path, "{sample}/{sample}_9_netmhc.xsl"),
        netMHC_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhc.xsl"),
        netMHC_11 = os.path.join(peptide_path, "{sample}/{sample}_11_netmhc.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 9  -xlsfile {output.netMHC_9};
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 10 -xlsfile {output.netMHC_10};
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 11 -xlsfile {output.netMHC_11};
        """

rule combine_netMHC:
    input:
        netMHC_9 = os.path.join(peptide_path, "{sample}/{sample}_9_netmhc.xsl"),
        netMHC_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhc.xsl"),
        netMHC_11 = os.path.join(peptide_path, "{sample}/{sample}_11_netmhc.xsl"),
    output:
        netMHC = os.path.join(peptide_path, "{sample}/{sample}_netMHC.out")
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        cat {input.netMHC_11} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2, "{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netMHC};
        cat {input.netMHC_10} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2, "{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netMHC}
        cat {input.netMHC_9} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2, "{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netMHC}
        """

rule netMHCstab:
    input:
        peptides = os.path.join(peptide_path, "{sample}/{sample}_peptides.in")
    output:
        netMHCstab_9  = os.path.join(peptide_path, "{sample}/{sample}_9_netmhcstab.xsl"),
        netMHCstab_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhcstab.xsl"),
        netMHCstab_11 = os.path.join(peptide_path, "{sample}/{sample}_11_netmhcstab.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 9  -xlsfile {output.netMHCstab_9};
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 10 -xlsfile {output.netMHCstab_10};
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 11 -xlsfile {output.netMHCstab_11};
        """

rule combine_netMHCstab:
    input:
        netMHCstab_9 = os.path.join(peptide_path, "{sample}/{sample}_9_netmhcstab.xsl"),
        netMHCstab_10 = os.path.join(peptide_path, "{sample}/{sample}_10_netmhcstab.xsl"),
        netMHCstab_11 = os.path.join(peptide_path, "{sample}/{sample}_11_netmhcstab.xsl"),
    output:
        netMHCstab = os.path.join(peptide_path, "{sample}/{sample}_netMHCstab.out")
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        cat {input.netMHCstab_11} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2,"{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netMHCstab};
        cat {input.netMHCstab_10} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2,"{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netMHCstab}
        cat {input.netMHCstab_9} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2,"{params.hla}", $5, $6, $7, $8, $9}}' >> {output.netMHCstab}
        """

rule filter_blast:
    input:
        os.path.join(peptide_path, "{sample}/{sample}_peptides.in")
    output:
        os.path.join(peptide_path, "{sample}/{sample}.epitopes.annotated.tsv")
    params:
        sample = "{sample}",
        outdir = os.path.join(peptide_path,"{sample}")
    shell:
        """
        python sequence_similarity.py -o {params.outdir} -s {params.sample} -i {input}
        """

rule prep_WT_netMHC:
    input:
        os.path.join(peptide_path, "{sample}/{sample}.epitopes.annotated.tsv")
    output:
        os.path.join(peptide_path, "{sample}/{sample}.WTmatch.in")
    shell:
        """
        cat {input} | awk '{{print ">"$1 "\\n" $2}}' > {output}
        """

rule run_WT_netMHC:
    input:
        peptides = os.path.join(peptide_path, "{sample}/{sample}.WTmatch.in")
    output:
        netMHC_9 = os.path.join(peptide_path, "{sample}/{sample}.WTmatch_9_netMHC.out"),
        netMHC_10 = os.path.join(peptide_path, "{sample}/{sample}.WTmatch_10_netMHC.out"),
        netMHC_11 = os.path.join(peptide_path, "{sample}/{sample}.WTmatch_11_netMHC.out"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -l 9 -xls -xlsfile {output.netMHC_9};
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -l 10 -xls -xlsfile {output.netMHC_10};
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -l 11 -xls -xlsfile {output.netMHC_11};
        """

rule combine_outputs:
    input:
        netMHC_9 = os.path.join(peptide_path, "{sample}/{sample}.WTmatch_9_netMHC.out"),
        netMHC_10 = os.path.join(peptide_path, "{sample}/{sample}.WTmatch_10_netMHC.out"),
        netMHC_11 = os.path.join(peptide_path, "{sample}/{sample}.WTmatch_11_netMHC.out"),
    output:
        netMHC = os.path.join(peptide_path, "{sample}/{sample}.WTmatch_netMHC.out")
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        cat {input.netMHC_11} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2, "{params.hla}", $3, $5, $6, $7, $8, $9}}' >> {output.netMHC};
        cat {input.netMHC_10} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2, "{params.hla}", $3, $5, $6, $7, $8, $9}}' >> {output.netMHC};
        cat {input.netMHC_9} | sed '/Allele/d' | sed '/HLA/d' | awk {{'print $2, "{params.hla}", $3, $5, $6, $7, $8, $9}}' >> {output.netMHC}
        """
