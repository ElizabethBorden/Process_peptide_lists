# Setting up filenames here:
from os.path import join
configfile: "carreno.config.json"

# Files
sample = "All"

# Path to files
path = "/scratch/eknodel/carreno_data/"

rule all:
    input:
        expand(path + "{sample}_lukszaprogram.in", sample=sample),
        expand(path + "{sample}_lukszablast.in", sample=sample, path=path),
        expand(path + "{sample}_lukszablast.out", sample=sample, path=path),
        expand(path + "neoantigen_fitness_{sample}.txt", sample=sample, path=path),
        expand(path + "Carreno_{sample}_final_terms.in", sample=sample, path=path)

rule prep_luksza_input:
    input:
        netMHC = os.path.join(path, "{sample}.WTmatch_netMHC.out")
    output:
        luksza = os.path.join(path, "{sample}_lukszaprogram.in")
    params:
        directory = "/scratch/eknodel/carreno_data/",
        sample = sample
    shell:
        """
        Rscript Prep_luksza_input.R {params.directory} {params.sample}
        """

rule format_blast_input:
    input:
        os.path.join(path, "{sample}_lukszaprogram.in")
    output:
        os.path.join(path, "{sample}_lukszablast.in")
    shell:
        """
        cat {input} | awk '{{print ">sample|"$1"|MUT|"$2 "\\n" $5 "\\n" ">sample|"$1"|WT|"$2 "\\n" $4}}' > {output};        sed -i -e 1,4d {output}
        """

rule blast:
    input:
        os.path.join(path, "{sample}_lukszablast.in")
    output:
        os.path.join(path, "{sample}_lukszablast.out")
    shell:
        """
        blastp -query {input} -db ~/Cohen_melanoma/luksza_program/iedb.fasta -outfmt 5 -evalue 100000000  -gapopen 11 -gapextend 1 > {output}
        """
        #blastp -query {input} -db ~/Luksza_programs/Input/iedb.fasta -outfmt 5 -evalue 100000000  -gapopen 11 -gapextend 1 > {output} # with restricted list

rule luksza:
    input:
        blast = os.path.join(path, "{sample}_lukszablast.out"),
        neoantigens = os.path.join(path, "{sample}_lukszaprogram.in")
    output:
        os.path.join(path, "neoantigen_fitness_{sample}.txt")
    params:
        directory = "/scratch/eknodel/consortium_data/"
    shell:
        """
        cd ~/Luksza_programs;
        python src/main.py {input.neoantigens} {params.directory} 26 4.86936 {output} {input.blast}
        """

rule prep_final_file:
    input:
        netMHC = os.path.join(path, "{sample}.WTmatch_netMHC.out")
    output:
        luksza = os.path.join(path, "Carreno_{sample}_final_terms.in")
    params:
        directory = "/scratch/eknodel/carreno_data/",
        sample = sample
    shell:
        """
        Rscript Prep_final_file.R {params.directory} {params.sample}
        """
