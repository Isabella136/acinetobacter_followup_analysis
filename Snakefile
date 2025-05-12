wildcard_constraints:
    biosample = r"SAM[EAN]+\d+",
    run = r"SRR\d+"

configfile: "/".join((workflow.basedir, "config/config.json"))

with open("/".join((workflow.basedir, config["accession_list"])), 'r') as accession_file:
    run_accessions = accession_file.read().split('\n')

with open("/".join((workflow.basedir, config["biosample_list"])), 'r') as biosample_file:
    biosamples = biosample_file.read().split('\n')

rule all:
    input:
        # Reference sequence from NCBI
        "reference.fasta",

        # Reads
        expand(
            "biosample/{biosample}/reads/{run}_1.fastq",
            zip,
            biosample = biosamples,
            run = run_accessions
        ),

        # mlst
        "/".join((config["output_dir"], "abaumannii_mlst.tsv")),

rule make_sample_dir:
    output: temp("{biosample}/text.txt")

    shell: "mkdir -p {wildcards.biosample} | cd {wildcards.biosample} | touch {output}"

rule get_ref:
    output: "reference.fasta"

    conda: "envs/ncbi-downloads.yaml"

    shell: """
        datasets download genome taxon 696749 \
        --filename reference.zip

        unzip -q reference.zip -d reference
        rm reference.zip

        mv reference/ncbi_dataset/data/* reference/
        mv reference/*/*.fna reference/
        rm -r reference/*/

        cat reference/*.fna > {output}
        rm -r reference/
        """

rule get_runs:
    input: rules.make_sample_dir.output

    output:
        fastq_1 = "biosample/{biosample}/reads/{run}_1.fastq",
        fastq_2 = "biosample/{biosample}/reads/{run}_2.fastq"

    conda: "envs/ncbi-downloads.yaml"

    shell: """
        cd {wildcards.biosample}/
        mkdir -p reads
        cd reads
        fasterq-dump --split-files {wildcards.run}
        """   
    
rule spades_assembly:
    input: rules.get_runs.output

    params: "/".join((config["output_dir"], "{biosample}/{run}"))

    output: "/".join((config["output_dir"], "{biosample}/{run}/contigs.fasta"))

    conda: "envs/assembly.yaml"

    resources:
        slurm_extra = "--qos=high"

    shell: "spades.py -1 {input[0]} -2 {input[1]} -o {params} -t {threads}"

rule mlst_analysis:
    input: 
        expand("/fs/cbcb-lab/mpop/projects/Acinetobacter_followup/{bio}/{run}/contigs.fasta",
            zip,
            run=run_accessions,
            bio=biosamples)

    params: config["output_dir"]

    output: "/".join((config["output_dir"], "abaumannii_mlst.tsv"))

    conda: "envs/mlst.yaml"

    shell: """
        cd {params}
        mlst */contigs.fasta > {output}
        """