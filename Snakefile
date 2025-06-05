wildcard_constraints:
    biosample = r"SAM[EAN]+\d+",
    run = r"SRR\d+"

configfile: "/".join((workflow.basedir, "config/config.json"))

with open("/".join((workflow.basedir, config["accession_list"])), 'r') as accession_file:
    run_accessions = accession_file.read().split('\n')

with open("/".join((workflow.basedir, config["biosample_list"])), 'r') as biosample_file:
    biosamples = biosample_file.read().split('\n')

with open("/".join((workflow.basedir, config["representative_accession_list"])), 'r') as accession_file:
    representative_run_accessions = accession_file.read().split('\n')

with open("/".join((workflow.basedir, config["representative_biosample_list"])), 'r') as biosample_file:
    representative_biosamples = biosample_file.read().split('\n')

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

        # Alignments
        expand(
            "/".join((config["output_dir"], "{biosample}/{run}/card_blast")),
            zip,
            biosample = biosamples,
            run = run_accessions
        ),

        # Alignments
        expand(
            "/".join((config["output_dir"], "{biosample}/{run}/wallace_blast")),
            zip,
            biosample = biosamples,
            run = run_accessions
        ),

        # mlst
        "/".join((config["output_dir"], "abaumannii_mlst.tsv")),

        # roary
        "/".join((config["output_dir"], "roary/summary_statistics.txt")),

        # representative roary
        "/".join((config["output_dir"], "representative_roary/summary_statistics.txt")),

rule make_sample_dir:
    output: temp("biosample/{biosample}/text.txt")

    shell: "mkdir -p biosample/{wildcards.biosample} | touch {output}"

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

rule quality_check:
    input: rules.get_runs.output

    output:
        fastq_1 = "biosample/{biosample}/reads/{run}_1_fastqc.html",
        fastq_2 = "biosample/{biosample}/reads/{run}_2_fastqc.html"

    shell: """
        module load fastqc
        fastqc biosample/{wildcards.biosample}/reads/*.fastq
    """
    
rule spades_assembly:
    input: rules.get_runs.output

    params: "/".join((config["output_dir"], "{biosample}/{run}/spades"))

    output: "/".join((config["output_dir"], "{biosample}/{run}/spades/contigs.fasta"))

    conda: "envs/assembly.yaml"

    resources:
        slurm_extra = "--qos=high"

    shell: """
        spades.py -1 {input[0]} -2 {input[1]} -o {params} -t {threads} || \
        spades.py -1 {input[0]} -2 {input[1]} -o {params} -t {threads} --phred-offset 33"""

rule mlst_analysis:
    input: 
        expand("../../{bio}/{run}/spades/contigs.fasta",
            zip,
            run=run_accessions,
            bio=biosamples)

    params: config["output_dir"]

    output: "/".join((config["output_dir"], "abaumannii_mlst.tsv"))

    conda: "envs/mlst.yaml"

    shell: """
        cd {params}
        mlst */*/contigs.fasta > {output}
        """

rule get_acb_complex_proteins:
    input: "/".join((workflow.basedir, "gene_id.txt"))

    output: "proteins.fasta"

    conda: "envs/ncbi-downloads.yaml"

    shell: """
        datasets download gene gene-id --inputfile {input} \
        --include protein --filename proteins.zip

        unzip -q proteins.zip -d proteins
        rm proteins.zip

        mv proteins/ncbi_dataset/data/* proteins/
        rm -r proteins/*/

        cat proteins/*.faa > {output}
        rm -r proteins/
        """

rule create_acb_complex_protein_database:
    input: rules.get_acb_complex_proteins.output

    output: "protein_database.fasta"

    conda: "envs/cd-hit.yaml"

    resources:
        slurm_extra = "--qos=high"

    shell: "cd-hit -i {input} -o {output} -s 0.8 -M 64000 -T {threads}"

rule sample_gene_prediction:
    input: 
        db = rules.create_acb_complex_protein_database.output,
        ct = rules.spades_assembly.output

    params: 
        outdir = "/".join((config["output_dir"], "{biosample}/{run}/prokka")),
        prefix = "{biosample}_{run}_results"

    output: "/".join((config["output_dir"], "{biosample}/{run}/prokka/{biosample}_{run}_results.gff"))

    conda: "envs/prokka.yaml"

    shell: "prokka --proteins {input.db} --outdir {params.outdir} --force\
        --prefix {params.prefix} {input.ct}"

rule reference_gene_prediction:
    input: 
        db = rules.create_acb_complex_protein_database.output,
        rf = rules.get_ref.output

    params: 
        outdir = "/".join((config["output_dir"], "reference/CP001921/prokka")),
        prefix = "ref_results"

    output: "/".join((config["output_dir"], "reference/CP001921/prokka/ref_results.gff"))

    conda: "envs/prokka.yaml"

    shell: "prokka --proteins {input.db} --outdir {params.outdir} --force\
        --prefix {params.prefix} {input.rf}"

rule pangenome_estimation:
    input:
        expand("../../{bio}/{run}/prokka/{bio}_{run}_results.gff",
            zip,
            run=run_accessions,
            bio=biosamples),
        rules.reference_gene_prediction.output

    params: config["output_dir"]

    output: "/".join((config["output_dir"], "roary/summary_statistics.txt"))

    conda: "envs/roary.yaml"

    resources:
        slurm_extra = "--qos=highmem"

    shell: """
        cd {params}
        rm -r roary || true
        roary -p {threads} -f roary -i 85 -s -v -e --mafft -n */*/prokka/*_results.gff
        mv roary_* roary || true
        """

rule pangenome_estimation_represeantative_only:
    input:
        expand("../../{bio}/{run}/prokka/{bio}_{run}_results.gff",
            zip,
            run=representative_run_accessions,
            bio=representative_biosamples),
        rules.reference_gene_prediction.output

    params: config["output_dir"]

    output: "/".join((config["output_dir"], "representative_roary/summary_statistics.txt"))

    conda: "envs/roary.yaml"

    resources:
        slurm_extra = "--qos=highmem"

    shell: """
        cd {params}
        rm -r representative_roary || true
        roary -p {threads} -f representative_roary -i 85 -s -v -e --mafft -n {input}
        mv representative_roary_* representative_roary || true
        """

rule blast_to_card:
    input: rules.spades_assembly.output

    params: "../CARD/card_homolog_and_knockout_models.fasta"

    output: "/".join((config["output_dir"], "{biosample}/{run}/card_blast"))

    conda: "envs/blast.yaml"

    shell: "blastn -query {params} -subject {input} -perc_identity 80.0 -outfmt 6 > {output}"

rule blast_to_resistant_wallace:
    input: rules.spades_assembly.output

    params: "../Wallace/wallace_genes_exclusive_resistant.fasta"

    output: "/".join((config["output_dir"], "{biosample}/{run}/wallace_blast"))

    conda: "envs/blast.yaml"

    shell: "blastn -query {params} -subject {input} -perc_identity 80.0 -outfmt 6 > {output}"