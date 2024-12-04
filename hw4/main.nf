#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samples_file = "./samples.csv"
params.fastqc_out = "./test_output/fastqc"
params.fastp_out = "./test_output/fastp"
params.spades_out = "./test_output/spades"
params.quast_out = "./test_output/quast"
params.prokka_out = "./test_output/prokka"
params.abricate_out = "./test_output/abricate"
params.prefix = "SRR31122807"
params.prefix_trimmed = "SRR31122807_trim"

// Load sample CSV
Channel
    .fromPath(params.samples_file)
    .splitCsv(header: true)
    .map { row -> 
        [
            sample_id: row.sample_id,
            read_1: file(row.read_1), 
            read_2: file(row.read_2)
        ] 
    }
    .set { samples }

// Workflow
workflow {
    samples | each { sample ->
        process fastqc_raw {
            input:
            path read_1 = sample.read_1
            path read_2 = sample.read_2
            
            output:
            path "${params.fastqc_out}/${sample.sample_id}_1_fastqc.html"
            path "${params.fastqc_out}/${sample.sample_id}_1_fastqc.zip"
            path "${params.fastqc_out}/${sample.sample_id}_2_fastqc.html"
            path "${params.fastqc_out}/${sample.sample_id}_2_fastqc.zip"

            script:
            """
            fastqc ${read_1} ${read_2} -o ${params.fastqc_out}
            """
        }

        process fastqc_trimming {
            input:
            path read_1 = sample.read_1
            path read_2 = sample.read_2
            
            output:
            path "${params.fastp_out}/${params.prefix_trimmed}_1.fastq.gz"
            path "${params.fastp_out}/${params.prefix_trimmed}_2.fastq.gz"

            script:
            """
            fastp -i ${read_1} -I ${read_2} -o ${output} -O ${output}
            """
        }

        process fastqc_trimmed {
            input:
            path trimmed_1 = fastqc_trimming.out[0]
            path trimmed_2 = fastqc_trimming.out[1]

            output:
            path "${params.fastqc_out}/${params.prefix_trimmed}_1_fastqc.html"
            path "${params.fastqc_out}/${params.prefix_trimmed}_1_fastqc.zip"
            path "${params.fastqc_out}/${params.prefix_trimmed}_2_fastqc.html"
            path "${params.fastqc_out}/${params.prefix_trimmed}_2_fastqc.zip"

            script:
            """
            fastqc ${trimmed_1} ${trimmed_2} -o ${params.fastqc_out}
            """
        }

        process spades {
            input:
            path trimmed_1 = fastqc_trimming.out[0]
            path trimmed_2 = fastqc_trimming.out[1]

            output:
            path "${params.spades_out}/contigs.fasta"
            path "${params.spades_out}/scaffolds.fasta"

            script:
            """
            spades.py -1 ${trimmed_1} -2 ${trimmed_2} -o ${params.spades_out}
            """
        }

        process quast {
            input:
            path contigs = spades.out[0]
            path scaffolds = spades.out[1]

            output:
            path "${params.quast_out}/report.html"
            path "${params.quast_out}/icarus.html"

            script:
            """
            quast.py -r ${contigs} ${scaffolds} -o ${params.quast_out}
            """
        }

        process prokka {
            input:
            path contigs = spades.out[0]
            path scaffolds = spades.out[1]

            output:
            path "${params.prokka_out}/PROKKA_12042024.gbf"
            path "${params.prokka_out}/PROKKA_12042024.gff"
            path "${params.prokka_out}/PROKKA_12042024.fsa"
            path "${params.prokka_out}/PROKKA_12042024.fna"

            script:
            """
            prokka --kingdom Viruses ${contigs} ${scaffolds} -o ${params.prokka_out} --force
            """
        }

        process abricate {
            input:
            path prokka_gbf = prokka.out[0]
            path prokka_gff = prokka.out[1]
            path prokka_fsa = prokka.out[2]
            path prokka_fna = prokka.out[3]

            output:
            path "${params.abricate_out}/results.tab"

            script:
            """
            abricate --db ncbi --db resfinder --quiet ${prokka_gbf} ${prokka_gff} ${prokka_fsa} ${prokka_fna} > ${params.abricate_out}/results.tab
            """
        }
    }
}