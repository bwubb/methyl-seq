process MULTIQC {

    tag { task.process.tokenize(':')[-1] }

    publishDir path:{ "${params.outdir}/${task.process.tokenize(':')[-1].toLowerCase()}" },mode:'copy',overwrite:true

    input:
    path fastqc_files

    output:
    path "multiqc_report.html",emit:report
    path "multiqc_report_data",emit:data

    script:
    def title='Methyl-Seq Raw FASTQ QC'
    if (task.process.contains('ALIGN')) {
        title='Methyl-Seq Alignment QC'
    } else if (task.process.contains('TRIM')) {
        title='Methyl-Seq Trimmed FASTQ QC'
    }
    """
    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=docker://quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1 \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        multiqc \
            . \
            --force \
            --title "${title}" \
            --filename multiqc_report.html
    """
}
