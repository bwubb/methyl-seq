process FASTQC {

    tag "${sample}"

    publishDir path:{ "${params.outdir}/${task.process.tokenize(':')[-1].toLowerCase()}" },mode:'copy',overwrite:true

    input:
    tuple val(sample),path(read1),path(read2)

    output:
    path "*_fastqc.html",emit:html
    path "*_fastqc.zip",emit:zip

    script:
    """
    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        fastqc \
            --threads ${task.cpus} \
            --outdir . \
            ${read1} \
            ${read2}
    """
}
