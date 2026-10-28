process TRIMGALORE {

    tag "${sample}"

    publishDir "${params.outdir}/trimmed",mode:'copy',overwrite:true

    input:
    tuple val(sample),path(read1),path(read2),val(trim)

    output:
    tuple val(sample),path("${sample}_val_1.fq.gz"),path("${sample}_val_2.fq.gz"),emit:reads
    path "*trimming_report.txt",emit:reports

    script:
    """
    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=docker://quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0 \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        trim_galore \
            --paired \
            --basename ${sample} \
            --adapter ${trim.adapter_r1} \
            --adapter2 ${trim.adapter_r2} \
            --quality ${trim.quality} \
            --length ${trim.min_length} \
            --cores ${task.cpus} \
            ${read1} \
            ${read2}
    """
}
