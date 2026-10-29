process ALIGNMENT_QC {

    tag "${sample}"

    publishDir "${params.outdir}/aligned/qc",mode:'copy',overwrite:true

    input:
    tuple val(sample),path(bam),path(bai)
    path fasta
    path biscuit_index_files

    output:
    path "${sample}.flagstat.txt",emit:flagstat
    path "${sample}.stats.txt",emit:stats
    path "${sample}.idxstats.txt",emit:idxstats
    path "${sample}.bsstrand.txt",emit:bsstrand
    path "${sample}.bsconv.txt",emit:bsconv

    script:
    """
    srun \\
        --ntasks=1 \\
        --cpus-per-task=${task.cpus} \\
        --mem-per-cpu=5632M \\
        --container-image=${params.samtools_container} \\
        --container-mounts=${projectDir}:${projectDir} \\
        --container-workdir=\$PWD \\
        bash -c '
            set -euo pipefail
            samtools flagstat -@ ${task.cpus} ${bam} > ${sample}.flagstat.txt
            samtools stats -@ ${task.cpus} ${bam} > ${sample}.stats.txt
            samtools idxstats ${bam} > ${sample}.idxstats.txt
        '

    srun \\
        --ntasks=1 \\
        --cpus-per-task=${task.cpus} \\
        --mem-per-cpu=5632M \\
        --container-image=${params.biscuit_container} \\
        --container-mounts=${projectDir}:${projectDir} \\
        --container-workdir=\$PWD \\
        bash -c '
            set -euo pipefail
            biscuit bsstrand ${fasta} ${bam} >/dev/null 2> ${sample}.bsstrand.txt
            {
                printf "CpA_R\\tCpA_C\\tCpC_R\\tCpC_C\\tCpG_R\\tCpG_C\\tCpT_R\\tCpT_C\\n"
                biscuit bsconv -p ${fasta} ${bam}
            } > ${sample}.bsconv.txt
        '
    """
}
