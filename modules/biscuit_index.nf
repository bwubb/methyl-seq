process BISCUIT_INDEX {

    tag "${fasta.simpleName}"

    publishDir "${params.outdir}/reference/biscuit",mode:'copy',overwrite:false

    input:
    path fasta

    output:
    path fasta,emit:fasta
    path "${fasta}.*",emit:index_files

    script:
    """
    set -euo pipefail
    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.samtools_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        samtools faidx ${fasta}

    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.biscuit_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        biscuit index ${fasta}
    """
}
