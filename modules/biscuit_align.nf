process BISCUIT_ALIGN {

    tag "${sample}"

    publishDir "${params.outdir}/aligned/biscuit",mode:'copy',overwrite:true

    input:
    tuple val(sample),path(read1),path(read2)
    path fasta
    path biscuit_index_files
    val align

    output:
    tuple val(sample),path("${sample}.bam"),path("${sample}.bam.bai"),emit:bam
    path "${sample}.align.log",emit:log

    script:
    def strand_mode=align.directional ? 1 : 0
    def align_cpus=32
    def dup_cpus=8
    def sort_cpus=4
    def rg="@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA"
    """
    set -euo pipefail
    echo "host=\$(hostname) start=\$(date -Is) dupsifter=1" > ${sample}.align.start
    unset SLURM_MEM_PER_NODE
    unset SLURM_MEM_PER_GPU
    mkfifo ${sample}.sam.fifo
    mkfifo ${sample}.dup.fifo

    srun \
        --overlap \
        --ntasks=1 \
        --cpus-per-task=${sort_cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.samtools_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        samtools sort \
            -@ ${sort_cpus} \
            -O BAM \
            -o ${sample}.bam \
            ${sample}.dup.fifo &
    sort_pid=\$!

    srun \
        --overlap \
        --ntasks=1 \
        --cpus-per-task=${dup_cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.dupsifter_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        dupsifter ${fasta} ${sample}.sam.fifo \
        > ${sample}.dup.fifo &
    dup_pid=\$!

    srun \
        --overlap \
        --ntasks=1 \
        --cpus-per-task=${align_cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.biscuit_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        biscuit align \
            -@ ${align_cpus} \
            -b ${strand_mode} \
            -R "${rg}" \
            ${fasta} \
            ${read1} \
            ${read2} \
        > ${sample}.sam.fifo \
        2> ${sample}.align.log

    wait \$dup_pid
    wait \$sort_pid

    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.samtools_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        samtools index \
            -@ ${task.cpus} \
            ${sample}.bam
    """
}
