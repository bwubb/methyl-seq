process METHYLDACKEL {

    tag "${sample}"

    publishDir "${params.outdir}/methylation/methyldackel",mode:'copy',overwrite:true

    input:
    tuple val(sample),path(bam),path(bai)
    path fasta
    path fai

    output:
    path "${sample}_CpG.bedGraph",emit:cpg
    path "${sample}.spikein_methylation.tsv",emit:spikein
    path "${sample}.methyl_bounds.txt",emit:bounds
    path "${sample}.mbias.log",emit:mbias_log
    path "${sample}.mbias*.svg",emit:mbias_svg,optional:true

    script:
    def clip=params.methyl_end_clip
    def bounds="${clip},${clip},${clip},${clip}"
    def contigs=params.spikein_contigs.toString().split(',').collect { it.trim() }.join(' ')
    """
    set -euo pipefail
    printf 'end_clip_bases=%s\\nnOT=%s\\nnOB=%s\\n' '${clip}' '${bounds}' '${bounds}' > ${sample}.methyl_bounds.txt

    srun \\
        --ntasks=1 \\
        --cpus-per-task=${task.cpus} \\
        --mem-per-cpu=5632M \\
        --container-image=${params.methyldackel_container} \\
        --container-mounts=${projectDir}:${projectDir} \\
        --container-workdir=\$PWD \\
        bash -c '
            set -euo pipefail
            MethylDackel mbias \\
                ${fasta} \\
                ${bam} \\
                ${sample}.mbias \\
                > ${sample}.mbias.log 2>&1
            MethylDackel extract \\
                -@ ${task.cpus} \\
                -o ${sample} \\
                --mergeContext \\
                --minDepth ${params.methyldackel_min_depth} \\
                --nOT ${bounds} \\
                --nOB ${bounds} \\
                ${fasta} \\
                ${bam}
        '

    printf 'contig\\tmethylated_calls\\tunmethylated_calls\\ttotal_calls\\tpct_methylated\\n' > ${sample}.spikein_methylation.tsv
    for contig in ${contigs}; do
        awk -v c="\$contig" '\$1==c && NF>=6 { m+=\$5; u+=\$6 } END { t=m+u; if (t>0) p=100*m/t; else p=0; printf "%s\\t%d\\t%d\\t%d\\t%.4f\\n", c, m, u, t, p }' \\
            ${sample}_CpG.bedGraph >> ${sample}.spikein_methylation.tsv
    done
    """
}
