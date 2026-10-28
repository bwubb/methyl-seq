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
    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.samtools_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        bash -c '
            set -euo pipefail
            samtools flagstat -@ ${task.cpus} ${bam} > ${sample}.flagstat.txt
            samtools stats -@ ${task.cpus} ${bam} > ${sample}.stats.txt
            samtools idxstats ${bam} > ${sample}.idxstats.txt
        '

    cat > bsconv_sum.awk << 'END'
BEGIN { OFS="\t"; print "context","retention","conversion","pct_retention" }
\$1 ~ /^[0-9.]+\$/ { for (i=1;i<=8;i++) a[i]+=\$i }
END {
    split("CpA,CpC,CpG,CpT", ctx, ",")
    for (j=0;j<4;j++) {
        r=a[j*2+1]; c=a[j*2+2]; t=r+c
        p=(t>0)?100*r/t:0
        printf "%s\t%.0f\t%.0f\t%.4f\n", ctx[j+1], r+0, c+0, p
    }
}
END

    srun \
        --ntasks=1 \
        --cpus-per-task=${task.cpus} \
        --mem-per-cpu=5632M \
        --container-image=${params.biscuit_container} \
        --container-mounts=${projectDir}:${projectDir} \
        --container-workdir=\$PWD \
        bash -c '
            set -euo pipefail
            biscuit bsstrand ${fasta} ${bam} >/dev/null 2> ${sample}.bsstrand.txt
            biscuit bsconv -p ${fasta} ${bam} | awk -f bsconv_sum.awk > ${sample}.bsconv.txt
        '
    """
}
