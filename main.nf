nextflow.enable.dsl=2

include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { MULTIQC as MULTIQC_RAW } from './modules/multiqc'
include { MULTIQC as MULTIQC_TRIMMED } from './modules/multiqc'
include { MULTIQC as MULTIQC_ALIGNMENT } from './modules/multiqc'
include { TRIMGALORE } from './modules/trimgalore'
include { BISCUIT_INDEX } from './modules/biscuit_index'
include { BISCUIT_ALIGN } from './modules/biscuit_align'
include { ALIGNMENT_QC } from './modules/alignment_qc'
include { METHYLDACKEL } from './modules/methyldackel'


def isFastq(path) {
    def fastq_extensions=['.fastq.gz','.fq.gz']
    def name=path.toString().toLowerCase()
    return fastq_extensions.any { name.endsWith(it) }
}


def classifyFastqHeader(path) {
    def stream=path.newInputStream()
    int b1=stream.read()
    int b2=stream.read()
    stream.close()

    if (b1<0) {
        return 'empty'
    }
    if (b1==0x1f && b2==0x8b) {
        return 'ok'
    }
    if (b1==0) {
        return 'leading_zeros'
    }
    def h1=String.format('%02x',b1)
    def h2=(b2<0) ? 'eof' : String.format('%02x',b2)
    return "not_gzip:${h1},${h2}"
}


def loadLibraryPreps(path) {
    def yaml=new org.yaml.snakeyaml.Yaml()
    return yaml.load(file(path,checkIfExists:true).text)
}


def resolveTrimParams() {
    def preps=loadLibraryPreps("${projectDir}/${params.library_preps}")
    if (!preps.containsKey(params.library_prep)) {
        error "Unknown library_prep '${params.library_prep}'. Known: ${preps.keySet().join(',')}"
    }

    def prep=preps[params.library_prep]
    def adapter_r1=params.adapter_r1 ?: prep.adapters.r1
    def adapter_r2=params.adapter_r2 ?: prep.adapters.r2
    def quality=params.trim_quality ?: prep.trimming.quality
    def min_length=params.trim_length ?: prep.trimming.min_length

    if (!adapter_r1 || !adapter_r2) {
        error "No adapter sequences for library_prep '${params.library_prep}'. Set --adapter_r1 and --adapter_r2 or add them to ${params.library_preps}."
    }

    return [
        adapter_r1:adapter_r1,
        adapter_r2:adapter_r2,
        quality:quality,
        min_length:min_length
    ]
}


def resolveAlignParams() {
    def preps=loadLibraryPreps("${projectDir}/${params.library_preps}")
    if (!preps.containsKey(params.library_prep)) {
        error "Unknown library_prep '${params.library_prep}'. Known: ${preps.keySet().join(',')}"
    }

    def prep=preps[params.library_prep]
    def aligner=params.aligner ?: prep.alignment.aligner
    def directional=(params.directional!=null) ? params.directional : prep.alignment.directional

    if (aligner!='biscuit') {
        error "aligner '${aligner}' is not implemented. v1 is BISCUIT_ALIGN; Bismark/bwa-meth come later."
    }

    return [
        aligner:aligner,
        directional:directional
    ]
}


def findFai(indexFiles) {
    def files=(indexFiles instanceof Iterable) ? indexFiles : [indexFiles]
    def fai=files.find { it.name.endsWith('.fai') }
    if (!fai) {
        error "No .fai next to the FASTA. Run samtools faidx on it before methylation calling."
    }
    return fai
}


def existingBiscuitIndex(fasta) {
    def fa=file(fasta,checkIfExists:true)
    def idx=fa.parent.listFiles().findAll { f ->
        f.name.startsWith(fa.name + '.')
    }
    def ready=idx.any { it.name.endsWith('.dau.bwt') || it.name.endsWith('.par.bwt') || it.name.endsWith('.bis.bwt') }
    return [
        fasta:fa,
        index_files:idx,
        ready:ready
    ]
}


def writeSamplesheetQc(valid,invalid) {
    def report_dir=file("${params.outdir}/samplesheet_qc")
    report_dir.mkdirs()

    def valid_lines=['sample\tfastq_1\tfastq_2']
    valid.each { row ->
        valid_lines << "${row[0]}\t${row[3]}\t${row[4]}"
    }
    file("${report_dir}/valid_samples.tsv").text=valid_lines.join('\n') + '\n'

    def invalid_lines=['sample\tfastq_1\tfastq_2\tfastq_1_status\tfastq_2_status']
    invalid.each { row ->
        invalid_lines << "${row[0]}\t${row[1]}\t${row[2]}\t${row[3]}\t${row[4]}"
    }
    file("${report_dir}/invalid_samples.tsv").text=invalid_lines.join('\n') + '\n'
}


def validateRow(row) {
    def expected_columns=['sample','fastq_1','fastq_2'] as Set
    def columns=row.keySet().findAll { it } as Set
    def extra=columns - expected_columns
    def missing=expected_columns - columns

    if (missing) {
        error "Samplesheet missing required columns: ${missing.join(',')}"
    }

    if (extra) {
        error "Samplesheet has unexpected columns: ${extra.join(',')}"
    }

    if (!row.sample) {
        error "Missing sample ID in samplesheet"
    }

    def sample=row.sample.toString().trim()
    def fastq_1=row.fastq_1?.toString()?.trim()
    def fastq_2=row.fastq_2?.toString()?.trim()

    if (!fastq_1 || !fastq_2) {
        error "Missing FASTQ path for sample ${sample}"
    }

    if (fastq_1==fastq_2) {
        error "R1 and R2 paths are identical for sample ${sample}: ${fastq_1}"
    }

    if (!isFastq(fastq_1)) {
        error "Malformed R1 extension for sample ${sample}: ${fastq_1}"
    }

    if (!isFastq(fastq_2)) {
        error "Malformed R2 extension for sample ${sample}: ${fastq_2}"
    }

    if (!fastq_1.contains('_R1_')) {
        error "FASTQ labeled fastq_1 does not look like R1 for sample ${sample}: ${fastq_1}"
    }

    if (!fastq_2.contains('_R2_')) {
        error "FASTQ labeled fastq_2 does not look like R2 for sample ${sample}: ${fastq_2}"
    }

    tuple(
        sample,
        file(fastq_1,checkIfExists:true),
        file(fastq_2,checkIfExists:true),
        fastq_1,
        fastq_2
    )
}


workflow {

    def trim=resolveTrimParams()
    def align=resolveAlignParams()
    log.info "library_prep=${params.library_prep} adapter_r1=${trim.adapter_r1} adapter_r2=${trim.adapter_r2} trim_quality=${trim.quality} trim_length=${trim.min_length}"
    log.info "aligner=${align.aligner} directional=${align.directional} fasta=${params.fasta}"

    /*
     * Read and validate samplesheet.
     * Gzip-header failures are reported and dropped; remaining
     * samples continue to raw QC, Trim Galore, and trimmed QC.
     */
    samples_ch=Channel
        .fromPath(params.input,checkIfExists:true)
        .splitCsv(
            header:true,
            sep:'\t'
        )
        .map { row ->
            validateRow(row)
        }
        .toList()
        .flatMap { rows ->
            if (rows.isEmpty()) {
                error "No samples found in ${params.input}"
            }

            def duplicates=rows
                .collect { it[0] }
                .countBy { it }
                .findAll { sample,n -> n > 1 }
                .keySet()

            if (duplicates) {
                error "Duplicate sample IDs in samplesheet: ${duplicates.join(',')}"
            }

            def valid=[]
            def invalid=[]
            rows.each { sample,r1,r2,p1,p2 ->
                def s1=classifyFastqHeader(r1)
                def s2=classifyFastqHeader(r2)
                if (s1=='ok' && s2=='ok') {
                    valid << [sample,r1,r2,p1,p2]
                } else {
                    invalid << [sample,p1,p2,s1,s2]
                }
            }

            writeSamplesheetQc(valid,invalid)

            log.info "FASTQ gzip header check: ${valid.size()} valid, ${invalid.size()} invalid"
            invalid.each { row ->
                log.warn "Invalid FASTQ sample ${row[0]}: R1=${row[3]} R2=${row[4]}"
            }

            if (valid.isEmpty()) {
                error "No valid gzip FASTQs. See ${params.outdir}/samplesheet_qc/invalid_samples.tsv"
            }

            valid.collect { row ->
                tuple(row[0],row[1],row[2])
            }
        }


    /*
     * Raw FASTQ quality control
     */
    FASTQC_RAW(samples_ch)

    MULTIQC_RAW(
        FASTQC_RAW.out.zip
            .flatten()
            .collect()
    )


    /*
     * Adapter / quality trim. Kit adapters come from
     * conf/library_preps.yml via params.library_prep.
     * No fixed end clipping unless later QC / M-bias says so.
     */
    TRIMGALORE(samples_ch.combine(Channel.value(trim)))

    FASTQC_TRIMMED(TRIMGALORE.out.reads)

    MULTIQC_TRIMMED(
        FASTQC_TRIMMED.out.zip
            .flatten()
            .collect()
    )


    /*
     * BISCUIT index once, then align trimmed reads.
     * Reuse an index already sitting next to the FASTA
     * (e.g. the Betty srun you started by hand).
     */
    def biscuit_idx=existingBiscuitIndex(params.fasta)
    if (biscuit_idx.ready) {
        log.info "Using existing BISCUIT index next to ${biscuit_idx.fasta}"
        ch_fasta=Channel.value(biscuit_idx.fasta)
        ch_index_files=Channel.value(biscuit_idx.index_files)
        ch_fai=Channel.value(findFai(biscuit_idx.index_files))
    } else {
        log.info "No BISCUIT index found; running BISCUIT_INDEX"
        BISCUIT_INDEX(Channel.fromPath(params.fasta,checkIfExists:true))
        ch_fasta=BISCUIT_INDEX.out.fasta
        ch_index_files=BISCUIT_INDEX.out.index_files
        ch_fai=BISCUIT_INDEX.out.index_files.map { idx -> findFai(idx) }
    }

    BISCUIT_ALIGN(
        TRIMGALORE.out.reads,
        ch_fasta,
        ch_index_files,
        Channel.value(align)
    )

    ALIGNMENT_QC(
        BISCUIT_ALIGN.out.bam,
        ch_fasta,
        ch_index_files
    )

    METHYLDACKEL(
        BISCUIT_ALIGN.out.bam,
        ch_fasta,
        ch_fai
    )

    MULTIQC_ALIGNMENT(
        ALIGNMENT_QC.out.flagstat
            .mix(ALIGNMENT_QC.out.stats)
            .mix(ALIGNMENT_QC.out.idxstats)
            .mix(ALIGNMENT_QC.out.bsstrand)
            .mix(ALIGNMENT_QC.out.bsconv)
            .collect()
    )
}
