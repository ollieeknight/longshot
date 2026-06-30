#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// ─── Imports ──────────────────────────────────────────────────────────────────
include { PREPROCESS              } from './subworkflows/preprocess'
include { ALIGN                   } from './subworkflows/align'
include { CLASSIFY                } from './subworkflows/classify'
include { QUANTIFY                } from './subworkflows/quantify'
include { EXPORT                  } from './subworkflows/export'
include { COLLECT_INSTRUMENT_STATS } from './modules/qc'
include { MULTIQC                 } from './modules/qc'


// ─── Helpers ──────────────────────────────────────────────────────────────────
def check_file(path, label) {
    if (path == null) error "Required parameter missing: ${label}"
    if (!file(path).exists()) error "File not found — ${label}: ${path}"
}


def preflight_samplesheet(String path) {
    def required = ['experiment', 'library_id', 'run_id', 'bam']
    def lines    = new File(path).readLines()
    if (lines.isEmpty()) error "Samplesheet is empty: ${path}"

    def headers = lines[0].split(',').collect { it.trim() }
    def missing = required - headers
    if (missing) error "Samplesheet missing required columns: ${missing.join(', ')}"

    def rows = lines.tail().withIndex().collect { line, i ->
        def vals = line.split(',').collect { it.trim() }
        if (vals.size() != headers.size())
            error "Samplesheet row ${i + 2}: expected ${headers.size()} fields, got ${vals.size()}"
        [headers, vals].transpose().collectEntries()
    }

    rows.each { row ->
        if (!file(row.bam).exists())
            error "BAM not found for ${row.experiment}/${row.library_id}: ${row.bam}"

        // Validate optional short-read barcodes if provided and not null/NA
        if (row.containsKey('shortread_barcodes') && row.shortread_barcodes) {
            def sr = row.shortread_barcodes.trim()
            def is_null_value = sr == "" || sr.equalsIgnoreCase("null") || sr.equalsIgnoreCase("na") || sr.equalsIgnoreCase("none")
            if (!is_null_value && !file(sr).exists()) {
                error "Short-read barcodes file not found for ${row.experiment}/${row.library_id}: ${sr}"
            }
        }
    }

    // Detect duplicate (experiment, library_id, run_id) combinations
    rows.groupBy { "${it.experiment}__${it.library_id}__${it.run_id}" }.each { key, group ->
        if (group.size() > 1)
            error "Duplicate samplesheet entry: ${key.replace('__', '/')}"
    }
}


// ─── Main Workflow ────────────────────────────────────────────────────────────
workflow {

    // ── 0. Pre-flight ─────────────────────────────────────────────────────────
    if (!params.samplesheet)      error "Please provide --samplesheet <path>"
    if (!params.ref_fasta)        error "Please provide --ref_fasta <path>"
    if (!params.ref_gtf)          error "Please provide --ref_gtf <path>"
    if (!params.adapter_primers)   error "Please provide --adapter_primers <path>"
    if (!params.tenx_3kit_primers) error "Please provide --tenx_3kit_primers <path>"
    if (!params.tenx_whitelist)   error "Please provide --tenx_whitelist <path>"
    if (!params.cage_peaks)       error "Please provide --cage_peaks <path>"
    if (!params.polya_list)       error "Please provide --polya_list <path>"
    if (!params.container_sqanti3) error "Please provide --container_sqanti3 <path to SIF>"

    check_file(params.samplesheet,       "Samplesheet")
    check_file(params.ref_fasta,         "Reference FASTA")
    check_file(params.ref_gtf,           "Reference GTF")
    check_file(params.adapter_primers,    "Adapter primers FASTA")
    check_file(params.tenx_3kit_primers, "10x 3kit primers FASTA")
    check_file(params.tenx_whitelist,    "10x barcode whitelist")
    check_file(params.cage_peaks,        "CAGE peaks BED")
    check_file(params.polya_list,        "PolyA motif list")
    // Validate per-sample adapter kit shorthands resolve to real files
    ['mas8', 'mas12', 'mas16'].each { kit ->
        check_file("${projectDir}/assets/adapters/${kit}_primers.fasta", "Adapter primers FASTA (${kit})")
    }
    preflight_samplesheet(params.samplesheet)

    // ── 1. Parse samplesheet ─────────────────────────────────────────────────
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def bam_file  = file(row.bam)
            def stats_dir = bam_file.parent.parent.resolve("statistics")
            
            def sr_val = row.containsKey('shortread_barcodes') && row.shortread_barcodes ? row.shortread_barcodes.trim() : ""
            def sr_barcodes = (sr_val == "" || sr_val.equalsIgnoreCase("null") || sr_val.equalsIgnoreCase("na") || sr_val.equalsIgnoreCase("none")) ? null : sr_val
            
            def idx_val = row.containsKey('10x_index') && row['10x_index'] ? row['10x_index'].trim() : ""
            def tenx_index = (idx_val == "" || idx_val.equalsIgnoreCase("null") || idx_val.equalsIgnoreCase("na") || idx_val.equalsIgnoreCase("none")) ? null : idx_val

            // Dynamic Chemistry Detection (5prime vs 3prime)
            def chem_val = row.containsKey('chemistry') && row.chemistry ? row.chemistry.trim() : ""
            def chemistry = params.chemistry
            if (chem_val != "" && !chem_val.equalsIgnoreCase("null") && !chem_val.equalsIgnoreCase("na") && !chem_val.equalsIgnoreCase("none")) {
                chemistry = chem_val
            } else if (tenx_index) {
                if (tenx_index.startsWith("SI-TT") || tenx_index.startsWith("SI-TN")) {
                    chemistry = "5prime"
                } else if (tenx_index.startsWith("SI-GA") || tenx_index.startsWith("SI-NA")) {
                    chemistry = "3prime"
                }
            }

            // Adapter kit resolution (mas8 / mas12 / mas16 shorthand or custom path)
            def kit_val = row.containsKey('adapter_kit') && row.adapter_kit ? row.adapter_kit.trim() : ""
            def is_null_kit = kit_val == "" || kit_val.equalsIgnoreCase("null") || kit_val.equalsIgnoreCase("na") || kit_val.equalsIgnoreCase("none")
            def adapter_file = params.adapter_primers
            if (!is_null_kit) {
                if      (kit_val == 'mas8')  adapter_file = "${projectDir}/assets/adapters/mas8_primers.fasta"
                else if (kit_val == 'mas12') adapter_file = "${projectDir}/assets/adapters/mas12_primers.fasta"
                else if (kit_val == 'mas16') adapter_file = "${projectDir}/assets/adapters/mas16_primers.fasta"
                else                         adapter_file = kit_val
            }

            def meta = [
                id:                  "${row.experiment}_${row.library_id}_${row.run_id}",
                sample_id:           "${row.experiment}_${row.library_id}",
                experiment:          row.experiment,
                library_id:          row.library_id,
                run_id:              row.run_id,
                shortread_barcodes:  sr_barcodes,
                tenx_index:          tenx_index,
                chemistry:           chemistry,
                adapter_file:        adapter_file
            ]
            [ meta, bam_file, stats_dir.exists() ? stats_dir : null ]
        }
        .set { ch_parsed }

    ch_parsed
        .map { meta, bam, stats_dir -> [ meta, bam ] }
        .set { ch_raw_bam }

    ch_parsed
        .filter { meta, bam, stats_dir -> stats_dir != null }
        .map    { meta, bam, stats_dir -> [ meta.run_id, stats_dir ] }
        .unique { it[0] }
        .map    { run_id, stats_dir -> [ [id: run_id], stats_dir ] }
        .set { ch_instrument_stats }

    // ── 2. Collect instrument-level QC from raw run directories ──────────────
    COLLECT_INSTRUMENT_STATS(ch_instrument_stats)

    // ── 3. Per-lane preprocessing: HiFi BAM → FLTNC ──────────────────────────
    PREPROCESS(ch_raw_bam)

    // ── 4. Per-sample alignment: merge FLTNC → CRAM ───────────────────────────
    ALIGN(PREPROCESS.out.fltnc_bam)

    // ── 5. Joint Isoform Discovery & SQANTI3 Classification ──────────────────
    CLASSIFY(ALIGN.out.aligned_bam)

    // ── 6. Per-Library Quantification ────────────────────────────────────────
    QUANTIFY(ALIGN.out.aligned_bam, CLASSIFY.out.rescued_gtf)

    // ── 7. Matrix Export & Saturation Curves ─────────────────────────────────
    EXPORT(QUANTIFY.out.counts_dir, CLASSIFY.out.filtered_class, CLASSIFY.out.rescued_gtf)

    // ── 8. Aggregate MultiQC ──────────────────────────────────────────────────
    Channel.empty()
        .mix( COLLECT_INSTRUMENT_STATS.out.stats.flatten() )
        .mix( PREPROCESS.out.lima_reports.flatten() )
        .mix( PREPROCESS.out.flagstat )
        .mix( PREPROCESS.out.skera_logs )
        .mix( PREPROCESS.out.refine_summaries )
        .mix( ALIGN.out.flagstat )
        .mix( ALIGN.out.cramino.flatten() )
        .mix( PREPROCESS.out.cramino_reports.flatten() )
        .mix( EXPORT.out.saturation )
        .mix( PREPROCESS.out.versions.flatten() )
        .mix( ALIGN.out.versions.flatten() )
        .mix( CLASSIFY.out.versions.flatten() )
        .mix( QUANTIFY.out.versions.flatten() )
        .mix( EXPORT.out.versions.flatten() )
        .collect()
        .set { ch_multiqc_reports }

    MULTIQC(ch_multiqc_reports)
}
