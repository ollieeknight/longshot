process SKERA_SPLIT {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_skera}"

    input:
    tuple val(meta), path(hifi_bam)

    output:
    tuple val(meta), path("${meta.id}_segmented.bam"), emit: segmented_bam
    path "versions.yml",                               emit: versions

    script:
    """
    skera split \\
        -j ${task.cpus} \\
        ${hifi_bam} \\
        ${params.mas16_primers} \\
        ${meta.id}_segmented.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skera: \$(skera --version 2>&1 | grep -oP '(?<=skera )\\S+')
    END_VERSIONS
    """
}


process LIMA_ISOSEQ {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_lima}"

    input:
    tuple val(meta), path(segmented_bam), path(primers)

    output:
    tuple val(meta), path("${meta.id}_fl.bam"), emit: fl_bam
    path "${meta.id}_fl.lima.*",                emit: lima_reports
    path "versions.yml",                        emit: versions

    script:
    """
    lima --isoseq \\
        -j ${task.cpus} \\
        ${segmented_bam} \\
        ${primers} \\
        fl.bam

    # Rename oriented BAM to a predictable name regardless of primer labels in FASTA
    # Exactly one non-unassigned oriented BAM must exist for --isoseq mode
    NBAM=\$(find . -maxdepth 1 -name "fl.*.bam" ! -name "fl.bam" ! -name "*unassigned*" | wc -l)
    if [ "\$NBAM" -ne 1 ]; then
        echo "ERROR [LIMA_ISOSEQ]: Expected exactly 1 oriented BAM output from lima, found \$NBAM" >&2
        find . -maxdepth 1 -name "fl.*.bam" >&2
        exit 1
    fi
    find . -maxdepth 1 -name "fl.*.bam" ! -name "fl.bam" ! -name "*unassigned*" \
        | xargs -I{} mv {} ${meta.id}_fl.bam
    samtools index ${meta.id}_fl.bam

    # Rename lima reports
    for f in fl.lima.*; do
        mv "\$f" "${meta.id}_\$f"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lima: \$(lima --version 2>&1 | grep -oP '(?<=lima )\\S+')
    END_VERSIONS
    """
}


process LIMA_MULTIPLEX {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_lima}"

    input:
    tuple val(meta), path(segmented_bam), path(primers)

    output:
    tuple val(meta), path("fl.*.bam"),   emit: split_bams
    path "${meta.id}_fl.lima.*",        emit: lima_reports
    path "versions.yml",                emit: versions

    script:
    """
    lima --isoseq \\
        -j ${task.cpus} \\
        ${segmented_bam} \\
        ${primers} \\
        fl.bam

    # Rename lima reports
    for f in fl.lima.*; do
        mv "\$f" "${meta.id}_\$f"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lima: \$(lima --version 2>&1 | grep -oP '(?<=lima )\\S+')
    END_VERSIONS
    """
}


process ISOSEQ_TAG {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_isoseq}"

    input:
    tuple val(meta), path(fl_bam)

    output:
    tuple val(meta), path("${meta.id}_flt.bam"), emit: flt_bam
    path "versions.yml",                         emit: versions

    script:
    """
    isoseq tag \\
        --design T-12U-16B \\
        -j ${task.cpus} \\
        ${fl_bam} \\
        ${meta.id}_flt.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$(isoseq --version 2>&1 | grep -oP '(?<=isoseq )\\S+')
    END_VERSIONS
    """
}


process ISOSEQ_REFINE {
    tag "${meta.id}"
    label 'process_high'
    container "${params.container_isoseq}"

    input:
    tuple val(meta), path(flt_bam)

    output:
    tuple val(meta), path("${meta.id}_fltnc.bam"), emit: fltnc_bam
    path "versions.yml",                           emit: versions

    script:
    """
    isoseq refine \\
        --require-polya \\
        -j ${task.cpus} \\
        ${flt_bam} \\
        ${params.tenx_3kit_primers} \\
        ${meta.id}_fltnc.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$(isoseq --version 2>&1 | grep -oP '(?<=isoseq )\\S+')
    END_VERSIONS
    """
}


process DETECT_SAMPLE_INDICES {
    tag "${meta.id}"
    label 'process_low'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(raw_bam), val(requested_indices)

    output:
    tuple val(meta), path("${meta.id}_index_report.tsv"), emit: index_report

    script:
    """
    samtools view ${raw_bam} | head -n 20000 | awk '{print ">read_"NR"\\n"\$10}' > temp.fasta

    python3 -c "
import os
import csv
import sys

def rev_comp(seq):
    tb = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(tb)[::-1]

index_db = {}
index_dir = '${projectDir}/assets/indexes'
for f in os.listdir(index_dir):
    if f.endswith('.csv') and 'truseq' not in f.lower():
        with open(os.path.join(index_dir, f), 'r') as fh:
            reader = csv.reader(fh)
            for row in reader:
                if len(row) >= 2:
                    index_id = row[0]
                    sequences = [s.strip().upper() for s in row[1:] if s.strip()]
                    index_db[index_id] = sequences

reads = []
if os.path.exists('temp.fasta'):
    with open('temp.fasta', 'r') as fh:
        current_seq = ''
        for line in fh:
            if line.startswith('>'):
                if current_seq:
                    reads.append(current_seq)
                current_seq = ''
            else:
                current_seq += line.strip().upper()
        if current_seq:
            reads.append(current_seq)

matches = {idx: 0 for idx in index_db}
for read in reads:
    read_sub = read[:100]
    for idx_id, seqs in index_db.items():
        found = False
        for seq in seqs:
            if seq in read_sub or rev_comp(seq) in read_sub:
                found = True
                break
        if found:
            matches[idx_id] += 1
            break

total_reads = len(reads) if reads else 1
index_report = []
for idx_id, cnt in matches.items():
    pct = (cnt / total_reads) * 100
    if pct >= 0.5:
        index_report.append(f'{idx_id}\\t{cnt}\\t{pct:.2f}%')

with open('${meta.id}_index_report.tsv', 'w') as fh:
    fh.write('index_id\\treads\\tpercentage\\n')
    for r in index_report:
        fh.write(r + '\\n')

requested = [r.strip() for r in '${requested_indices}'.split(',') if r.strip()]
if requested:
    print(f'Validating requested indices: {requested}')
    for req in requested:
        cnt = matches.get(req, 0)
        pct = (cnt / total_reads) * 100
        print(f'Requested index {req}: {cnt} reads ({pct:.2f}%)')
        if cnt == 0:
            sys.stderr.write(f'ERROR: Specified 10x index \"{req}\" was not detected in the raw reads. Please check your samplesheet index settings!\\n')
            sys.exit(1)
"
    """
}


process CONSTRUCT_MULTIPLEX_PRIMERS {
    tag "${meta.id}"
    label 'process_low'
    container "${params.container_samtools}"

    input:
    tuple val(meta), val(index_mappings) // index_mappings is a list of [library_id, 10x_index]

    output:
    tuple val(meta), path("${meta.id}_multiplex_primers.fasta"), emit: primers

    script:
    """
    python3 -c "
import os
import csv

# 1. Parse standard primers file to get 5p/3p core sequences
# Auto-detect whether it is 5' kit or 3' kit
primers_file = '${params.tenx_5kit_primers}'
is_5prime = '5kit' in primers_file.lower()

p_5p = ''
p_3p = ''
with open(primers_file, 'r') as f:
    current_label = ''
    for line in f:
        if line.startswith('>'):
            current_label = line.strip()[1:]
        else:
            if '5p' in current_label:
                p_5p += line.strip()
            elif '3p' in current_label:
                p_3p += line.strip()

# 2. Parse all index databases
index_db = {}
index_dir = '${projectDir}/assets/indexes'
for f in os.listdir(index_dir):
    if f.endswith('.csv') and 'truseq' not in f.lower():
        with open(os.path.join(index_dir, f), 'r') as fh:
            reader = csv.reader(fh)
            for row in reader:
                if len(row) >= 2:
                    index_id = row[0]
                    sequences = [s.strip().upper() for s in row[1:] if s.strip()]
                    index_db[index_id] = sequences

# 3. Construct custom multiplex FASTA
# We format as:
# >[library_id]_[1-4]_5p
# [primer_sequence]
# >[library_id]_3p
# [primer_sequence]
output_file = '${meta.id}_multiplex_primers.fasta'
mappings = [m.split(':') for m in '${index_mappings}'.replace('[','').replace(']','').split(',') if m.strip()]

with open(output_file, 'w') as fh:
    for lib_id, idx_id in mappings:
        lib_id = lib_id.strip()
        idx_id = idx_id.strip()
        seqs = index_db.get(idx_id, [])
        
        if is_5prime:
            # For 5prime GEX, index is appended to 5p primer
            for idx, seq in enumerate(seqs):
                fh.write(f'>{lib_id}_{idx+1}_5p\\n{p_5p}{seq}\\n')
            fh.write(f'>{lib_id}_3p\\n{p_3p}\\n')
        else:
            # For 3prime GEX, index is appended to 3p primer
            fh.write(f'>{lib_id}_5p\\n{p_5p}\\n')
            for idx, seq in enumerate(seqs):
                fh.write(f'>{lib_id}_{idx+1}_3p\\n{p_3p}{seq}\\n')
"
    """
}


process MERGE_INDEX_BAMS {
    tag "${meta.library_id}"
    label 'process_medium'
    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.library_id}_fl.bam"), emit: fl_bam

    script:
    """
    samtools merge -f -@ ${task.cpus} ${meta.library_id}_fl.bam ${bams}
    samtools index -@ ${task.cpus} ${meta.library_id}_fl.bam
    """
}
