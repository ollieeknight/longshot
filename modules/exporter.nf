process EXPORT_LIBRARY_MTX {
    tag "${meta.sample_id}"
    label 'process_medium'
    container "${params.container_multiqc}" // Python 3 and pandas are installed in multiqc container
    publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc_export" }, mode: 'copy'

    input:
    tuple val(meta), path(counts_dir), path(filtered_class)

    output:
    tuple val(meta), path("gene/"), path("transcript/"), emit: mtx_export
    path "versions.yml",                                  emit: versions

    script:
    """
    python3 -c "
import os
import gzip
import shutil
import pandas as pd

counts_dir = '${counts_dir}'
filtered_class = '${filtered_class}'
ref_gtf = '${params.ref_gtf}'
library_id = '${meta.library_id}'

# Create export folders
os.makedirs('gene', exist_ok=True)
os.makedirs('transcript', exist_ok=True)

# 1. Parse Ensembl ID to Gene Symbol mapping from GTF
gene_map = {}
if os.path.exists(ref_gtf):
    print('Parsing GTF for Gene Symbol mapping...')
    open_fun = gzip.open if ref_gtf.endswith('.gz') else open
    with open_fun(ref_gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attrs = fields[8]
            gene_id = ''
            gene_name = ''
            # Extract gene_id
            if 'gene_id \"' in attrs:
                gene_id = attrs.split('gene_id \"')[1].split('\"')[0]
            # Extract gene_name
            if 'gene_name \"' in attrs:
                gene_name = attrs.split('gene_name \"')[1].split('\"')[0]
            if gene_id and gene_name:
                gene_map[gene_id] = gene_name
                gene_map[gene_id.split('.')[0]] = gene_name

# 2. Parse SQANTI3 filtering classification
sqanti_meta = {}
if os.path.exists(filtered_class):
    print('Parsing SQANTI3 classifications...')
    df_sq = pd.read_csv(filtered_class, sep='\\t')
    df_sq['isoform_seurat'] = df_sq['isoform'].str.replace('_', '-')
    for _, row in df_sq.iterrows():
        sqanti_meta[row['isoform']] = {
            'associated_gene': str(row.get('associated_gene', '')),
            'associated_transcript': str(row.get('associated_transcript', '')),
            'structural_category': str(row.get('structural_category', ''))
        }
        sqanti_meta[row['isoform_seurat']] = sqanti_meta[row['isoform']]

# 3. Export Gene Matrix
print('Exporting Gene Matrix...')
raw_gene_mtx = os.path.join(counts_dir, f'{library_id}.gene_grouped_barcode_counts.matrix.mtx')
raw_gene_bc = os.path.join(counts_dir, f'{library_id}.gene_grouped_barcode_counts.barcodes.tsv')
raw_gene_feat = os.path.join(counts_dir, f'{library_id}.gene_grouped_barcode_counts.features.tsv')

if os.path.exists(raw_gene_mtx):
    with open(raw_gene_mtx, 'rb') as f_in, gzip.open('gene/matrix.mtx.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    with open(raw_gene_bc, 'rb') as f_in, gzip.open('gene/barcodes.tsv.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    with open(raw_gene_feat, 'r') as f_in, gzip.open('gene/features.tsv.gz', 'wt') as f_out:
        for line in f_in:
            gene_id = line.strip()
            gene_symbol = gene_map.get(gene_id, gene_map.get(gene_id.split('.')[0], gene_id))
            f_out.write(f'{gene_id}\\t{gene_symbol}\\tGene Expression\\n')

# 4. Export Transcript Matrix
print('Exporting Transcript Matrix...')
raw_iso_mtx = os.path.join(counts_dir, f'{library_id}.transcript_grouped_barcode_counts.matrix.mtx')
raw_iso_bc = os.path.join(counts_dir, f'{library_id}.transcript_grouped_barcode_counts.barcodes.tsv')
raw_iso_feat = os.path.join(counts_dir, f'{library_id}.transcript_grouped_barcode_counts.features.tsv')

if os.path.exists(raw_iso_mtx):
    with open(raw_iso_mtx, 'rb') as f_in, gzip.open('transcript/matrix.mtx.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    with open(raw_iso_bc, 'rb') as f_in, gzip.open('transcript/barcodes.tsv.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    with open(raw_iso_feat, 'r') as f_in, gzip.open('transcript/features.tsv.gz', 'wt') as f_out:
        for line in f_in:
            tx_id = line.strip()
            tx_seurat = tx_id.replace('_', '-')
            meta = sqanti_meta.get(tx_id, sqanti_meta.get(tx_seurat, {}))
            assoc_gene = meta.get('associated_gene', '')
            assoc_tx = meta.get('associated_transcript', '')
            struct_cat = meta.get('structural_category', '')
            gene_symbol = gene_map.get(assoc_gene, gene_map.get(assoc_gene.split('.')[0], assoc_gene))
            if not gene_symbol:
                gene_symbol = assoc_gene
            f_out.write(f'{tx_seurat}\\t{gene_symbol}\\tGene Expression\\t{assoc_gene}\\t{assoc_tx}\\t{struct_cat}\\n')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oP '(?<=Python )\\S+')
    END_VERSIONS
    """
}


process GENERATE_SHARED_CATALOG {
    tag "${experiment}"
    label 'process_medium'
    container "${params.container_multiqc}"
    publishDir { "${params.outdir}/${experiment}/joint" }, mode: 'copy'

    input:
    tuple val(experiment), path(filtered_gtf), path(filtered_class)

    output:
    path "shared_isoform_catalog.tsv.gz", emit: catalog
    path "shared_isoform_map.tsv.gz",     emit: map
    path "versions.yml",                  emit: versions

    script:
    """
    python3 -c "
import os
import gzip
import pandas as pd

# Load classification and replace potential problematic NAs
df = pd.read_csv('${filtered_class}', sep='\\t')
df = df.replace('NA', '').fillna('')

# Predict cell-barcode/isoform name dashes conversion
df['isoform'] = df['isoform'].str.replace('_', '-')

# Export cleaned catalog
df.to_csv('shared_isoform_catalog.tsv.gz', sep='\\t', index=False, compression='gzip')

# Build mapping catalog
gene_map = {}
ref_gtf = '${params.ref_gtf}'
if os.path.exists(ref_gtf):
    open_fun = gzip.open if ref_gtf.endswith('.gz') else open
    with open_fun(ref_gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attrs = fields[8]
            gene_id = ''
            gene_name = ''
            if 'gene_id \"' in attrs:
                gene_id = attrs.split('gene_id \"')[1].split('\"')[0]
            if 'gene_name \"' in attrs:
                gene_name = attrs.split('gene_name \"')[1].split('\"')[0]
            if gene_id and gene_name:
                gene_map[gene_id] = gene_name
                gene_map[gene_id.split('.')[0]] = gene_name

# Write shared mapping map
df_map = df[['isoform', 'associated_gene', 'associated_transcript', 'structural_category']].copy()
df_map['associated_gene_symbol'] = df_map['associated_gene'].apply(lambda x: gene_map.get(x, gene_map.get(x.split('.')[0], x)))
df_map.to_csv('shared_isoform_map.tsv.gz', sep='\\t', index=False, compression='gzip')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oP '(?<=Python )\\S+')
    END_VERSIONS
    """
}


process CALCULATE_SATURATION {
    tag "${meta.sample_id}"
    label 'process_medium'
    container "${params.container_multiqc}"
    publishDir { "${params.outdir}/${meta.experiment}/${meta.library_id}/qc" }, mode: 'copy'

    input:
    tuple val(meta), path(counts_dir)

    output:
    path "${meta.sample_id}_saturation.tsv", emit: report
    path "versions.yml",                     emit: versions

    script:
    """
    python3 -c "
import os
import random
import pandas as pd

# Search for read_to_transcript or read_assignments in counts_dir
tsv_file = None
for f in os.listdir('${counts_dir}'):
    if 'read_to_transcript' in f or 'read_assignments' in f or f.endswith('.read_to_transcript.tsv'):
        tsv_file = os.path.join('${counts_dir}', f)
        break

if not tsv_file or not os.path.exists(tsv_file):
    print('No read assignments TSV found, writing dummy saturation table.')
    with open('${meta.sample_id}_saturation.tsv', 'w') as f:
        f.write('fraction\\ttotal_reads\\tgenes\\tisoforms\\tcells\\tumis\\n')
        f.write('1.0\\t0\\t0\\t0\\t0\\t0\\n')
else:
    print(f'Reading assignments from {tsv_file}...')
    # Load columns: read_id (or read), barcode (or CB), UMI (or XM), transcript_id, gene_id
    df = pd.read_csv(tsv_file, sep='\\t')
    
    # Auto-detect column names
    col_bc = next((c for c in df.columns if c in ['barcode', 'CB', 'cell_barcode']), None)
    col_umi = next((c for c in df.columns if c in ['umi', 'XM', 'UMI']), None)
    col_tx = next((c for c in df.columns if c in ['transcript_id', 'transcript', 'tx_id']), None)
    col_gene = next((c for c in df.columns if c in ['gene_id', 'gene']), None)
    
    if not col_bc or not col_tx:
        print('Required columns not found, writing dummy saturation.')
        with open('${meta.sample_id}_saturation.tsv', 'w') as f:
            f.write('fraction\\ttotal_reads\\tgenes\\tisoforms\\tcells\\tumis\\n')
    else:
        results = []
        n_reads = len(df)
        fractions = [0.1, 0.25, 0.5, 0.75, 1.0]
        
        for frac in fractions:
            if frac == 1.0:
                df_sub = df
            else:
                df_sub = df.sample(frac=frac, random_state=42)
            
            total_reads = len(df_sub)
            n_genes = df_sub[col_gene].nunique() if col_gene else 0
            n_isoforms = df_sub[col_tx].nunique()
            n_cells = df_sub[col_bc].nunique()
            n_umis = df_sub[col_umi].nunique() if col_umi else 0
            
            results.append({
                'fraction': frac,
                'total_reads': total_reads,
                'genes': n_genes,
                'isoforms': n_isoforms,
                'cells': n_cells,
                'umis': n_umis
            })
            
        df_res = pd.DataFrame(results)
        df_res.to_csv('${meta.sample_id}_saturation.tsv', sep='\\t', index=False)
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oP '(?<=Python )\\S+')
    END_VERSIONS
    """
}
