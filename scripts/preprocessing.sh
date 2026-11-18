#!/usr/bin/env bash

set -euo pipefail

# Project directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# Samples
SAMPLES=(
    SRR41141592
    SRR41141593
    SRR41141594
    SRR41141595
    SRR41141596
    SRR41141597
)

# Thread settings
THREADS=2
SORT_THREADS=1

# Directories
SRA_DIR="data/raw/sra"
FASTQ_DIR="data/raw/fastq"
TRIMMED_DIR="data/processed/trimmed"
TMP_DIR="data/raw/tmp"
REFERENCE_DIR="data/reference"
INDEX_DIR="$REFERENCE_DIR/hisat2_index"
QC_DIR="results/qc"
RAW_FASTQC_DIR="$QC_DIR/fastqc_raw"
TRIMMED_FASTQC_DIR="$QC_DIR/fastqc_trimmed"
RAW_MULTIQC_DIR="$QC_DIR/multiqc_raw"
TRIMMED_MULTIQC_DIR="$QC_DIR/multiqc_trimmed"
FASTP_DIR="results/fastp"
ALIGN_DIR="results/alignment"
ALIGN_QC_DIR="$ALIGN_DIR/qc"
COUNTS_DIR="results/counts"
LOG_DIR="logs"

# Check existing BAMs
if compgen -G "$ALIGN_DIR/*.sorted.bam" > /dev/null; then
    echo "ERROR: Existing sorted BAM files were detected in $ALIGN_DIR."
    echo "Run this pipeline only in a clean project directory."
    exit 1
fi

mkdir -p \
    "$SRA_DIR" "$FASTQ_DIR" "$TRIMMED_DIR" "$TMP_DIR" \
    "$REFERENCE_DIR" "$INDEX_DIR" "$RAW_FASTQC_DIR" \
    "$TRIMMED_FASTQC_DIR" "$RAW_MULTIQC_DIR" \
    "$TRIMMED_MULTIQC_DIR" "$FASTP_DIR" \
    "$ALIGN_DIR" "$ALIGN_QC_DIR" "$COUNTS_DIR" "$LOG_DIR"

echo "===== 1. DOWNLOAD SRA DATA ====="

for sample in "${SAMPLES[@]}"; do
    if [[ ! -s "$SRA_DIR/$sample/$sample.sra" ]]; then
        prefetch "$sample" -O "$SRA_DIR"
    fi
done

echo "===== 2. CONVERT SRA TO PAIRED FASTQ ====="

for sample in "${SAMPLES[@]}"; do
    mkdir -p "$FASTQ_DIR/$sample"

    if [[ ! -s "$FASTQ_DIR/$sample/${sample}_1.fastq" ||
          ! -s "$FASTQ_DIR/$sample/${sample}_2.fastq" ]]; then
        fasterq-dump "$SRA_DIR/$sample/$sample.sra" \
            --split-files \
            -O "$FASTQ_DIR/$sample" \
            -t "$TMP_DIR" \
            -e 4
    fi
done

echo "===== 3. DOWNLOAD REFERENCE ANNOTATION ====="

GTF_GZ="$REFERENCE_DIR/Homo_sapiens.GRCh38.92.gtf.gz"
GTF="$REFERENCE_DIR/Homo_sapiens.GRCh38.92.gtf"

if [[ ! -s "$GTF_GZ" ]]; then
    curl -fL --retry 3 \
        "https://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz" \
        -o "$GTF_GZ"
fi

if [[ ! -s "$GTF" ]]; then
    gzip -dc "$GTF_GZ" > "$GTF"
fi

echo "===== 4. DOWNLOAD PREBUILT HISAT2 INDEX ====="

INDEX_ARCHIVE="$INDEX_DIR/grch38_tran.tar.gz"
INDEX_PREFIX="$INDEX_DIR/grch38_tran/genome_tran"

INDEX_FILES=(
    "${INDEX_PREFIX}.1.ht2"
    "${INDEX_PREFIX}.2.ht2"
    "${INDEX_PREFIX}.3.ht2"
    "${INDEX_PREFIX}.4.ht2"
    "${INDEX_PREFIX}.5.ht2"
    "${INDEX_PREFIX}.6.ht2"
    "${INDEX_PREFIX}.7.ht2"
    "${INDEX_PREFIX}.8.ht2"
)

index_complete=true

for index_file in "${INDEX_FILES[@]}"; do
    if [[ ! -s "$index_file" ]]; then
        index_complete=false
        break
    fi
done

if [[ "$index_complete" == false ]]; then
    if [[ ! -s "$INDEX_ARCHIVE" ]]; then
        curl -fL --retry 3 \
            "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz" \
            -o "$INDEX_ARCHIVE"
    fi

    gzip -t "$INDEX_ARCHIVE"
    tar -xzf "$INDEX_ARCHIVE" -C "$INDEX_DIR"
fi

for index_file in "${INDEX_FILES[@]}"; do
    if [[ ! -s "$index_file" ]]; then
        echo "ERROR: HISAT2 index file is missing or empty: $index_file"
        exit 1
    fi
done

echo "===== 5. RAW READ QUALITY CONTROL ====="

raw_fastqs=()

for sample in "${SAMPLES[@]}"; do
    raw_fastqs+=(
        "$FASTQ_DIR/$sample/${sample}_1.fastq"
        "$FASTQ_DIR/$sample/${sample}_2.fastq"
    )
done

fastqc "${raw_fastqs[@]}" \
    --outdir "$RAW_FASTQC_DIR" \
    -t "$THREADS"

multiqc "$RAW_FASTQC_DIR" \
    --outdir "$RAW_MULTIQC_DIR"

echo "===== 6. QUALITY FILTERING AND ADAPTER TRIMMING ====="

for sample in "${SAMPLES[@]}"; do
    mkdir -p "$TRIMMED_DIR/$sample"

    fastp \
        -i "$FASTQ_DIR/$sample/${sample}_1.fastq" \
        -I "$FASTQ_DIR/$sample/${sample}_2.fastq" \
        -o "$TRIMMED_DIR/$sample/${sample}_1.fastq.gz" \
        -O "$TRIMMED_DIR/$sample/${sample}_2.fastq.gz" \
        -h "$FASTP_DIR/${sample}_fastp.html" \
        -j "$FASTP_DIR/${sample}_fastp.json" \
        -w "$THREADS"
done

echo "===== 7. POST-TRIMMING QUALITY CONTROL ====="

trimmed_fastqs=()

for sample in "${SAMPLES[@]}"; do
    trimmed_fastqs+=(
        "$TRIMMED_DIR/$sample/${sample}_1.fastq.gz"
        "$TRIMMED_DIR/$sample/${sample}_2.fastq.gz"
    )
done

fastqc "${trimmed_fastqs[@]}" \
    --outdir "$TRIMMED_FASTQC_DIR" \
    -t "$THREADS"

multiqc "$TRIMMED_FASTQC_DIR" \
    --outdir "$TRIMMED_MULTIQC_DIR"

echo "===== 8. CONVERT GTF TO BED12 ====="

python "scripts/gtf_to_bed12.py"

echo "===== 9. ALIGNMENT WITH HISAT2 AND BAM SORTING ====="

for i in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$i]}"

    if (( i < 3 )); then
        group="OSM"
        replicate=$((i + 1))
    else
        group="Control"
        replicate=$((i - 2))
    fi

    sample_id="${group}_${replicate}"

    hisat2 \
        -p "$THREADS" \
        -x "$INDEX_PREFIX" \
        -1 "$TRIMMED_DIR/$sample/${sample}_1.fastq.gz" \
        -2 "$TRIMMED_DIR/$sample/${sample}_2.fastq.gz" \
        2> "$LOG_DIR/${sample_id}_hisat2.log" |
    samtools sort \
        -@ "$SORT_THREADS" \
        -o "$ALIGN_DIR/${sample_id}.sorted.bam" -
done

echo "===== 10. ALIGNMENT QUALITY CONTROL ====="

for bam in "$ALIGN_DIR"/*.sorted.bam; do
    sample="$(basename "$bam" .sorted.bam)"

    samtools quickcheck -v "$bam" \
        > "$ALIGN_QC_DIR/${sample}.quickcheck.txt" 2>&1

    samtools flagstat "$bam" \
        > "$ALIGN_QC_DIR/${sample}.flagstat.txt"

    samtools stats "$bam" \
        > "$ALIGN_QC_DIR/${sample}.stats.txt"
done

echo "===== 11. STRANDEDNESS ASSESSMENT ====="

infer_experiment.py \
    -i "$ALIGN_DIR/OSM_1.sorted.bam" \
    -r "$REFERENCE_DIR/Homo_sapiens.GRCh38.92.bed12" \
    -s 200000 -q 30

infer_experiment.py \
    -i "$ALIGN_DIR/Control_1.sorted.bam" \
    -r "$REFERENCE_DIR/Homo_sapiens.GRCh38.92.bed12" \
    -s 200000 -q 30

echo "===== 12. GENE-LEVEL QUANTIFICATION ====="

featureCounts \
    -T 2 \
    -p \
    --countReadPairs \
    -s 0 \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o "$COUNTS_DIR/gene_counts.txt" \
    "$ALIGN_DIR/OSM_1.sorted.bam" \
    "$ALIGN_DIR/OSM_2.sorted.bam" \
    "$ALIGN_DIR/OSM_3.sorted.bam" \
    "$ALIGN_DIR/Control_1.sorted.bam" \
    "$ALIGN_DIR/Control_2.sorted.bam" \
    "$ALIGN_DIR/Control_3.sorted.bam"

echo "===== PREPROCESSING COMPLETED ====="
echo "Gene count matrix: $COUNTS_DIR/gene_counts.txt"
