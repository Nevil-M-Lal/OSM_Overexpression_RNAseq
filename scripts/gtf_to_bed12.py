from collections import defaultdict

gtf_file = "data/reference/Homo_sapiens.GRCh38.92.gtf"
bed_file = "data/reference/Homo_sapiens.GRCh38.92.bed12"

transcripts = defaultdict(list)
exon_count = 0

with open(gtf_file) as gtf:
    for line_number, line in enumerate(gtf, start=1):

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            raise ValueError(f"Invalid GTF record at line {line_number}")

        if fields[2] != "exon":
            continue

        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]

        if start < 1 or end < start:
            raise ValueError(f"Invalid exon coordinates at line {line_number}")

        if strand not in ("+", "-"):
            raise ValueError(f"Invalid strand at line {line_number}")

        attributes = {}

        for item in fields[8].strip().split(";"):
            item = item.strip()

            if not item:
                continue

            parts = item.split(None, 1)

            if len(parts) == 2:
                attributes[parts[0]] = parts[1].strip().strip('"')

        transcript_id = attributes.get("transcript_id")

        if not transcript_id:
            raise ValueError(
                f"Missing transcript_id at GTF line {line_number}"
            )

        transcripts[transcript_id].append(
            (chrom, start, end, strand)
        )

        exon_count += 1


bed_rows = []

for transcript_id, exons in transcripts.items():

    chromosomes = {exon[0] for exon in exons}
    strands = {exon[3] for exon in exons}

    if len(chromosomes) != 1 or len(strands) != 1:
        raise ValueError(
            f"Inconsistent chromosome or strand: {transcript_id}"
        )

    chrom = exons[0][0]
    strand = exons[0][3]

    exons.sort(key=lambda x: x[1])

    chrom_start = min(exon[1] for exon in exons) - 1
    chrom_end = max(exon[2] for exon in exons)

    block_sizes = [
        exon[2] - exon[1] + 1
        for exon in exons
    ]

    block_starts = [
        exon[1] - 1 - chrom_start
        for exon in exons
    ]

    block_count = len(exons)

    bed_row = [
        chrom,
        str(chrom_start),
        str(chrom_end),
        transcript_id,
        "0",
        strand,
        str(chrom_start),
        str(chrom_end),
        "0",
        str(block_count),
        ",".join(map(str, block_sizes)) + ",",
        ",".join(map(str, block_starts)) + ","
    ]

    bed_rows.append(bed_row)


bed_rows.sort(key=lambda row: (row[0], int(row[1]), row[3]))

with open(bed_file, "w") as output:
    for row in bed_rows:
        output.write("\t".join(row) + "\n")

print("GTF-to-BED12 conversion completed.")
print(f"Input exon records: {exon_count:,}")
print(f"Transcript records generated: {len(bed_rows):,}")
print(f"Output: {bed_file}")
