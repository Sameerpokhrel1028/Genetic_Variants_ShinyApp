#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
VCF -> window-based variant count matrix (TSV)

Required:
  -v   input multi-sample VCF (bgzipped .vcf.gz strongly recommended)
  -o   output TSV path

Optional:
  -w   window size in bp (default: 10000 = 10 kb)
  -t   threads for bcftools when supported (default: 4)
       NOTE: bcftools query may not support --threads in some builds.

What it counts:
  - For each sample, counts NON-REFERENCE genotypes per window.
    GT contains any ALT allele: 0/1, 1/1, 0|1, 1|0, 1/2, 2/2, etc.
    Ignores 0/0 and ./. (and missing/blank GT)

Output format:
  Chrom   Start   End   Sample1   Sample2   ...

Examples:
  bash scripts/parse_vcf_to_matrix.sh -v my.vcf.gz -o my_10kb.tsv
  bash scripts/parse_vcf_to_matrix.sh -v my.vcf.gz -w 100000 -t 16 -o my_100kb.tsv
EOF
}

VCF=""
OUT=""
WINDOW_SIZE=10000
THREADS=4

while getopts ":v:o:w:t:h" opt; do
  case $opt in
    v) VCF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    w) WINDOW_SIZE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$VCF" || -z "$OUT" ]]; then
  echo "ERROR: -v and -o are required." >&2
  usage
  exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
  echo "ERROR: bcftools not found in PATH."
  echo "If on HPC, try: module load BCFtools"
  exit 1
fi

if ! command -v bedtools >/dev/null 2>&1; then
  echo "ERROR: bedtools not found in PATH."
  echo "If on HPC, try: module load BEDTools"
  exit 1
fi

echo "Using: $(bcftools --version | head -n1)"
echo "Using: $(bedtools --version)"

if [[ "$VCF" == *.gz ]]; then
  if [[ ! -f "${VCF}.tbi" && ! -f "${VCF}.csi" ]]; then
    echo "NOTE: $VCF does not appear indexed (.tbi/.csi not found)."
    echo "      Index helps performance: tabix -p vcf $VCF"
  fi
fi

WORKDIR="$(mktemp -d -p . parse_vcf_to_matrix.XXXXXX)"
cleanup() { rm -rf "$WORKDIR"; }
trap cleanup EXIT

CONTIGS="$WORKDIR/contigs.genome"
WINDOWS="$WORKDIR/genome_windows.bed"
SAMPLES="$WORKDIR/sample_list.txt"
COUNTS_LIST="$WORKDIR/counts_files.txt"

echo "[1/6] Extract contig lengths from VCF header..."
bcftools view -h "$VCF" \
  | awk -F'[=,<>]' '
      /^##contig=/{
        id=""; len="";
        for(i=1;i<=NF;i++){
          if($i=="ID")     id=$(i+1);
          if($i=="length") len=$(i+1);
        }
        if(id!="" && len ~ /^[0-9]+$/) print id "\t" len;
      }' > "$CONTIGS"

if [[ ! -s "$CONTIGS" ]]; then
  echo "ERROR: contigs.genome is empty. Failed to parse contig lengths from VCF header." >&2
  echo "       Your VCF should contain lines like: ##contig=<ID=chr01,length=12345>" >&2
  echo "       Debug tip: bcftools view -h \"$VCF\" | grep '^##contig' | head" >&2
  exit 1
fi

echo "      Parsed contigs: $(wc -l < "$CONTIGS")"
echo "      Example: $(head -n 1 "$CONTIGS")"

echo "[2/6] Create windows (w=${WINDOW_SIZE} bp)..."
bedtools makewindows -g "$CONTIGS" -w "$WINDOW_SIZE" > "$WINDOWS"

echo "[3/6] Sample list..."
bcftools query -l "$VCF" > "$SAMPLES"

if [[ ! -s "$SAMPLES" ]]; then
  echo "ERROR: No samples found in VCF: $VCF" >&2
  exit 1
fi

echo "[4/6] Count per-sample non-ref variants per window..."
: > "$COUNTS_LIST"

while read -r SAMPLE; do
  echo "  - $SAMPLE"

  SAMPLE_BED="$WORKDIR/${SAMPLE}.bed"
  SAMPLE_COUNTS="$WORKDIR/${SAMPLE}_counts.tsv"

  # NOTE: Some bcftools builds do NOT support --threads for `bcftools query`.
  bcftools query -s "$SAMPLE" -f '%CHROM\t%POS\t[%GT]\n' "$VCF" \
    | awk '($3 ~ /(^|[\/|])[1-9]/) {print $1, $2-1, $2}' OFS="\t" \
    > "$SAMPLE_BED"

  bedtools coverage -a "$WINDOWS" -b "$SAMPLE_BED" -counts \
    | cut -f4 > "$SAMPLE_COUNTS"

  echo "$SAMPLE_COUNTS" >> "$COUNTS_LIST"
done < "$SAMPLES"

echo "[5/6] Assemble TSV matrix..."
{
  printf "Chrom\tStart\tEnd"
  while read -r SAMPLE; do printf "\t%s" "$SAMPLE"; done < "$SAMPLES"
  printf "\n"
} > "$OUT"

paste "$WINDOWS" $(cat "$COUNTS_LIST") >> "$OUT"

echo "[6/6] Done."
echo "Output: $OUT"
