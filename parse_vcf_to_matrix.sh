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
  -t   threads for bcftools (default: 4)

What it counts:
  - For each sample, counts NON-REFERENCE genotypes per window.
    GT contains any ALT allele: 0/1, 1/1, 0|1, 1|0, 1/2, 2/2, etc.
    Ignores 0/0 and ./.

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

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found in PATH"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH"; exit 1; }

if [[ "$VCF" == *.gz ]]; then
  if [[ ! -f "${VCF}.tbi" && ! -f "${VCF}.csi" ]]; then
    echo "NOTE: $VCF does not appear indexed (.tbi/.csi not found)."
    echo "      Index helps performance: tabix -p vcf $VCF"
  fi
fi

echo "[1/6] Extract contig lengths from VCF header..."
bcftools view -h "$VCF" \
  | awk -F'[=,]' '
      /^##contig=/{
        id=""; len="";
        for(i=1;i<=NF;i++){
          if($i=="ID"){id=$(i+1)}
          if($i=="length"){len=$(i+1)}
        }
        if(id!="" && len!=""){print id"\t"len}
      }' \
  > contigs.genome

echo "[2/6] Create windows (w=${WINDOW_SIZE} bp)..."
bedtools makewindows -g contigs.genome -w "$WINDOW_SIZE" > genome_windows.bed

echo "[3/6] Sample list..."
bcftools query -l "$VCF" > sample_list.txt

echo "[4/6] Count per-sample non-ref variants per window..."
rm -f counts_files.txt
while read -r SAMPLE; do
  echo "  - $SAMPLE"

  bcftools query --threads "$THREADS" -s "$SAMPLE" -f '%CHROM\t%POS\t[%GT]\n' "$VCF" \
    | awk '$3 ~ /(^|[\/|])[1-9]/ {print $1, $2-1, $2}' OFS="\t" \
    > "${SAMPLE}.bed"

  bedtools coverage -a genome_windows.bed -b "${SAMPLE}.bed" -counts \
    | cut -f4 > "${SAMPLE}_counts.tsv"

  echo "${SAMPLE}_counts.tsv" >> counts_files.txt
  rm -f "${SAMPLE}.bed"
done < sample_list.txt

echo "[5/6] Assemble TSV matrix..."
{
  printf "Chrom\tStart\tEnd"
  while read -r SAMPLE; do printf "\t%s" "$SAMPLE"; done < sample_list.txt
  printf "\n"
} > "$OUT"

paste genome_windows.bed $(cat counts_files.txt) >> "$OUT"

echo "[6/6] Cleanup..."
rm -f contigs.genome genome_windows.bed sample_list.txt counts_files.txt *_counts.tsv

echo "Done: $OUT"
