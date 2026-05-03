#!/usr/bin/env bash
set -euo pipefail

VCF=""
OUT=""
WINDOW_SIZE=10000

while getopts ":v:o:w:h" opt; do
  case $opt in
    v) VCF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    w) WINDOW_SIZE="$OPTARG" ;;
    h) echo "Usage: bash parse_vcf_to_matrix.sh -v input.vcf.gz -o output.tsv [-w 10000]"; exit 0 ;;
    *) echo "Usage: bash parse_vcf_to_matrix.sh -v input.vcf.gz -o output.tsv [-w 10000]"; exit 1 ;;
  esac
done

[[ -z "$VCF" || -z "$OUT" ]] && { echo "ERROR: -v and -o required"; exit 1; }

WORKDIR="$(mktemp -d -p . parse_vcf_to_matrix.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

CONTIGS="$WORKDIR/contigs.genome"
WINDOWS="$WORKDIR/windows.bed"
SAMPLES="$WORKDIR/samples.txt"
COUNTS_LIST="$WORKDIR/counts.list"

bcftools view -h "$VCF" | awk -F'[=,<>]' '
/^##contig=/{
  id=""; len="";
  for(i=1;i<=NF;i++){
    if($i=="ID") id=$(i+1);
    if($i=="length") len=$(i+1);
  }
  if(id!="" && len ~ /^[0-9]+$/) print id "\t" len;
}' > "$CONTIGS"

bedtools makewindows -g "$CONTIGS" -w "$WINDOW_SIZE" > "$WINDOWS"
bcftools query -l "$VCF" > "$SAMPLES"

: > "$COUNTS_LIST"

while read -r SAMPLE; do
  SAMPLE_BED="$WORKDIR/${SAMPLE}.bed"
  SAMPLE_COUNTS="$WORKDIR/${SAMPLE}.counts"

  bcftools query -s "$SAMPLE" -f '%CHROM\t%POS\t[%GT]\n' "$VCF" \
    | awk '($3 ~ /(^|[\/|])[1-9]/) {print $1, $2-1, $2}' OFS="\t" \
    > "$SAMPLE_BED"

  bedtools coverage -a "$WINDOWS" -b "$SAMPLE_BED" -counts | cut -f4 > "$SAMPLE_COUNTS"
  echo "$SAMPLE_COUNTS" >> "$COUNTS_LIST"
done < "$SAMPLES"

{
  printf "Chrom\tStart\tEnd"
  while read -r SAMPLE; do printf "\t%s" "$SAMPLE"; done < "$SAMPLES"
  printf "\n"
} > "$OUT"

paste "$WINDOWS" $(cat "$COUNTS_LIST") >> "$OUT"

echo "Done: $OUT"
