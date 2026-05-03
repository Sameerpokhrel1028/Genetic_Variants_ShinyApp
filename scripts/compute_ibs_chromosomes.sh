#!/usr/bin/env bash
set -euo pipefail

VCF=""
OUT=""
MAX_MISS="0.20"
REGION_MB=10
PLINK_MEM="16000"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--vcf) VCF="$2"; shift 2 ;;
    -o|--out) OUT="$2"; shift 2 ;;
    -m|--max-miss) MAX_MISS="$2"; shift 2 ;;
    --region-mb) REGION_MB="$2"; shift 2 ;;
    --memory) PLINK_MEM="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: bash compute_ibs_scope.sh --vcf input.vcf.gz --out ibs.scope.tsv [--region-mb 10]"
      exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

[[ -z "$VCF" || -z "$OUT" ]] && { echo "ERROR: --vcf and --out required"; exit 1; }

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

CLEAN0="$TMPDIR/clean_snps.biallelic.vcf.gz"
CLEAN="$TMPDIR/clean_snps.filtered.vcf.gz"

bcftools view -m2 -M2 -v snps -Oz -o "$CLEAN0" "$VCF"
bcftools index -t "$CLEAN0"

bcftools view -i "F_MISSING<$MAX_MISS" -Oz -o "$CLEAN" "$CLEAN0"
bcftools index -t "$CLEAN"

echo -e "Scope\tSample1\tSample2\tDistance" > "$OUT"

compute_scope() {
  local SCOPE="$1"
  local REGION="$2"

  local VCF_USE="$CLEAN"
  local SCOPE_VCF="$TMPDIR/${SCOPE}.vcf.gz"
  local PLAIN="$TMPDIR/${SCOPE}.vcf"
  local PFX="$TMPDIR/plink_${SCOPE}"
  local DFX="$TMPDIR/dist_${SCOPE}"

  if [[ -n "$REGION" ]]; then
    bcftools view -r "$REGION" -Oz -o "$SCOPE_VCF" "$CLEAN" || return 0
    bcftools index -t "$SCOPE_VCF"
    [[ "$(bcftools view -H "$SCOPE_VCF" | wc -l)" -eq 0 ]] && return 0
    VCF_USE="$SCOPE_VCF"
  fi

  gunzip -c "$VCF_USE" > "$PLAIN"

  plink --vcf "$PLAIN" --double-id --allow-extra-chr \
    --make-bed --out "$PFX" --memory "$PLINK_MEM" >/dev/null

  [[ ! -s "${PFX}.bim" ]] && return 0

  plink --bfile "$PFX" --allow-extra-chr \
    --distance 1-ibs square --out "$DFX" --memory "$PLINK_MEM" >/dev/null

  MAT=""
  IDF=""
  MODE=""

  if [[ -s "${DFX}.dist" ]]; then
    MAT="${DFX}.dist"; IDF="${DFX}.dist.id"; MODE="dist"
  elif [[ -s "${DFX}.mdist" ]]; then
    MAT="${DFX}.mdist"; IDF="${DFX}.mdist.id"; MODE="dist"
  elif [[ -s "${DFX}.mibs" ]]; then
    MAT="${DFX}.mibs"; IDF="${DFX}.mibs.id"; MODE="mibs"
  else
    echo "WARNING: no PLINK distance output for $SCOPE" >&2
    return 0
  fi

  python3 - <<PY
import pandas as pd
import numpy as np

scope = "${SCOPE}"
mode = "${MODE}"
mat = pd.read_csv("${MAT}", sep=r"\s+", header=None)
ids = pd.read_csv("${IDF}", sep=r"\s+", header=None)
names = ids.iloc[:,1] if ids.shape[1] > 1 else ids.iloc[:,0]
names = names.astype(str).tolist()

if mode == "mibs":
    mat = 1.0 - mat

mat = mat.clip(lower=0)
np.fill_diagonal(mat.values, 0.0)

with open("${OUT}", "a") as f:
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            f.write(f"{scope}\t{names[i]}\t{names[j]}\t{mat.iat[i,j]}\n")
PY
}

compute_scope "Genome" ""

CHRS=$(bcftools view -h "$CLEAN" | grep '^##contig' | sed -E 's/.*ID=([^,]+).*/\1/')

for CHR in $CHRS; do
  LEN=$(bcftools view -h "$CLEAN" | awk -v c="$CHR" -F'[=,<>]' '
  /^##contig=/{
    id=""; len="";
    for(i=1;i<=NF;i++){
      if($i=="ID") id=$(i+1);
      if($i=="length") len=$(i+1);
    }
    if(id==c) print len;
  }')

  compute_scope "$CHR" "$CHR"

  STEP=$(( REGION_MB * 1000000 ))
  START=1
  while [[ "$START" -le "$LEN" ]]; do
    END=$(( START + STEP - 1 ))
    [[ "$END" -gt "$LEN" ]] && END="$LEN"

    SCOPE="${CHR}:${START}-${END}"
    REGION="${CHR}:${START}-${END}"
    compute_scope "$SCOPE" "$REGION"

    START=$(( END + 1 ))
  done
done

echo "Done: $OUT"
