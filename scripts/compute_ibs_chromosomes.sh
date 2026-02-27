#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Compute IBS-derived distance (1-IBS) genome-wide and per chromosome, write ONE long TSV.

Required:
  -v, --vcf       Input multi-sample VCF (.vcf.gz recommended)
  -o, --out       Output long-format TSV

Optional:
  -m, --max-miss  Max missingness filter using F_MISSING [default: 0.20]
  --memory        PLINK memory in MB [default: 16000]

Output TSV columns:
  Scope   Sample1   Sample2   Distance

Scope values:
  Genome and all contigs in VCF header (e.g., chr01..chr20)

Example:
  bash scripts/compute_ibs_scope.sh \
    --vcf /scratch/.../pangenome_v3.vcf.gz \
    --out ibs_scope.tsv
EOF
}

VCF=""
OUT=""
MAX_MISS="0.20"
PLINK_MEM="16000"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--vcf) VCF="$2"; shift 2 ;;
    -o|--out) OUT="$2"; shift 2 ;;
    -m|--max-miss) MAX_MISS="$2"; shift 2 ;;
    --memory) PLINK_MEM="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$VCF" || -z "$OUT" ]]; then
  echo "ERROR: --vcf and --out are required." >&2
  usage
  exit 1
fi

for cmd in bcftools plink python3 gunzip; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found in PATH." >&2; exit 1; }
done

if [[ ! -f "$VCF" ]]; then
  echo "ERROR: VCF not found: $VCF" >&2
  exit 1
fi

# --------------------------------------
# Check for risky sample names
# --------------------------------------
echo "Checking sample names..."
if bcftools query -l "$VCF" | grep -qE '[[:space:]/:()\[\]]'; then
  echo "WARNING: Sample names contain spaces or special characters (/ : ( ) [ ]). PLINK may fail."
  echo "Tip: sanitize names with bcftools reheader before running."
fi
echo "--------------------------------------"

TMPDIR=$(mktemp -d)

cleanup() {
  if [[ $? -eq 0 ]]; then
    rm -rf "$TMPDIR"
  else
    echo "" >&2
    echo "ERROR occurred. Temp kept at: $TMPDIR" >&2
    echo "Logs to inspect:" >&2
    echo "  $TMPDIR/*plink*.log" >&2
    echo "  $TMPDIR/*ibs*.log" >&2
    echo "" >&2
  fi
}
trap cleanup EXIT

CLEAN0="$TMPDIR/clean_snps.biallelic.vcf.gz"
CLEAN="$TMPDIR/clean_snps.filtered.vcf.gz"

echo "--------------------------------------"
echo "Variant Summary (raw VCF)"
echo "--------------------------------------"
echo "Total variants:        $(bcftools view -H "$VCF" | wc -l)"
echo "Total SNPs:            $(bcftools view -H -v snps "$VCF" | wc -l)"
echo "Total INDELs:          $(bcftools view -H -v indels "$VCF" | wc -l)"
echo "Multiallelic variants: $(bcftools view -H -m3 "$VCF" | wc -l)"
echo "Biallelic SNPs:        $(bcftools view -H -m2 -M2 -v snps "$VCF" | wc -l)"
echo "--------------------------------------"

echo "[1/4] Filtering clean SNP set (biallelic SNPs + missingness < $MAX_MISS)..."
bcftools view -m2 -M2 -v snps -Oz -o "$CLEAN0" "$VCF"
bcftools index -t "$CLEAN0"

bcftools view -i "F_MISSING<$MAX_MISS" -Oz -o "$CLEAN" "$CLEAN0"
bcftools index -t "$CLEAN"

KEPT=$(bcftools view -H "$CLEAN" | wc -l)
echo "SNPs kept after missingness filter: $KEPT"
if [[ "$KEPT" -eq 0 ]]; then
  echo "ERROR: 0 SNPs remain after filtering. Try a higher --max-miss (e.g. 0.5) or check VCF." >&2
  exit 1
fi
echo "--------------------------------------"

echo -e "Scope\tSample1\tSample2\tDistance" > "$OUT"

compute_scope() {
  local SCOPE="$1"
  local REGION="$2"

  local VCF_USE="$CLEAN"
  local SCOPE_VCFGZ="$TMPDIR/scope_${SCOPE}.vcf.gz"
  local INVCF="$TMPDIR/input_${SCOPE}.vcf"

  # Unique prefixes per scope (avoid collisions)
  local PLINK_PREFIX="$TMPDIR/plink_${SCOPE}"
  local PLINK_OUT="$TMPDIR/ibs_${SCOPE}"

  # Make scope VCF if region requested
  if [[ -n "$REGION" ]]; then
    bcftools view -r "$REGION" -Oz -o "$SCOPE_VCFGZ" "$CLEAN"
    bcftools index -t "$SCOPE_VCFGZ"
    VCF_USE="$SCOPE_VCFGZ"

    # Check scope has variants
    local N
    N=$(bcftools view -H "$VCF_USE" | wc -l)
    if [[ "$N" -eq 0 ]]; then
      echo "  NOTE: $SCOPE has 0 variants after filtering; skipping." >&2
      return 0
    fi
  fi

  # PLINK 1.9 often behaves better with plain VCF than .vcf.gz
  gunzip -c "$VCF_USE" > "$INVCF"

  echo "    PLINK import ($SCOPE)..."
  plink --vcf "$INVCF" \
    --double-id --allow-extra-chr \
    --make-bed --out "$PLINK_PREFIX" \
    --memory "$PLINK_MEM"

  if [[ ! -s "${PLINK_PREFIX}.bim" ]]; then
    echo "ERROR: PLINK produced empty .bim (no variants imported) for scope: $SCOPE" >&2
    echo "Check log: ${PLINK_PREFIX}.log" >&2
    exit 1
  fi

  echo "    PLINK distance ($SCOPE)..."
  plink --bfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --distance 1-ibs square \
    --out "$PLINK_OUT" \
    --memory "$PLINK_MEM"

  # Determine which files PLINK wrote (.dist OR .mdist OR .mibs)
  local MAT=""
  local IDF=""
  local MODE=""   # "dist"/"mdist" already distance, "mibs" similarity

  if [[ -s "${PLINK_OUT}.dist" && -s "${PLINK_OUT}.dist.id" ]]; then
    MAT="${PLINK_OUT}.dist"
    IDF="${PLINK_OUT}.dist.id"
    MODE="dist"
  elif [[ -s "${PLINK_OUT}.mdist" && -s "${PLINK_OUT}.mdist.id" ]]; then
    MAT="${PLINK_OUT}.mdist"
    IDF="${PLINK_OUT}.mdist.id"
    MODE="mdist"
  elif [[ -s "${PLINK_OUT}.mibs" && -s "${PLINK_OUT}.mibs.id" ]]; then
    MAT="${PLINK_OUT}.mibs"
    IDF="${PLINK_OUT}.mibs.id"
    MODE="mibs"
  else
    echo "ERROR: PLINK did not produce distance outputs for scope: $SCOPE" >&2
    echo "Checked for:" >&2
    echo "  ${PLINK_OUT}.dist  / ${PLINK_OUT}.dist.id" >&2
    echo "  ${PLINK_OUT}.mdist / ${PLINK_OUT}.mdist.id" >&2
    echo "  ${PLINK_OUT}.mibs  / ${PLINK_OUT}.mibs.id" >&2
    echo "See logs:" >&2
    echo "  ${PLINK_PREFIX}.log" >&2
    echo "  ${PLINK_OUT}.log" >&2
    exit 1
  fi

  # Append upper triangle to OUT in long format
  python3 - <<PY
import pandas as pd
import numpy as np

scope = "${SCOPE}"
mode  = "${MODE}"

mat = pd.read_csv("${MAT}", sep=r"\\s+", header=None)
ids = pd.read_csv("${IDF}", sep=r"\\s+", header=None)

names = ids.iloc[:,1] if ids.shape[1] > 1 else ids.iloc[:,0]
names = names.astype(str).tolist()

# Only mibs is similarity; convert to distance
if mode == "mibs":
    mat = 1.0 - mat

mat = mat.clip(lower=0)
np.fill_diagonal(mat.values, 0.0)

with open("${OUT}", "a") as f:
    n = len(names)
    for i in range(n):
        for j in range(i+1, n):
            f.write(f"{scope}\\t{names[i]}\\t{names[j]}\\t{mat.iat[i,j]}\\n")
PY
}

echo "[2/4] Genome-wide distance (1-IBS)..."
compute_scope "Genome" ""

echo "[3/4] Chromosome-wise distance (1-IBS)..."
CHRS=$(bcftools view -h "$CLEAN" | grep '^##contig' | sed -E 's/.*ID=([^,]+).*/\1/')

for CHR in $CHRS; do
  echo "  Processing $CHR"
  compute_scope "$CHR" "$CHR"
done

echo "--------------------------------------"
echo "Finished."
echo "Output TSV: $OUT"
echo "--------------------------------------"
