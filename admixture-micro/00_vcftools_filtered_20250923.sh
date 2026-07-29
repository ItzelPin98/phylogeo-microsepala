#!/bin/bash

# Definir rutas base
VCF_DIR="/home/ivalu/RADseq_micro/micro_test_20250604/micro_c85_m95_outfiles"
OUT_DIR="/home/ivalu/RADseq_micro/micro_test_20250604/micro_c85_m95_outfiles/vcftools_filtered_20250923"

# Crear carpeta de salida si no existe
mkdir -p "$OUT_DIR"

# 0) Filtrar solo SNPs
vcftools \
  --vcf "$VCF_DIR/micro_c85_m95.vcf" \
  --remove-indels \
  --recode \
  --out "$OUT_DIR/micro0923_SNPs"

# 1) Transformar de VCF a 012 (sin filtrar aún)
vcftools \
  --vcf "$OUT_DIR/micro0923_SNPs.recode.vcf" \
  --012 \
  --out "$OUT_DIR/micro0923_012"

# 2) Filtrar por profundidad mínima promedio >=20 y MAF >=0.02 y máx. 30% datos faltantes
vcftools \
  --vcf "$OUT_DIR/micro0923_SNPs.recode.vcf" \
  --min-meanDP 21 \
  --maf 0.02 \
  --max-missing 0.7 \
  --recode \
  --out "$OUT_DIR/microv0923_DP20_MAF002_miss30"

# 3) Calcular proporción de datos faltantes por individuo
vcftools \
  --vcf "$OUT_DIR/microv0923_DP20_MAF002_miss30.recode.vcf" \
  --missing-indv \
  --out "$OUT_DIR/micro0923_missing"

# 4) Crear archivo con individuos con >10% de datos faltantes
awk '$5 > 0.10 {print $1}' "$OUT_DIR/micro0923_missing.imiss" > "$VCF_DIR/ind_missHigh_micro0923.txt"

# 5) Quitar individuos con >10% de datos faltantes
vcftools \
  --vcf "$OUT_DIR/microv0923_DP20_MAF002_miss30.recode.vcf" \
  --remove "$VCF_DIR/ind_missHigh_micro0923.txt" \
  --recode \
  --out "$OUT_DIR/micro0923_miss10_filtered"

# 6) Transformar en matriz 012 final
vcftools \
  --vcf "$OUT_DIR/micro0923_miss10_filtered.recode.vcf" \
  --012 \
  --out "$OUT_DIR/micro0923_final012_filtered"
  

