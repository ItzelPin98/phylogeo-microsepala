#!/bin/bash

#---------------------------------------
# Script para correr ADMIXTURE excluyendo outgroups (miss30)
#---------------------------------------

# Programas
VCFTOOLS=/opt/homebrew/bin/vcftools
PLINK=/opt/homebrew/bin/plink
ADMIXTURE=/Applications/admixture_macosx-1.3.0/admixture

# Rutas de entrada y salida
VCF="/Users/itziliwi/Desktop/admixture_test_20250710/01_datos/microv3_DP20_MAF002_miss30.recode.vcf"
OUTDIR="/Users/itziliwi/Desktop/admixture_test_20250710/02_out/admixture_out_filtered_20250717"
PREFIX="microv3_allSamples_filtered"

# Crear carpeta de salida
mkdir -p "$OUTDIR"

# Ruta al archivo de outgroups (FID IID por línea)
OUTGROUPS_FILE="/Users/itziliwi/Desktop/admixture_test_20250710/01_datos/outgroups_deleteAdmixture.txt"

# 1. Convertir VCF a formato PED/MAP (PLINK)
$VCFTOOLS \
  --vcf "$VCF" \
  --plink \
  --out "$OUTDIR/$PREFIX"

# 2. Convertir PED/MAP a BED/BIM/FAM
$PLINK \
  --file "$OUTDIR/$PREFIX" \
  --make-bed \
  --out "$OUTDIR/${PREFIX}_tmp"

# 3. Remover outgroups
$PLINK \
  --bfile "$OUTDIR/${PREFIX}_tmp" \
  --remove "$OUTGROUPS_FILE" \
  --make-bed \
  --out "$OUTDIR/$PREFIX"

# 4. Ejecutar ADMIXTURE para K=1 a K=10 sin réplicas
cd "$OUTDIR"
for K in {1..20}
do
    $ADMIXTURE --cv ${PREFIX}.bed $K | tee log_K${K}.out
done

# 5. Extraer errores CV
grep -h "CV error" log*.out > CV_errors.txt
echo "¡Corridas finalizadas! Resultados en $OUTDIR"
