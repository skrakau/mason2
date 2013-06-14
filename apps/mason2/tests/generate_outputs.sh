#!/bin/sh
#
# We create variants from a randomly generated genome.

VARIATOR=../../../../../../seqan-trunk-build/debug/bin/mason_variator
MATERIALIZER=../../../../../../seqan-trunk-build/debug/bin/mason_materializer

# ============================================================
# mason_variator
# ============================================================

echo "${VARIATOR} -n 2 -if random.fasta -ov random_var1.vcf -of random_var1.fasta --snp-rate 0.001 --small-indel-rate 0.001 --sv-indel-rate 0.001 --sv-inversion-rate 0.001 --sv-translocation-rate 0.001 --sv-duplication-rate 0.001 --min-sv-size 50 --max-sv-size 100 --methylation-levels --meth-fasta-out random_var1_meth.fasta --out-breakpoints random_var1_bp.txt >random_var1.vcf.stdout 2>random_var1.vcf.stderr"
${VARIATOR} -n 2 -if random.fasta -ov random_var1.vcf -of random_var1.fasta --snp-rate 0.001 --small-indel-rate 0.001 --sv-indel-rate 0.001 --sv-inversion-rate 0.001 --sv-translocation-rate 0.001 --sv-duplication-rate 0.001 --min-sv-size 50 --max-sv-size 100 --methylation-levels --meth-fasta-out random_var1_meth.fasta --out-breakpoints random_var1_bp.txt >random_var1.vcf.stdout 2>random_var1.vcf.stderr

# ============================================================
# mason_materializer
# ============================================================

echo "${MATERIALIZER} -if random.fasta -iv random_var1.vcf -of materializer.random_var1.fasta >materializer.random_var1.stdout 2>materializer.random_var1.stderr"
${MATERIALIZER} -if random.fasta -iv random_var1.vcf -of materializer.random_var1.fasta >materializer.random_var1.stdout 2>materializer.random_var1.stderr
rm materializer.random_var1.fasta  # we'll compare against variator output
