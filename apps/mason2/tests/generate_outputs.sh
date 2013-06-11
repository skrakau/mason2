#!/bin/sh
#
# We create variants from a randomly generated genome.

VARIATOR=../../../../../../seqan-trunk-build/debug-clang/bin/mason_variator

# ============================================================
# Simulate Variants
# ============================================================

echo "${VARIATOR} -n 2 -if random.fasta -ov random_var1.vcf -of random_var2.vcf --snp-rate 0.001 --small-indel-rate 0.001 --sv-indel-rate 0.001 --sv-inversion-rate 0.001 --sv-translocation-rate 0.001 --sv-duplication-rate 0.001 --min-sv-len 50 --max-sv-len 100 >random_var1.vcf.stdout 2>random_var1.vcf.stderr"
${VARIATOR} -n 2 -if random.fasta -ov random_var1.vcf -of random_var1.fasta --snp-rate 0.001 --small-indel-rate 0.001 --sv-indel-rate 0.001 --sv-inversion-rate 0.001 --sv-translocation-rate 0.001 --sv-duplication-rate 0.001 --min-sv-size 50 --max-sv-size 100 >random_var1.vcf.stdout 2>random_var1.vcf.stderr
