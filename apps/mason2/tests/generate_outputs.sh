#!/bin/sh
#
# We create variants from a randomly generated genome.

GENOME=../../../../../../seqan-trunk-build/debug/bin/mason_genome
VARIATOR=../../../../../../seqan-trunk-build/debug/bin/mason_variator
MATERIALIZER=../../../../../../seqan-trunk-build/debug/bin/mason_materializer
SIMULATOR=../../../../../../seqan-trunk-build/debug/bin/mason_simulator

# ============================================================
# mason_genome
# ============================================================

echo "${GENOME} -l 1000 -o genome.test1.fasta >genome.test1.stdout 2>genome.test1.stderr"
${GENOME} -l 1000 -o genome.test1.fasta >genome.test1.stdout 2>genome.test1.stderr
echo "${GENOME} -s 1 -l 1000 -l 100 -o genome.test2.fasta >genome.test2.stdout 2>genome.test2.stderr"
${GENOME} -s 1 -l 1000 -l 100 -o genome.test2.fasta >genome.test2.stdout 2>genome.test2.stderr

# ============================================================
# mason_variator
# ============================================================

echo "${VARIATOR} -n 2 -ir random.fasta -ov random_var1.vcf -of random_var1.fasta --snp-rate 0.001 --small-indel-rate 0.001 --sv-indel-rate 0.001 --sv-inversion-rate 0.001 --sv-translocation-rate 0.001 --sv-duplication-rate 0.001 --min-sv-size 50 --max-sv-size 100 --methylation-levels --meth-fasta-out random_var1_meth.fasta --out-breakpoints random_var1_bp.txt >random_var1.vcf.stdout 2>random_var1.vcf.stderr"
${VARIATOR} -n 2 -ir random.fasta -ov random_var1.vcf -of random_var1.fasta --snp-rate 0.001 --small-indel-rate 0.001 --sv-indel-rate 0.001 --sv-inversion-rate 0.001 --sv-translocation-rate 0.001 --sv-duplication-rate 0.001 --min-sv-size 50 --max-sv-size 100 --methylation-levels --meth-fasta-out random_var1_meth.fasta --out-breakpoints random_var1_bp.txt >random_var1.vcf.stdout 2>random_var1.vcf.stderr

# ============================================================
# mason_materializer
# ============================================================

echo "${MATERIALIZER} -ir random.fasta -iv random_var1.vcf -o materializer.random_var1.fasta >materializer.random_var1.stdout 2>materializer.random_var1.stderr"
${MATERIALIZER} -ir random.fasta -iv random_var1.vcf -o materializer.random_var1.fasta >materializer.random_var1.stdout 2>materializer.random_var1.stderr
rm -f materializer.random_var1.fasta  # we'll compare against variator output

# ============================================================
# mason_simulator
# ============================================================

# Without VCF variants, FASTQ output, with SAM alignments, paired-end
echo "${SIMULATOR} -n 1000 -ir random.fasta -o simulator.left1.fq -or simulator.right1.fq -oa simulator.out1.sam >simulator.out1.stdout 2>simulator.out1.stderr"
${SIMULATOR} -n 1000 -ir random.fasta -o simulator.left1.fq -or simulator.right1.fq -oa simulator.out1.sam >simulator.out1.stdout 2>simulator.out1.stderr
# With VCF variants, FASTQ output, with SAM alignment, paired-end
echo "${SIMULATOR} -n 1000 -ir random.fasta -iv random_var1.vcf -o simulator.left2.fq -or simulator.right2.fq -oa simulator.out2.sam >simulator.out2.stdout 2>simulator.out2.stderr"
${SIMULATOR} -n 1000 -ir random.fasta -iv random_var1.vcf -o simulator.left2.fq -or simulator.right2.fq -oa simulator.out2.sam >simulator.out2.stdout 2>simulator.out2.stderr
# Without VCF variants, FASTA output, no SAM alignments, paired-end
echo "${SIMULATOR} -n 1000 -ir random.fasta -o simulator.left3.fa -or simulator.right3.fa >simulator.out3.stdout 2>simulator.out3.stderr"
${SIMULATOR} -n 1000 -ir random.fasta -o simulator.left3.fa -or simulator.right3.fa >simulator.out3.stdout 2>simulator.out3.stderr
# Without VCF variants, FASTA output, no SAM alignments, single-end
echo "${SIMULATOR} -n 1000 -ir random.fasta -o simulator.left4.fa >simulator.out4.stdout 2>simulator.out4.stderr"
${SIMULATOR} -n 1000 -ir random.fasta -o simulator.left4.fa >simulator.out4.stdout 2>simulator.out4.stderr
