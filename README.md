Usage
```
usage: SNP_pN_pS_calculate.v4.py [-h] -f FFN -b VCF -o OUTDIR [-d MINDEPTH]
                                 [-c SNPCOVERAGE]

optional arguments:
  -h, --help            show this help message and exit
  -f FFN, --ffn FFN     ffn path from Prodigal
  -b VCF, --vcf VCF     vcf file from bcftools
  -o OUTDIR, --outdir OUTDIR
                        name of outdir
  -d MINDEPTH, --minDepth MINDEPTH
                        minimum depth, default 20
  -c SNPCOVERAGE, --SNPcoverage SNPCOVERAGE
                        SNP coverage, default 0.05
```

Example
The sorted bam file is obtained by mapping the reads to genome.
```
samtools faidx YB4_M4.fa
bcftools mpileup --threads 8 --annotate DP4 -d 20 -q 20 -Q 30 -O b -f YB4_M4.fa YB4_M4.sort.bam -o YB4_M4.bcf.gz
bcftools filter -e 'TYPE=="INDEL" || DP < 20 || (DP4[*:0]+DP4[*:1]+DP4[*:2]+DP4[*:3] < 20)' --threads 6 -g 10 YB4_M4.bcf.gz | bcftools view -H > YB4_M4.filter.vcf
## add bin names for vcf file
sh modified_vcf.sh YB4_M4.filter.vcf YB4_M4 > YB4_M4.filter.vcf.f
## predicte the CDS
mkdir prodigal_out
prodigal -i YB4_M4.fa -f gff -a prodigal_out/YB4_M4.faa -o prodigal_out/YB4_M4.gff -d prodigal_out/YB4_M4.ffn
## calculate the PNPS
SNP_pN_pS_calculate.v4.py -f Hybrid_Normal_MAGs.prodigal -b YB4_M4.filter.vcf.f -o pnps
```
