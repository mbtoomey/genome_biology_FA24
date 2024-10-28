vcftools --gzvcf flicker_chr3.phased.HZ.recode.vcf --plink --out flicker_chr3_outputPlinkformat

plink --file flicker_chr3_outputPlinkformat --make-bed --out flicker_chr3_HZ_bed