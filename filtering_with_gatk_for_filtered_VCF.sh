
#create an index file using samtools

brew install samtools

samtools faidx ./lyrata.fasta.fa

#create fasta dictionary file
#to do this need to download picard 

#https://github.com/broadinstitute/picard/releases

mv picard.jar gatk-4.5.0.0

java -jar ./picard.jar SortVcf \
INPUT=Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz \
OUTPUT=sorted.vcf.gz



#now run code
gatk SelectVariants -R lyrata.fasta.fa -V sorted.vcf.gz -sn KAG-01tl -sn KAG-02tl -sn KAG-03tl -sn KAG-04tl -sn KAG-05tl -sn LIC-01tl -sn LIC-02tl -sn LIC-03tl -sn MOD-01tl -sn MOD-02tl -sn MOD-03tl -sn MOD-04tl -sn FRE-01tl -sn FRE-02tl -sn FRE-03tl -sn FRE-04tl -sn FRE-05tl -sn FRE-06tl -sn FRE-07tl -sn FRE-08tl -sn OCH-01tl -sn OCH-02tl -sn OCH-03tl -sn OCH-04tl -sn OCH-05tl -sn OCH-06tl -sn OCH-07tl -sn OCH-08tl -sn HAB-01tl -sn HAB-02tl -sn HAB-03tl -sn PIZ-03dl -sn PIZ-04dl -sn PIZ-05dl -sn PIZ-06dl -sn PIZ-08dl -sn PIZ-09dl -sn PIZ-11dl -sn KEH-01tl -sn KEH-04tl -sn KEH-05tl -sn KEH-06tl -sn KEH-07tl -sn KEH-08tl -sn KEH-09tl -sn KEH-10tl -sn BZD-01tl -sn BZD-02tl -sn BZD-03tl -sn BZD-04tl --allow-nonoverlapping-command-line-samples -O Filtered_VCF.vcf
