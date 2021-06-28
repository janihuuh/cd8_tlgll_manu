# File access
#sudo chown -R pipeuser vartrix/

# Link files
#find bam -name "*possorted_genome_bam.bam" | perl -p -e 's/^.*/ln -s ..\/$& xxx $&/g' | perl -p -e 's/xxx.bam./analysis_bam\//g' | perl -p -e 's|(.+)\/|$1_|' | sh
#find bam -name "*possorted_genome_bam.bam.bai" | perl -p -e 's/^.*/ln -s ..\/$& xxx $&/g' | perl -p -e 's/xxx.bam./analysis_bam\//g' | perl -p -e 's|(.+)\/|$1_|' | sh
#
#find barcode -name "*barcodes.tsv.gz" | perl -p -e 's/^.*/ln -s ..\/$& xxx $&/g' | perl -p -e 's/xxx.barcode./analysis_barcode\//g' | perl -p -e 's|(.+)\/|$1_|' | sh
#gzip -d --force *

# Set directories
IN_FOLD="/csc/mustjoki/vartrix/t_lgll/"

FILES=${IN_FOLD}/analysis_bam/*bam
for SAMPLE in $FILES
do
  OUT=${SAMPLE}
  OUT=${OUT##*/}
  OUT=${OUT}
  OUT=${OUT/_possorted_genome_bam.bam/}

  if [ ! -f ${OUT}.alt.out ]; then
    echo grun.py -n vartrix-${OUT} -q seq_hugemem2.q -c \"/fs/vault/pipelines/gatk/bin/vartrix-1.1.0/vartrix-v1.1.0-x86_64-linux/vartrix --bam ${SAMPLE} --cell-barcodes analysis_barcode/${OUT}_barcodes.tsv --fasta /csc/mustjoki2/bioinformatics/transcriptome/data/homo_sapiens_v82_star/homo_sapiens_82.fa --out-matrix ${OUT}.alt.out --ref-matrix ${OUT}.ref.out --threads 12 --vcf ../../scrnaseq_vartrix/LGL/cosmic/cosmic_v86.len10.nonN.nodot.vcf --out-variants ${OUT}.vars --scoring-method coverage --umi\"

    #echo /fs/vault/pipelines/gatk/bin/vartrix-1.1.0/vartrix-v1.1.0-x86_64-linux/vartrix --bam ${SAMPLE} --cell-barcodes analysis_barcode/${OUT}_barcodes.tsv --fasta /csc/mustjoki2/bioinformatics/transcriptome/data/homo_sapiens_v82_star/homo_sapiens_82.fa --out-matrix ${OUT}.alt.out --ref-matrix ${OUT}.ref.out --threads 12 --vcf ../../scrnaseq_vartrix/LGL/cosmic/cosmic_v86.len10.nonN.nodot.vcf --out-variants ${OUT}.vars --scoring-method coverage --umi
  fi

done

#/fs/vault/pipelines/gatk/bin/vartrix-1.1.0/vartrix-v1.1.0-x86_64-linux/vartrix --bam ./bam/SI-GA-D6_outs_possorted_genome_bam.bam --cell-barcodes barcode/SI-GA-D6_outs_possorted_genome_bam_barcodes.tsv --fasta /fs/vault/pipelines/rnaseq/data/homo_sapiens_v82_star/homo_sapiens_82.fa --out-matrix SI-GA-D6_outs_possorted_genome_bam.alt.out --ref-matrix SI-GA-D6_outs_possorted_genome_bam.ref.out --threads 12 --vcf cosmic/cosmic_v86.len10.nonN.nodot.vcf --out-variants SI-GA-D6_outs_possorted_genome_bam.vars --scoring-method coverage --umi

# Split BAM by matching ref
/csc/mustjoki2/bioinformatics/general/bin/bamtools/bin/bamtools split -in ../analysis_bam/HRUH1248_CD3neg_possorted_genome_bam.bam -reference
/csc/mustjoki2/bioinformatics/general/bin/bamtools/bin/bamtools split -in ../analysis_bam/HRUH1248_CD3pos_possorted_genome_bam.bam -reference;

# move REF split BAMs
mkdir bam_split
cd bam_split
mv ../analysis_bam/*REF*

# Index all BAMs
ls *.bam | xargs -n1 -P5 samtools index
