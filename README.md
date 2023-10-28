# EP2_AulaSomatica

# 1. Instalando sratoolskit

`` 
brew install sratoolkit
``

# 2. Parallel fastq-dump

``
pip install parallel-fastq-dump
``

# 3. Não deu certo? Fazer isso antes:

``wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz
``

``tar -zxvf sratoolkit.3.0.0-ubuntu64.tar.gz
``

``export PATH=$PATH://workspace/somaticoEP1/sratoolkit.3.0.0-ubuntu64/bin/
``

``echo "Aexyo" | sratoolkit.3.0.0-ubuntu64/bin/vdb-config -i
``

``time parallel-fastq-dump --sra-id SRR8856724 \
--threads 4 \
--outdir ./ \
--split-files \
--gzip``

# 4. Referências do Genoma hg19 

``wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf``

``wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx``

``wget -c  https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf``

``wget -c  https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx``

*Arquivo no formato FASTA do genoma humano hg19*

``wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz``

# 5. BWA 

``brew install bwa``

*BWA index do chr9.fa.gz*

``gunzip chr9.fa.gz``

``bwa index chr9.fa``

``brew install samtools``

``samtools faidx chr9.fa``

*Combinar bwa + samtools view e sort*

``NOME=WP312; Biblioteca=Nextera; Plataforma=illumina;
bwa mem -t 10 -M -R "@RG\tID:$NOME\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" chr9.fa SRR8856724_1.fastq.gz SRR8856724_2.fastq.gz | samtools view -F4 -Sbu -@2 - | samtools sort -m4G -@2 -o WP312_sorted.bam``

*Remover duplicata de PCR*

``samtools rmdup WP312_sorted.bam WP312_sorted_rmdup.bam``

# 6. Instalar bedtools

``brew install bedtools``

*Gerar BED do arquivo BAM*

``bedtools bamtobed -i WP312_sorted_rmdup.bam > WP312_sorted_rmdup.bed``

``bedtools merge -i WP312_sorted_rmdup.bed > WP312_sorted_rmdup_merged.bed``

``bedtools sort -i WP312_sorted_rmdup_merged.bed > WP312_sorted_rmdup_merged_sorted.bed``

*Cobertura Média*

``bedtools coverage -a WP312_sorted_rmdup_merged_sorted.bed \
-b WP312_sorted_rmdup.bam -mean \
> WP312_coverageBed.bed``

*Filtro por total de reads >=20*

``cat WP312_coverageBed.bed | \
awk -F "\t" '{if($4>=20){print}}' \
> WP312_coverageBed20x.bed``

# 7. Instalação do GATK

``wget -c https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip``

*Descompactar*

``unzip gatk-4.2.2.0.zip``

*Gerar arquivo .dict*
``./gatk-4.2.2.0/gatk CreateSequenceDictionary -R chr9.fa -O chr9.dict``

*Gerar interval_list do chr9*

``./gatk-4.2.2.0/gatk ScatterIntervalsByNs -R chr9.fa -O chr9.interval_list -OT ACGT``

*Converter Bed para Interval_list*

``./gatk-4.2.2.0/gatk BedToIntervalList -I WP312_coverageBed20x.bed \
-O WP312_coverageBed20x.interval_list -SD chr9.dict``

# 8. Adicionar chr nos VCFs do Gnomad e PoN

``grep "\#" af-only-gnomad.raw.sites.vcf > af-only-gnomad.raw.sites.chr.vcf
grep  "^9" af-only-gnomad.raw.sites.vcf |  awk '{print("chr"$0)}' >> af-only-gnomad.raw.sites.chr.vcf``

*Indexing*

``bgzip af-only-gnomad.raw.sites.chr.vcf
tabix -p vcf af-only-gnomad.raw.sites.chr.vcf.gz``


 ``grep "\#" Mutect2-WGS-panel-b37.vcf > Mutect2-WGS-panel-b37.chr.vcf 
grep  "^9" Mutect2-WGS-panel-b37.vcf |  awk '{print("chr"$0)}' >> Mutect2-WGS-panel-b37.chr.vcf ``

``bgzip Mutect2-WGS-panel-b37.chr.vcf 
tabix -p vcf Mutect2-WGS-panel-b37.chr.vcf.gz``

# 9. GATK4

``./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I WP312_sorted_rmdup.bam  \
	-V af-only-gnomad.raw.sites.chr.vcf.gz  \
	-L WP312_coverageBed20x.interval_list \
	-O WP312.table``

 *CalculateContamination*

 ``./gatk-4.2.2.0/gatk CalculateContamination \
-I WP312.table \
-O WP312.contamination.table``

*GATK4 - MuTect2 Call*

``./gatk-4.2.2.0/gatk Mutect2 \
  -R chr9.fa \
  -I WP312_sorted_rmdup_F4.bam \
  --germline-resource af-only-gnomad.hg38.vcf.gz  \
  --panel-of-normals 1000g_pon.hg38.vcf.gz \
  -L WP312_coverageBed20x.interval_list \
  -O WP312.somatic.pon.vcf.gz``

*GATK4 - MuTect2 FilterMutectCalls*

  ``./gatk-4.2.2.0/gatk FilterMutectCalls \
-R chr9.fa \
-V WP312.somatic.pon.vcf.gz \
--contamination-table WP312.contamination.table \
-O WP312.filtered.pon.vcf.gz``

 # 10. Instalar vcftools

``brew install vcftools``

 *Arquivos do drive*
WP312.filtered.vcf.gz
WP312.filtered.vcf.gz.tbi
 
 # pegando apenas o cabeçalho
``zgrep "\#" WP312.filtered.vcf.gz > header.txt``

``zgrep -v "\#" WP312.filtered.vcf.gz | awk '{print("chr"$0)}' > variants.txt``

``cat header.txt variants.txt > WP312.filtered.chr.vcf``

 *Rodar o vcf compare*

 ``vcf-compare WP312.filtered.pon.vcf.gz ../WP312.filtered.chr.vcf.gz``







