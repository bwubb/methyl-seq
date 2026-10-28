#BSUB -P BISCUIT
#BSUB -R "rusage[mem=8000]"
#BSUB -n 25
#BSUB -q standard
#BSUB -R "span[hosts=1]"
#set -ue pipefail

THR=25

module purge
module load biscuit
module load samtools/1.7
module load bbmap
module load dupsifter
ADAPTERS=/hpcf/apps/bbmap/install/37.28/resources/adapters.fa

### FASTA
FASTA=GRCh38.primary_assembly.genome.pUC19.Lambda.CarT.HHV6B.fa
FILEMAP=../Metadata/cfDNA_filemap.txt

# Find FASTQs
SAMP=%s
FASTQS=$(grep -P "${SAMP}\t" $FILEMAP | awk '{print "../Raw/"$2}')

# Check if R3 present: UMI
printf "\n\n\n"
echo "SAMPLE####################"
echo $SAMP
R1=$( ls $FASTQS | grep "_R1" )
R2=$(  ls $FASTQS | grep "_R3" )
if (( ${#R2} == 0 ));
then
        R2=$(  ls $FASTQS | grep "_R2" )
fi

echo "READ1####"
ls -lh $R1
echo "READ2####"
ls -lh $R2

# Check if BAM exists already
BAM="../Processed/BAM/${SAMP}.bam"
SAM="../Processed/BAM/${SAMP}.sam"
SORTEDBAM="../Processed/BAM/${SAMP}_sorted.bam"

# Merge fastqs
MR1="${TMPDIR}/Read1_${SAMP}.fastq"
MR2="${TMPDIR}/Read2_${SAMP}.fastq"
cat $R1 | gunzip > $MR1
cat $R2 | gunzip > $MR2

# Trim fastq of adaptors
TR1="${TMPDIR}/Read1_Trim_${SAMP}.fastq"
TR2="${TMPDIR}/Read2_Trim_${SAMP}.fastq"

bbduk.sh -Xmx32g threads=$THR \
in1=$MR1 in2=$MR2 \
out1=$TR1 out2=$TR2 \
ref=$ADAPTERS ktrim=r ordered \
k=23 mink=11 hdist=1 tbo tpe \
rcomp=f \
ow=t \
minlen=30

# Map with BWA meth
DUP="../Processed/BAM/${SAMP}dup.bam"
biscuit align -@ $THR -R "@RG\tID:${SAMP}${i}\tSM:${SAMP}\tPL:ILLUMINA\tLB:LIB" $FASTA $TR1 $TR2 | dupsifter $FASTA | samtools sort -@ $THR -T $TMPDIR -o $DUP -O BAM -

samtools index $DUP
biscuit pileup -@ $THR -o "../Processed/Mutations/${SAMP}_pileup.vcf" $FASTA $DUP
bgzip -@ $THR "../Processed/Mutations/${SAMP}_pileup.vcf"
tabix -p vcf "../Processed/Mutations/${SAMP}_pileup.vcf.gz"

biscuit vcf2bed -t cg -k 1 -c "../Processed/Mutations/${SAMP}_pileup.vcf.gz" > "../Processed/Mutations/${SAMP}_methyl_data.bed"
biscuit mergecg $FASTA "../Processed/Mutations/${SAMP}_methyl_data.bed"
bgzip "../Processed/Mutations/${SAMP}_methyl_data.bed"
tabix -p bed "../Processed/Mutations/${SAMP}_methyl_data.bed.gz"

biscuit vcf2bed -t snp -k 3 "../Processed/Mutations/${SAMP}_pileup.vcf.gz" > "../Processed/Mutations/${SAMP}_snp_data.bed"
bgzip "../Processed/Mutations/${SAMP}_snp_data.bed"
tabix -p bed "../Processed/Mutations/${SAMP}_snp_data.bed.gz"

# Call methylation
module load methyldackel/0.5.1
MethylDackel extract -o "../Processed/Betas/${SAMP}" --mergeContext --minDepth 1 --OT 10,140,10,140 --OB 10,140,10,140 $FASTA $DUP






