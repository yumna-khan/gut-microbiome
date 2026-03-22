DB=~/scratch/kraken2_db
INPUT_DIR=~/scratch/asg3/trimmed_reads
OUTPUT_DIR=~/scratch/asg3/kraken2_output

mkdir -p ${OUTPUT_DIR}

for R1 in ${INPUT_DIR}/*_1_trimmed.fastq; do

    SAMPLE=$(basename ${R1} _1_trimmed.fastq)

    kraken2 \
        --db ${DB} \
        --paired \
        --threads 4 \
        --confidence 0.15 \
        --output ${OUTPUT_DIR}/${SAMPLE}.kraken \
        --report ${OUTPUT_DIR}/${SAMPLE}.report \
        ${INPUT_DIR}/${SAMPLE}_1_trimmed.fastq \
        ${INPUT_DIR}/${SAMPLE}_2_trimmed.fastq

done

echo "All samples complete!"
